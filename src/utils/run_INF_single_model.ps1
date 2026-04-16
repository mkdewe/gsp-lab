<#
.SYNOPSIS
Oblicza INF dla jednego modelu względem wszystkich plików *solution*.pdb z tego samego katalogu pdb/.
Zapisuje tylko najlepszy wynik (najwyższy INF_all) - analogicznie do select_best_for_group w gsp_compute.py.

.PARAMETER ModelPath
Ścieżka do pliku .pdb modelu (wymagana).

.PARAMETER OutDir
Katalog wynikowy. Domyślnie: $ROOT\results\INF_all\single_model.

.PARAMETER BestOnly
Jeśli $true (domyślnie), zapisuje tylko najlepszy wynik. Jeśli $false, zapisuje wszystkie wyniki.

.PARAMETER SummaryFile
Ścieżka do pliku CSV z wynikami. Domyślnie: $OutDir\best_results.csv.

.EXAMPLE
.\run_INF_single_model.ps1 -ModelPath "D:\Studia\Projekty\gsp-lab\data\finished\CR1149\pdb\CR1149_TS235_3_refined.pdb"
#>

param(
    [Parameter(Mandatory=$true)]
    [string]$ModelPath,
    [string]$OutDir,
    [bool]$BestOnly = $true,
    [string]$SummaryFile
)

$ModelPath = $ModelPath -replace '^"|"$', ''

if (-not (Test-Path $ModelPath)) {
    Write-Error "Plik modelu nie istnieje: $ModelPath"
    exit 1
}

# Ustal ROOT (gsp-lab)
$ROOT = $ModelPath
while ($ROOT -notlike '*\gsp-lab' -and $ROOT -ne '') {
    $ROOT = Split-Path $ROOT -Parent
}
if ($ROOT -eq '') { Write-Error "Nie można odnaleźć katalogu gsp-lab w ścieżce."; exit 1 }
$ROOT = $ROOT.TrimEnd('\')

$pdbDir        = Split-Path $ModelPath -Parent
$caseName      = Split-Path (Split-Path $pdbDir -Parent) -Leaf
$modelFileName = Split-Path $ModelPath -Leaf

if (-not $OutDir) { $OutDir = Join-Path $ROOT "results\INF_all\single_model" }
if (-not (Test-Path $OutDir)) { New-Item -ItemType Directory -Path $OutDir -Force | Out-Null }
$logDir = Join-Path $OutDir "logs"
if (-not (Test-Path $logDir)) { New-Item -ItemType Directory -Path $logDir -Force | Out-Null }

if ($BestOnly) {
    if (-not $SummaryFile) { $SummaryFile = Join-Path $OutDir "best_results.csv" }
    $csvFile = $SummaryFile
} else {
    $csvFile = Join-Path $OutDir "results.csv"
}
$logFile       = Join-Path $logDir "run.log"
$dockerOutFile = Join-Path $logDir "docker_output.log"

function Log($msg) {
    "$((Get-Date).ToString('yyyy-MM-dd HH:mm:ss')) `t $msg" | Add-Content -Path $logFile -Encoding utf8
}

function ToInvariant([object] $x) {
    if ($null -eq $x) { return "" }
    if ($x -is [string]) { return $x }
    try { return $x.ToString([System.Globalization.CultureInfo]::InvariantCulture) } catch { return [string]$x }
}

function Parse-JsonFromOutput([object] $outputLines) {
    $joined = if ($outputLines -is [array]) { [string]::Join("`n", $outputLines) } else { [string]$outputLines }
    $start = $joined.IndexOf('{'); $end = $joined.LastIndexOf('}')
    if ($start -ge 0 -and $end -ge $start) {
        $jsonText = $joined.Substring($start, $end - $start + 1)
        try   { return $jsonText | ConvertFrom-Json }
        catch { Log "Parse-JsonFromOutput: ConvertFrom-Json failed"; return $null }
    }
    Log "Parse-JsonFromOutput: no JSON block found"
    return $null
}

# --- walidacja ---
if (-not $ModelPath.StartsWith($ROOT)) { Write-Error "Ścieżka modelu nie zaczyna się od ROOT ($ROOT)"; exit 1 }

# solutions = wszystkie pliki zawierające "solution" w nazwie (tak samo jak collect_files w gsp_compute.py)
$solutions = Get-ChildItem $pdbDir -File -Filter "*.pdb" | Where-Object { $_.Name -match "solution" }

if ($solutions.Count -eq 0) {
    Write-Host "UWAGA: Brak plików *solution*.pdb w katalogu $pdbDir dla modelu $modelFileName"
    exit 0
}

Log "===== START ====="
Log "Model:       $ModelPath"
Log "pdbDir:      $pdbDir"
Log "Solutions:   $($solutions.Count) -> $($solutions.Name -join ', ')"
Log "BestOnly:    $BestOnly"
Log "CSV file:    $csvFile"

if (-not (Test-Path $csvFile)) {
    if ($BestOnly) {
        "case_id,model_name,best_solution_name,best_INF_all,INF_WC,INF_nWC" | Out-File $csvFile -Encoding utf8
    } else {
        "case_id,model_name,solution_name,INF_all,INF_WC,INF_nWC" | Out-File $csvFile -Encoding utf8
    }
}

$relativeModel        = $ModelPath.Substring($ROOT.Length + 1).Replace('\', '/')
$modelPathInContainer = "/work/$relativeModel"

$bestInfAll      = -1.0
$bestSolutionName = ""
$bestResult      = $null

foreach ($sol in $solutions) {
    Log "Testuję: solution=$($sol.Name)"

    $relativeSol      = $sol.FullName.Substring($ROOT.Length + 1).Replace('\', '/')
    $solPathInContainer = "/work/$relativeSol"

    $dockerArgs = @(
        'run','--rm',
        '-v', ("{0}:/work" -f $ROOT),
        '-v', 'D:\Studia\Projekty\RNA_assessment:/rna_assessment',
        '-e', 'RNA_ASSESSMENT_ROOT=/rna_assessment',
        '-e', 'PYTHONPATH=/rna_assessment',
        'mcannotate-runtime',
        'python3', '/work/src/utils/INF.py',
        $solPathInContainer,
        $modelPathInContainer
    )

    Log ("Docker cmd: docker {0}" -f ($dockerArgs -join ' '))
    try   { $output = & docker @dockerArgs 2>&1; $exit = $LASTEXITCODE }
    catch { $output = $_.Exception.Message; $exit = 1 }

    $output | Out-File -FilePath $dockerOutFile -Encoding utf8 -Append
    Log "Docker exit: $exit"

    if ($exit -ne 0) { Log "Non-zero exit dla solution=$($sol.Name) -> pomijam"; continue }

    $jsonObj = Parse-JsonFromOutput $output
    if ($null -eq $jsonObj) { Log "Brak JSON dla solution=$($sol.Name) -> pomijam"; continue }

    $infAllVal = try {
        [double]::Parse((ToInvariant $jsonObj.inf_all), [System.Globalization.CultureInfo]::InvariantCulture)
    } catch { -1.0 }

    Log ("INF_all={0} dla solution={1}" -f $infAllVal, $sol.Name)
    Write-Host ("  {0} vs {1} -> INF_all = {2}" -f $modelFileName, $sol.Name, $infAllVal)

    if (-not $BestOnly) {
        $line = "{0},{1},{2},{3},{4},{5}" -f $caseName, $modelFileName, $sol.Name,
            (ToInvariant $jsonObj.inf_all), (ToInvariant $jsonObj.inf_wc), (ToInvariant $jsonObj.inf_nwc)
        $line | Out-File -Append $csvFile -Encoding utf8
        Log "Zapisano: $line"
    } else {
        if ($infAllVal -gt $bestInfAll) {
            $bestInfAll       = $infAllVal
            $bestSolutionName = $sol.Name
            $bestResult       = $jsonObj
            Log "Nowy najlepszy: solution=$($sol.Name) INF_all=$infAllVal"
        }
    }
}

if ($BestOnly) {
    if ($bestSolutionName -ne "") {
        $bestLine = "{0},{1},{2},{3},{4},{5}" -f `
            $caseName, $modelFileName, $bestSolutionName,
            (ToInvariant $bestResult.inf_all),
            (ToInvariant $bestResult.inf_wc),
            (ToInvariant $bestResult.inf_nwc)
        $bestLine | Out-File -Append $csvFile -Encoding utf8
        Log "Zapisano najlepszy wynik: $bestLine"
        Write-Host ("Najlepszy dla {0}: {1} -> INF_all = {2}" -f $modelFileName, $bestSolutionName, $bestInfAll)
    } else {
        Log "Brak jakiegokolwiek wyniku INF dla modelu $modelFileName"
        Write-Host "UWAGA: Brak wyników INF dla modelu $modelFileName"
    }
}

Log "===== KONIEC ====="
Write-Host "Wyniki: $csvFile"
