# run_INF_finished_to_INF_all.ps1
# PowerShell 7 — jeden plik, bez zewnętrznych zależności
#
# Logika doboru targetu (identyczna jak w oryginalnym skrypcie z marca):
#
# 1. Dla plików *_refined.pdb:
#    - szukaj predefiniowanego targetu: <base>_solution.pdb (np. 1_santalucia_1_solution.pdb)
#    - jeśli istnieje → użyj TYLKO jego
#    - jeśli nie → użyj globalnych solution(ów) i wybierz najlepszy
#
# 2. Dla pozostałych modeli (nie refined, nie solution):
#    - użyj globalnych solution(ów) i wybierz najlepszy
#
# Globalne solutions = pliki pasujące do wzorca *_solution_0.pdb lub solution_0.pdb
# (ogólne, nienazwane targety)

$ROOT   = "D:\Studia\Projekty\gsp-lab"
$DATA   = Join-Path $ROOT "data\finished"
$OUTDIR = Join-Path $ROOT "results\INF_all"

if (!(Test-Path $OUTDIR)) { New-Item -ItemType Directory -Path $OUTDIR | Out-Null }
$LOGDIR = Join-Path $OUTDIR "logs"
if (!(Test-Path $LOGDIR)) { New-Item -ItemType Directory -Path $LOGDIR | Out-Null }

Write-Host "START INF — FINISHED DATASET"
Write-Host "DATA DIR   = $DATA"
Write-Host "OUTPUT DIR = $OUTDIR"

$cases = Get-ChildItem $DATA -Directory | Select-Object -ExpandProperty Name
Write-Host ("Found {0} case(s)" -f $cases.Count)

$cases | ForEach-Object -Parallel {

    $case    = $_
    $caseLog = Join-Path $using:LOGDIR "$case.log"

    function Log($msg) {
        "$((Get-Date).ToString('yyyy-MM-dd HH:mm:ss')) `t $msg" | Add-Content -Path $caseLog -Encoding utf8
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

    # Uruchom Docker dla jednej pary (targetFullPath, modelFullPath)
    # Zwraca obiekt z inf_all/inf_wc/inf_nwc lub $null
    function Run-INF($targetFullPath, $modelFullPath) {
        $relTarget = $targetFullPath.Substring($using:ROOT.Length + 1).Replace('\', '/')
        $relModel  = $modelFullPath.Substring($using:ROOT.Length + 1).Replace('\', '/')

        $dockerArgs = @(
            'run', '--rm',
            '-v', ("{0}:/work" -f $using:ROOT),
            '-v', 'D:\Studia\Projekty\RNA_assessment:/rna_assessment',
            '-e', 'RNA_ASSESSMENT_ROOT=/rna_assessment',
            '-e', 'PYTHONPATH=/rna_assessment',
            'mcannotate-runtime',
            'python3', '/work/src/utils/INF.py',
            "/work/$relTarget",
            "/work/$relModel"
        )

        Log ("Docker: {0}" -f ($dockerArgs -join ' '))
        try   { $output = & docker @dockerArgs 2>&1; $exit = $LASTEXITCODE }
        catch { $output = $_.Exception.Message; $exit = 1 }

        Log ("exit={0}" -f $exit)
        if ($exit -ne 0) { Log ("stderr: {0}" -f ([string]::Join(' ', $output))); return $null }

        return Parse-JsonFromOutput $output
    }

    # Spośród listy targetów wybierz ten z najwyższym INF_all
    # Zwraca (bestResult, bestTargetFile) lub ($null, $null)
    function Best-INF($targetFiles, $modelFullPath) {
        $bestVal    = -1.0
        $bestResult = $null
        $bestTarget = $null

        foreach ($tgt in $targetFiles) {
            Log ("  testing target: {0}" -f $tgt.Name)
            $res = Run-INF $tgt.FullName $modelFullPath
            if ($null -eq $res) { Log "  -> no result"; continue }

            $val = try {
                [double]::Parse((ToInvariant $res.inf_all), [System.Globalization.CultureInfo]::InvariantCulture)
            } catch { -1.0 }

            Log ("  -> INF_all={0}" -f $val)
            if ($val -gt $bestVal) {
                $bestVal    = $val
                $bestResult = $res
                $bestTarget = $tgt
            }
        }
        return $bestResult, $bestTarget
    }

    function AlreadyRecorded($modelFileName, $summaryFile) {
        if (!(Test-Path $summaryFile)) { return $false }
        $pattern = ",{0}," -f [regex]::Escape($modelFileName)
        return Select-String -Path $summaryFile -Pattern $pattern -Quiet
    }

    function Write-Result($summaryFile, $caseId, $modelFile, $targetFile, $res) {
        if (!(Test-Path $summaryFile)) {
            "case_id,model_name,target_name,INF_all,INF_WC,INF_nWC" | Out-File $summaryFile -Encoding utf8
        }
        $line = "{0},{1},{2},{3},{4},{5}" -f `
            $caseId, $modelFile.Name, $targetFile.Name,
            (ToInvariant $res.inf_all),
            (ToInvariant $res.inf_wc),
            (ToInvariant $res.inf_nwc)
        $line | Out-File -Append $summaryFile -Encoding utf8
        Log ("Saved: {0}" -f $line)
        Write-Host ("  Saved: {0} | target={1} | INF_all={2}" -f $modelFile.Name, $targetFile.Name, (ToInvariant $res.inf_all))
    }

    # ------------------------------------------------------------------ #

    Log "START case: $case"

    $pdbDir      = Join-Path $using:DATA "$case\pdb"
    $summaryFile = Join-Path $using:OUTDIR "$case.csv"

    if (!(Test-Path $pdbDir)) {
        Log "SKIP (no pdb dir)"; Write-Host "SKIP (no pdb): $case"; return
    }

    $allPdbs = Get-ChildItem $pdbDir -File -Filter "*.pdb"

    # Wszystkie solution files
    $allSolutions = $allPdbs | Where-Object { $_.Name -match "solution" }

    # Globalne solutions: pasują do *_solution_0.pdb lub solution_0.pdb (ogólny target, nie powiązany z konkretnym labem)
    $globalSolutions = $allSolutions | Where-Object { $_.Name -match "_solution_0\.pdb$" -or $_.Name -match "^solution_0\.pdb$" }

    # Modele: wszystko co nie jest solution
    $models = $allPdbs | Where-Object { $_.Name -notmatch "solution" }

    if ($allSolutions.Count -eq 0) { Log "SKIP (no solutions)"; Write-Host "SKIP (no solutions): $case"; return }
    if ($models.Count -eq 0)       { Log "SKIP (no models)";    Write-Host "SKIP (no models): $case";    return }

    Log ("Global solutions ({0}): {1}" -f $globalSolutions.Count, ($globalSolutions.Name -join ', '))
    Log ("All solutions    ({0}): {1}" -f $allSolutions.Count,    ($allSolutions.Name    -join ', '))
    Log ("Models           ({0})"      -f $models.Count)
    Write-Host ("Processing: {0}  [{1} models, {2} solutions ({3} global)]" -f $case, $models.Count, $allSolutions.Count, $globalSolutions.Count)

    foreach ($mod in $models) {
        if (AlreadyRecorded $mod.Name $summaryFile) {
            Log ("SKIP already recorded: {0}" -f $mod.Name); continue
        }

        Log ("RUN model: {0}" -f $mod.Name)

        # Czy to plik refined?
        if ($mod.Name -match "^(.+)_refined\.pdb$") {
            $base = $matches[1]
            # Szukaj predefiniowanego targetu: <base>_solution.pdb
            $specificTarget = $allSolutions | Where-Object { $_.Name -eq "$base`_solution.pdb" } | Select-Object -First 1

            if ($null -ne $specificTarget) {
                # Predefiniowany target — używamy TYLKO jego, bez best-of
                Log ("  refined with specific target: {0}" -f $specificTarget.Name)
                $res = Run-INF $specificTarget.FullName $mod.FullName
                if ($null -ne $res) {
                    Write-Result $summaryFile $case $mod $specificTarget $res
                } else {
                    Log ("  no result for specific target {0}" -f $specificTarget.Name)
                    Write-Host ("  WARN: no result for {0}" -f $mod.Name)
                }
                continue
            }
            # Brak predefiniowanego targetu → fallback do globalnych (best-of)
            Log ("  refined but no specific target found, using global solutions")
        }

        # Modele bez predefiniowanego targetu → best-of spośród globalnych solutions
        if ($globalSolutions.Count -eq 0) {
            Log ("  no global solutions available, skipping {0}" -f $mod.Name)
            Write-Host ("  WARN: no global solutions for {0}" -f $mod.Name)
            continue
        }

        $bestResult, $bestTarget = Best-INF $globalSolutions $mod.FullName

        if ($null -ne $bestResult) {
            Write-Result $summaryFile $case $mod $bestTarget $bestResult
        } else {
            Log ("  no result for any global solution for {0}" -f $mod.Name)
            Write-Host ("  WARN: no result for {0}" -f $mod.Name)
        }
    }

    Log "DONE case: $case"

} -ThrottleLimit 13

Write-Host "DONE"