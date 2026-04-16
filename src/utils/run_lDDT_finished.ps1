<#
.SYNOPSIS
Oblicza lDDT dla jednego modelu względem wszystkich plików referencyjnych (_solution.pdb) z katalogu references.
Domyślnie zapisuje tylko najlepszy wynik (najwyższy lDDT) dla modelu.

.DESCRIPTION
Skrypt uruchamia kontener openstructure, dla każdej referencji wykonuje compare-structures --lddt -rna,
odczytuje wynik z pliku JSON i zapisuje do pliku CSV. Domyślnie zapisuje tylko najlepszy wynik (najwyższy lDDT).
Można też zachować wszystkie wyniki (parametr -BestOnly:$false) oraz zachować pliki JSON.

.PARAMETER ModelPath
Ścieżka do pliku .pdb modelu (wymagana).

.PARAMETER ReferencesDir
Ścieżka do katalogu z referencjami. Jeśli nie podano, zakłada katalog "references" obok "models".

.PARAMETER OutDir
Katalog wynikowy. Domyślnie: $ROOT\results\lDDT\single_model.

.PARAMETER KeepJson
Jeśli podany, nie usuwa tymczasowych plików JSON z wynikami.

.PARAMETER BestOnly
Jeśli $true (domyślnie), zapisuje tylko najlepszy wynik (najwyższy lDDT) dla modelu. Jeśli $false, zapisuje wszystkie wyniki.

.PARAMETER SummaryFile
Ścieżka do pliku CSV z podsumowaniem najlepszych wyników. Domyślnie: $OutDir\best_results.csv.
Używany tylko gdy -BestOnly jest włączony.

.EXAMPLE
.\run_lddt_single_model_vs_references.ps1 -ModelPath "D:\Studia\Projekty\gsp-lab\data\finished\CR1149\pdb\CR1149_TS235_3_refined.pdb"

.EXAMPLE
.\run_lddt_single_model_vs_references.ps1 -ModelPath "..." -KeepJson -BestOnly:$false
#>

param(
    [Parameter(Mandatory=$true)]
    [string]$ModelPath,
    [string]$ReferencesDir,
    [string]$OutDir,
    [switch]$KeepJson,
    [bool]$BestOnly = $true,   # zmieniono z [switch] na [bool] – teraz można ustawić domyślnie
    [string]$SummaryFile
)

# ----- Usuń ewentualne cudzysłowy z początku/końca (gdy podano interaktywnie) -----
$ModelPath = $ModelPath -replace '^"|"$', ''

# ----- Sprawdź, czy plik modelu istnieje -----
if (-not (Test-Path $ModelPath)) {
    Write-Error "Plik modelu nie istnieje: $ModelPath"
    exit 1
}

# ----- Ustal ścieżkę główną projektu (ROOT) -----
$ROOT = $ModelPath
while ($ROOT -notlike '*\gsp-lab' -and $ROOT -ne '') {
    $ROOT = Split-Path $ROOT -Parent
}
if ($ROOT -eq '') {
    Write-Error "Nie można odnaleźć katalogu głównego projektu (gsp-lab) w ścieżce."
    exit 1
}
$ROOT = $ROOT.TrimEnd('\')

# ----- Domyślne wartości -----
$caseName = Split-Path (Split-Path $ModelPath -Parent) -Parent | Split-Path -Leaf
$modelFileName = Split-Path $ModelPath -Leaf

if (-not $ReferencesDir) {
    $caseDir = Split-Path (Split-Path $ModelPath -Parent) -Parent
    $ReferencesDir = Join-Path $caseDir "references"
}

if (-not $OutDir) {
    $OutDir = Join-Path $ROOT "results\lDDT\single_model"
}

# ----- Utwórz katalogi wynikowe -----
if (-not (Test-Path $OutDir)) {
    New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
}
$logDir = Join-Path $OutDir "logs"
if (-not (Test-Path $logDir)) {
    New-Item -ItemType Directory -Path $logDir -Force | Out-Null
}

# ----- Usuń stare pliki tymczasowe -----
Get-ChildItem $OutDir -Filter "temp_*.json" | Remove-Item -Force

# ----- Pliki wyjściowe -----
if ($BestOnly) {
    if (-not $SummaryFile) {
        $SummaryFile = Join-Path $OutDir "best_results.csv"
    }
    $csvFile = $SummaryFile
} else {
    $csvFile = Join-Path $OutDir "results.csv"
}
$logFile = Join-Path $logDir "run.log"
$dockerOutFile = Join-Path $logDir "docker_output.log"

# ----- Funkcja logująca -----
function Log($msg) {
    $time = (Get-Date).ToString("yyyy-MM-dd HH:mm:ss")
    "$time `t $msg" | Add-Content -Path $logFile -Encoding utf8
}

# ----- Funkcja parsująca lDDT z pliku JSON -----
function Get-LddtFromJsonFile($jsonPath) {
    if (-not (Test-Path $jsonPath)) { return $null }
    try {
        $content = Get-Content $jsonPath -Raw -Encoding utf8
        $obj = $content | ConvertFrom-Json
        return $obj.lddt
    } catch {
        Log "Błąd parsowania JSON z pliku $jsonPath : $_"
        return $null
    }
}

# ----- Funkcja pomocnicza do zapisu CSV (z niezmienną kulturą) -----
function ToInvariant($x) {
    if ($null -eq $x) { return "" }
    if ($x -is [string]) { return $x }
    try { return $x.ToString([System.Globalization.CultureInfo]::InvariantCulture) } catch { return [string]$x }
}

# ----- Sprawdź, czy ścieżka modelu zaczyna się od ROOT -----
if (-not $ModelPath.StartsWith($ROOT)) {
    Log "BŁĄD: Ścieżka modelu nie zaczyna się od ROOT ($ROOT)"
    exit 1
}

# ----- Główna pętla -----
Log "===== START ====="
Log "Model: $ModelPath"
Log "References dir: $ReferencesDir"
Log "Output dir: $OutDir"
Log "ROOT: $ROOT"
Log "Keep JSON: $KeepJson"
Log "BestOnly: $BestOnly"
Log "CSV file: $csvFile"

# Sprawdź, czy katalog referencji istnieje
if (-not (Test-Path $ReferencesDir)) {
    Log "BŁĄD: Katalog referencji nie istnieje: $ReferencesDir"
    exit 1
}

# Pobierz wszystkie pliki referencyjne
$refFiles = Get-ChildItem $ReferencesDir -Filter "*_solution.pdb" -File
if ($refFiles.Count -eq 0) {
    Log "BŁĄD: Nie znaleziono plików *_solution.pdb w katalogu referencji."
    exit 1
}
Log "Znaleziono $($refFiles.Count) plików referencyjnych."

# Inicjalizuj plik CSV (nagłówek) tylko jeśli nie istnieje
if (-not (Test-Path $csvFile)) {
    if ($BestOnly) {
        "case_id,model_name,best_ref_name,best_lddt" | Out-File $csvFile -Encoding utf8
    } else {
        "case_id,model_name,ref_name,lddt" | Out-File $csvFile -Encoding utf8
    }
}

# Zmienne do przechowywania najlepszego wyniku
$bestLddt = -1.0
$bestRefName = ""

$volumeMount = ($ROOT -replace '\\', '/') + ':/data'

foreach ($ref in $refFiles) {
    $refName = $ref.Name

    Log "Przetwarzanie: $refName"

    # Generuj losowy identyfikator dla pliku tymczasowego
    $rand = Get-Random
    $tempJsonHost = Join-Path $OutDir "temp_$($caseName)_$rand.json"
    $tempJsonContainer = "/data/results/lDDT/single_model/temp_$($caseName)_$rand.json"

    # ----- Buduj ścieżki względne w kontenerze -----
    $relativeModel = $ModelPath.Substring($ROOT.Length + 1).Replace('\', '/')
    $modelPathInContainer = "/data/$relativeModel"

    if (-not $ref.FullName.StartsWith($ROOT)) {
        Log "BŁĄD: Referencja nie zaczyna się od ROOT: $($ref.FullName)"
        continue
    }
    $relativeRef = $ref.FullName.Substring($ROOT.Length + 1).Replace('\', '/')
    $refPathInContainer = "/data/$relativeRef"

    $dockerArgs = @(
        'run', '--rm',
        '-v', $volumeMount,
        'registry.scicore.unibas.ch/schwede/openstructure:latest',
        'compare-structures',
        '-r', $refPathInContainer,
        '-m', $modelPathInContainer,
        '--lddt',
        '--lddt-no-stereochecks',
        '-o', $tempJsonContainer
    )

    Log "Docker cmd: docker $($dockerArgs -join ' ')"
    try {
        $output = & docker @dockerArgs 2>&1
        $exit = $LASTEXITCODE
    } catch {
        $output = $_.Exception.Message
        $exit = 1
    }

    $output | Out-File -FilePath $dockerOutFile -Encoding utf8 -Append
    Log "Docker exit code: $exit"

    # Odczytaj wynik z pliku JSON
    if (Test-Path $tempJsonHost) {
        $lddt = Get-LddtFromJsonFile $tempJsonHost
        if ($null -ne $lddt) {
            $lddt_s = ToInvariant $lddt
            if (-not $BestOnly) {
                $csvLine = "$caseName,$modelFileName,$refName,$lddt_s"
                $csvLine | Out-File -Append $csvFile -Encoding utf8
                Log "Zapisano: $csvLine"
                Write-Host "OK: $caseName - $refName -> lDDT = $lddt"
            } else {
                if ($lddt -gt $bestLddt) {
                    $bestLddt = $lddt
                    $bestRefName = $refName
                    Log "Nowy najlepszy wynik: $refName -> $lddt (poprzedni: $bestLddt)"
                }
                Write-Host "OK: $caseName - $refName -> lDDT = $lddt"
            }
        } else {
            Log "Plik JSON istnieje, ale nie zawiera pola lddt."
        }
        if (-not $KeepJson) {
            Remove-Item $tempJsonHost -Force
        } else {
            Log "JSON zachowany: $tempJsonHost"
        }
    } else {
        Log "Plik JSON nie został utworzony: $tempJsonHost"
    }
}

# Jeśli tryb BestOnly i znaleziono co najmniej jeden wynik
if ($BestOnly -and $bestRefName -ne "") {
    $bestLine = "$caseName,$modelFileName,$bestRefName,$(ToInvariant $bestLddt)"
    $bestLine | Out-File -Append $csvFile -Encoding utf8
    Log "Zapisano najlepszy wynik: $bestLine"
    Write-Host "Najlepszy wynik dla ${modelFileName}: $bestRefName -> lDDT = $bestLddt"
} elseif ($BestOnly -and $bestRefName -eq "") {
    Log "Nie udało się uzyskać żadnego wyniku lDDT."
    Write-Host "UWAGA: Brak wyników lDDT dla modelu $modelFileName."
}

Log "===== KONIEC ====="
Write-Host "Wyniki zapisano w: $csvFile"
if ($KeepJson) {
    Write-Host "Pliki JSON zachowane w katalogu: $OutDir"
}