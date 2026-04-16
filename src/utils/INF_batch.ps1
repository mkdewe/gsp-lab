# INF_batch_fixed.ps1
$BaseDir = "D:\Studia\Projekty\gsp-lab"
$RNAAssessmentDir = "D:\Studia\Projekty\RNA_assessment"
$DockerImage = "mcannotate-runtime"

Write-Host "RNA Puzzles - Batch INF Calculation" -ForegroundColor Cyan
Write-Host "=======================================" -ForegroundColor Cyan

# Utwórz katalog wynikowy
$ResultsDir = "$BaseDir\results"
if (!(Test-Path $ResultsDir)) {
    New-Item -ItemType Directory -Path $ResultsDir -Force
}

# Sprawdź Docker
try {
    docker --version | Out-Null
    Write-Host "✓ Docker is available" -ForegroundColor Green
} catch {
    Write-Host "✗ Docker is not available" -ForegroundColor Red
    exit 1
}

# Sprawdź obraz Docker
$imageExists = docker images -q $DockerImage
if (!$imageExists) {
    Write-Host "✗ Docker image $DockerImage not found" -ForegroundColor Red
    exit 1
}
Write-Host "✓ Docker image $DockerImage found" -ForegroundColor Green

Write-Host "`nRunning batch processing..." -ForegroundColor Yellow

# Użyj heredoc dla kodu Python
$pythonCode = @'
import sys
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/work/src/utils")

from INF import compute_INF
import json
import csv
import os
from pathlib import Path

results = []
count = 0

data_dir = Path("/work/data/finished")
for puzzle_dir in data_dir.iterdir():
    if puzzle_dir.is_dir():
        puzzle_name = puzzle_dir.name
        pdb_dir = puzzle_dir / "pdb"
        
        if pdb_dir.exists():
            # Znajdź pliki solution
            solutions = list(pdb_dir.glob("*solution*.pdb"))
            if not solutions:
                continue
            
            # Znajdź pliki modeli
            models = []
            for pdb in pdb_dir.glob("*.pdb"):
                if "solution" not in pdb.name.lower():
                    models.append(pdb)
            
            if models:
                solution = solutions[0]
                
                for model in models:
                    count += 1
                    try:
                        print(f"[{count}] {puzzle_name}: {model.name}")
                        result = compute_INF(str(solution), str(model))
                        
                        result["puzzle"] = puzzle_name
                        result["model"] = model.name
                        result["solution"] = solution.name
                        
                        results.append(result)
                        
                        inf_val = result.get("inf_all", "N/A")
                        if inf_val != "N/A":
                            print(f"   INF-all: {inf_val:.4f}")
                        else:
                            print(f"   INF-all: N/A")
                            
                    except Exception as e:
                        print(f"   Error: {str(e)[:60]}")

# Zapisz wyniki
os.makedirs("/work/results", exist_ok=True)

# JSON
json_path = "/work/results/inf_results.json"
with open(json_path, "w") as f:
    json.dump(results, f, indent=2)

# CSV
if results:
    csv_path = "/work/results/inf_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["puzzle", "model", "solution", "rmsd", "inf_all", "inf_wc", "inf_nwc", "inf_stack"])
        for r in results:
            writer.writerow([
                r["puzzle"], r["model"], r["solution"],
                r.get("rmsd"), r.get("inf_all"), r.get("inf_wc"),
                r.get("inf_nwc"), r.get("inf_stack")
            ])

print(f"\n========== SUMMARY ==========")
print(f"Total processed: {count}")
print(f"Successful: {len(results)}")
print(f"Failed: {count - len(results)}")
print(f"Results saved to: /work/results/")
'@

# Zapisz kod Python do tymczasowego pliku
$tempPythonFile = [System.IO.Path]::GetTempFileName() + ".py"
$pythonCode | Out-File -FilePath $tempPythonFile -Encoding UTF8

# Polecenie Docker
$dockerCmd = @"
docker run --rm `
  -v "${BaseDir}:/work" `
  -v "${RNAAssessmentDir}:/rna_assessment" `
  -e RNA_ASSESSMENT_ROOT=/rna_assessment `
  ${DockerImage} `
  bash -c "
    ln -sf /opt/mcannotate/MC-Annotate /usr/local/bin/MC-Annotate 2>/dev/null || true
    ln -sf /opt/mcannotate/MC-Annotate /work/MC-Annotate 2>/dev/null || true
    ln -sf /opt/mcannotate/MC-Annotate /work/src/utils/MC-Annotate 2>/dev/null || true
    cd /work && python3 $(Get-Content $tempPythonFile | Out-String | ForEach-Object { $_.Replace('"', '\"').Replace('$', '\$') })
  "
"@

Invoke-Expression $dockerCmd

# Usuń tymczasowy plik
Remove-Item $tempPythonFile -Force -ErrorAction SilentlyContinue

Write-Host "`nResults:" -ForegroundColor Green
if (Test-Path "$ResultsDir\inf_results.json") {
    Write-Host "✓ JSON results: $ResultsDir\inf_results.json" -ForegroundColor Green
}
if (Test-Path "$ResultsDir\inf_results.csv") {
    Write-Host "✓ CSV results: $ResultsDir\inf_results.csv" -ForegroundColor Green
}

Write-Host "`n=======================================" -ForegroundColor Cyan
Write-Host "Batch processing completed!" -ForegroundColor Cyan