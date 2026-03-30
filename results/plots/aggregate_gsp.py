import csv
import os
from pathlib import Path
from collections import defaultdict

# --- KONFIGURACJA ŚCIEŻEK ---
BASE_DIR = Path("D:/Studia/Projekty/gsp-lab/results/gsp/finished")
TYPES_CSV = Path("D:/Studia/Projekty/gsp-lab/results/plots/RNA-types.csv")

# --- WCZYTAJ MAPOWANIE target -> type ---
type_map = {}
if TYPES_CSV.exists():
    with open(TYPES_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            target = row.get('target_label', '').strip()
            typ = row.get('type', '').strip()
            if target and typ:
                type_map[target] = typ
else:
    print(f"Ostrzeżenie: nie znaleziono pliku {TYPES_CSV}")

# --- STRUKTURA DANYCH: klucz = próg, wartość = lista wierszy ---
threshold_data = defaultdict(list)

# --- PRZETWARZANIE WSZYSTKICH PODKATALOGÓW ---
for target_dir in BASE_DIR.iterdir():
    if not target_dir.is_dir():
        continue
    target_name = target_dir.name

    # Pomiń katalog 'gGSP' (nie zawiera report.txt)
    if target_name == "gGSP":
        continue

    report_path = target_dir / "report.txt"
    if not report_path.exists():
        print(f"Brak report.txt w {target_dir}")
        continue

    # --- Ustal typ i źródło ---
    typ = type_map.get(target_name, "Unknown")
    if target_name.startswith("CR"):
        source = "CASP RNA"
    elif target_name.startswith("PZ"):
        source = "RNA Puzzles"
    else:
        source = "Unknown"

    # --- ODCZYTAJ ZAWARTOŚĆ PLIKU ---
    with open(report_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # --- ZNAJDŹ LINIĘ NAGŁÓWKA (zaczyna się od "ModelGroup") ---
    header_idx = None
    header_line = None
    for i, line in enumerate(lines):
        if line.strip().startswith("ModelGroup"):
            header_idx = i
            header_line = line
            break
    if header_idx is None:
        print(f"Brak nagłówka w {report_path}")
        continue

    # --- WYZNACZ POZYCJE KOLUMN (fixed‑width) ---
    # Nagłówek ma stałe szerokości, które są takie same dla całego pliku.
    # Szukamy początku każdej kolumny po spacji.
    # Uproszczenie: pierwsze trzy kolumny (ModelGroup, BestModel, BestSolution)
    # są oddzielone dużą liczbą spacji – znajdziemy ich granice.
    # Kolumny z progami są wyrównane do prawej, ale ich pozycje można odczytać
    # z nagłówka (np. "2.0" zaczyna się w kolumnie 150).
    # Dla bezpieczeństwa użyjemy indeksów opartych na długości nagłówka.

    # Dzielimy nagłówek na kolumny po spacji, aby poznać nazwy.
    # Zrobimy to tylko dla nazw, ale potem użyjemy ich do określenia pozycji.
    header_parts = header_line.split()
    # Oczekiwane nazwy: ModelGroup, BestModel, BestSolution, 2.0, 4.0, 5.0, 8.0, ...
    # W oryginalnym przykładzie jest ich 7.
    # Znajdźmy indeksy pierwszego i ostatniego wystąpienia każdej nazwy.
    # Ustalmy granice kolumn jako początek każdej nazwy w nagłówku.
    col_starts = []
    col_names = []
    pos = 0
    for part in header_parts:
        # Znajdź pozycję part w nagłówku od aktualnej pozycji
        start = header_line.find(part, pos)
        if start == -1:
            # może być wielokrotnie – szukaj od początku
            start = header_line.find(part)
        if start != -1:
            col_starts.append(start)
            col_names.append(part)
            pos = start + len(part)
    # Dodaj koniec linii jako koniec ostatniej kolumny
    col_starts.append(len(header_line))
    col_names.append('')  # dummy

    # Mapowanie nazwy progu na indeks kolumny (w liście col_starts)
    threshold_cols = {}
    for i, name in enumerate(col_names[:-1]):  # pomijamy dummy
        # Jeśli nazwa jest liczbą (próg)
        try:
            float(name)
            threshold_cols[name] = i
        except ValueError:
            # nie jest liczbą – to kolumny z tekstem (ModelGroup, BestModel, BestSolution)
            pass

    # --- PRZESUŃ WSKAŹNIK DO POCZĄTKU DANYCH ---
    # Po nagłówku jest linia z myślnikami
    data_start = header_idx + 1
    # Pomijaj linie składające się tylko z myślników
    while data_start < len(lines) and lines[data_start].strip().startswith('-'):
        data_start += 1

    # --- PRZETWARZAJ LINIE Z DANYMI ---
    for line in lines[data_start:]:
        line_stripped = line.rstrip('\n')
        if not line_stripped:
            continue
        # Zatrzymaj się, gdy napotkamy sekcję "=== FAILED MODELS ===" lub "Failed:"
        if line_stripped.startswith('=== FAILED MODELS ===') or line_stripped.startswith('Failed:'):
            break

        # Wyciągnij ModelGroup z pierwszej kolumny (używając pozycji)
        if len(col_starts) < 2:
            print(f"Błąd: nieprawidłowa liczba kolumn w {report_path}")
            continue
        start0 = col_starts[0]
        end0 = col_starts[1]
        if end0 > len(line_stripped):
            # Linia może być krótsza niż oczekiwana szerokość kolumny – spróbujmy do końca
            end0 = len(line_stripped)
        model_group = line_stripped[start0:end0].strip()

        # Dla każdego progu wyciągnij wartość
        for thresh, col_idx in threshold_cols.items():
            start = col_starts[col_idx]
            end = col_starts[col_idx+1]
            if end > len(line_stripped):
                end = len(line_stripped)
            val_str = line_stripped[start:end].strip()
            # Spróbuj przekonwertować na float
            try:
                score = float(val_str)
            except ValueError:
                # Pomijamy, jeśli nie da się przekonwertować (np. puste)
                continue
            threshold_data[thresh].append(
                (target_name, typ, source, model_group, score)
            )

# --- ZAPIS WYNIKÓW DO PLIKÓW CSV ---
for threshold, rows in threshold_data.items():
    output_filename = f"gGSP_{threshold}.csv"
    with open(output_filename, 'w', newline='', encoding='utf-8') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Target", "Type", "Source", "ModelGroup", "gGSP"])
        writer.writerows(rows)
    print(f"Zapisano {len(rows)} wierszy dla progu {threshold} do {output_filename}")

print("Przetwarzanie zakończone.")