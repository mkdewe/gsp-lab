"""
fill_inf_from_computed.py

Dwa tryby pracy:
  1. --download  : pobiera z rnapuzzles.org tylko pliki INF_all, INF_wc, INF_nwc, INF_stacking
  2. --fill      : uzupełnia pobrane pliki wartościami z wyników INF_all (naszych obliczeń)
  3. (domyślnie): robi oba kroki po kolei

Logika fill:
  - dopasowanie po Lab+Num (parsowane z nazwy pliku modelu: CR1107_TS029_1.pdb → TS029, 1)
  - zastępuje wartość w pobranym pliku naszą tylko jeśli:
      * pobrany plik ma 0 lub -1
      * AND nasza wartość jest > 0 (dla wc/nwc/stacking) lub wyższa niż pobrana (dla all)
"""

import os
import re
import sys
import argparse
import requests
import pandas as pd
from bs4 import BeautifulSoup

# ── ścieżki ──────────────────────────────────────────────────────────────────
PROJECT_DIR        = r"D:\Studia\Projekty\gsp-lab"
OUTPUT_BASE_DIR    = os.path.join(PROJECT_DIR, "results", "other-metrics", "rnapuzzles.github.io")
COMPUTED_INF_DIR   = os.path.join(PROJECT_DIR, "results", "INF_all")

BASE_URL = "https://www.rnapuzzles.org/table/2000/01/01/"
SUFFIXES = ["", "a", "b", "tRNA", "tBox", "Free", "Bound", "NMR", "noend"]

# Tylko te 4 kolumny nas interesują (dopasowanie case-insensitive)
WANTED_COLS = {"inf all", "inf wc", "inf nwc", "inf stacking"}

# Nazwy plików wyjściowych dla każdej kolumny
COL_FILENAME_MAP = {
    "inf all":      "INF_all.csv",
    "inf wc":       "INF_wc.csv",
    "inf nwc":      "INF_nwc.csv",
    "inf stacking": "INF_stacking.csv",
}

# Kolumna w naszym CSV odpowiadająca każdemu plikowi
COL_COMPUTED_MAP = {
    "inf all":      "best_INF_all",
    "inf wc":       "INF_WC",
    "inf nwc":      "INF_nWC",
    "inf stacking": None,   # nie liczymy stackingu — pomijamy fill
}

# ── helpers ───────────────────────────────────────────────────────────────────

def parse_lab_num_from_filename(model_name):
    """
    CR1107_TS029_1.pdb        → (TS029, 1)
    CR1107_TSR01_3_refined.pdb → (TSR01, 3)
    1_bujnicki_2.pdb           → (bujnicki, 2)  ← PZ style
    """
    stem = os.path.splitext(model_name)[0]
    # usuń suffix _refined / _rpr
    stem = re.sub(r'_(refined|rpr)$', '', stem, flags=re.IGNORECASE)

    # CR style: CR1107_TS029_1  lub  CR1107_TSR01_3
    m = re.search(r'_(TS[A-Z0-9]+)_(\d+)$', stem, re.IGNORECASE)
    if m:
        return m.group(1).upper(), str(int(m.group(2)))

    # PZ style: 1_bujnicki_2
    m = re.match(r'^\d+_([^_]+)_(\d+)$', stem)
    if m:
        return m.group(1).lower(), str(int(m.group(2)))

    return None, None


def is_placeholder(val):
    """Czy wartość jest placeholderem (0 lub -1) wymagającym uzupełnienia."""
    try:
        v = float(val)
        return v == 0.0 or v == -1.0
    except (ValueError, TypeError):
        return False


# ── DOWNLOAD ──────────────────────────────────────────────────────────────────

def fetch_table(puzzle_id):
    url = f"{BASE_URL}{puzzle_id}-3d.html"
    print(f"  GET {url}")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            print(f"    HTTP {r.status_code}")
            return None
        soup = BeautifulSoup(r.content, "html.parser")
        table = soup.find("table", class_="sortable")
        if not table:
            return None
        headers = [th.get_text(strip=True) for th in table.find_all("tr")[0].find_all("th")]
        rows = []
        for row in table.find_all("tr")[1:]:
            cols = [td.get_text(strip=True) for td in row.find_all("td")]
            if cols:
                rows.append(cols)
        if not rows:
            return None
        return pd.DataFrame(rows, columns=headers)
    except Exception as e:
        print(f"    ERROR: {e}")
        return None


def save_inf_columns(df, puzzle_id):
    puzzle_dir = os.path.join(OUTPUT_BASE_DIR, puzzle_id)
    os.makedirs(puzzle_dir, exist_ok=True)

    # upewnij się że Lab i Num są w df
    if "Lab" not in df.columns or "Num" not in df.columns:
        print(f"    Brak kolumn Lab/Num w tabeli dla {puzzle_id}")
        return

    saved = []
    for col in df.columns:
        col_key = col.strip().lower()
        if col_key not in WANTED_COLS:
            continue
        filename = COL_FILENAME_MAP[col_key]
        filepath = os.path.join(puzzle_dir, filename)
        df[["Lab", "Num", col]].to_csv(filepath, index=False)
        saved.append(filename)

    if saved:
        print(f"    Zapisano: {', '.join(saved)}")
    else:
        print(f"    Brak kolumn INF w tabeli dla {puzzle_id}")


def download_all():
    print("\n=== DOWNLOAD ===")

    # PZ series
    i = 1
    while True:
        base_id = f"PZ{i}"
        found = False
        for suffix in SUFFIXES:
            full_id = base_id + suffix
            df = fetch_table(full_id)
            if df is not None and not df.empty:
                print(f"  [{full_id}] znaleziono dane")
                save_inf_columns(df, full_id)
                found = True
        if not found:
            print(f"  Brak danych dla {base_id} — zatrzymuję PZ")
            break
        i += 1

    # CR series
    cr_ids = [
        "CR1107", "CR1108", "CR1116", "CR1117", "CR1126",
        "CR1128", "CR1136", "CR1138", "CR1149", "CR1156",
        "CR1189", "CR1190",
    ]
    for cr_id in cr_ids:
        df = fetch_table(cr_id)
        if df is not None and not df.empty:
            print(f"  [{cr_id}] znaleziono dane")
            save_inf_columns(df, cr_id)
        else:
            print(f"  [{cr_id}] brak danych")


# ── FILL ──────────────────────────────────────────────────────────────────────

def fill_case(case_id):
    computed_path = os.path.join(COMPUTED_INF_DIR, f"{case_id}.csv")
    puzzle_dir    = os.path.join(OUTPUT_BASE_DIR, case_id)

    if not os.path.isfile(computed_path):
        print(f"  [{case_id}] Brak pliku wyników: {computed_path}")
        return
    if not os.path.isdir(puzzle_dir):
        print(f"  [{case_id}] Brak katalogu pobranych danych: {puzzle_dir}")
        return

    try:
        df_comp = pd.read_csv(computed_path)
    except Exception as e:
        print(f"  [{case_id}] Błąd wczytywania {computed_path}: {e}")
        return

    # Zbuduj lookup: (lab_upper, num_str) → row z df_comp
    lookup = {}
    for _, row in df_comp.iterrows():
        lab, num = parse_lab_num_from_filename(str(row.get("model_name", "")))
        if lab and num:
            lookup[(lab.upper(), num)] = row

    if not lookup:
        print(f"  [{case_id}] Brak sparsowanych wpisów w pliku wyników")
        return

    for col_key, filename in COL_FILENAME_MAP.items():
        comp_col = COL_COMPUTED_MAP[col_key]
        if comp_col is None:
            continue  # INF_stacking — nie mamy własnych wyników, pomijamy

        filepath = os.path.join(puzzle_dir, filename)
        if not os.path.isfile(filepath):
            continue

        try:
            df_target = pd.read_csv(filepath)
        except Exception as e:
            print(f"  [{case_id}/{filename}] Błąd wczytywania: {e}")
            continue

        # nazwa kolumny z wartością (trzecia kolumna, po Lab i Num)
        value_cols = [c for c in df_target.columns if c not in ("Lab", "Num")]
        if not value_cols:
            continue
        vcol = value_cols[0]

        changes = 0
        for idx, row in df_target.iterrows():
            lab = str(row["Lab"]).strip().upper()
            num = str(row["Num"]).strip()
            key = (lab, num)

            if key not in lookup:
                continue

            current_val = row[vcol]
            if not is_placeholder(current_val):
                continue  # wartość już jest sensowna — nie nadpisujemy

            new_val_raw = lookup[key].get(comp_col)
            if new_val_raw is None or pd.isna(new_val_raw):
                continue

            try:
                new_val = float(new_val_raw)
            except (ValueError, TypeError):
                continue

            if new_val <= 0:
                continue  # nasza wartość też jest placeholderem — bez sensu nadpisywać

            df_target.at[idx, vcol] = round(new_val, 3)
            changes += 1

        if changes:
            df_target.to_csv(filepath, index=False)
            print(f"  [{case_id}/{filename}] Zaktualizowano {changes} wpisów")
        else:
            print(f"  [{case_id}/{filename}] Brak zmian")


def fill_all():
    print("\n=== FILL ===")
    if not os.path.isdir(COMPUTED_INF_DIR):
        print(f"Brak katalogu z wynikami: {COMPUTED_INF_DIR}")
        return

    cases = [
        os.path.splitext(f)[0]
        for f in os.listdir(COMPUTED_INF_DIR)
        if f.endswith(".csv")
    ]
    print(f"Znaleziono {len(cases)} case(ów) w INF_all")

    for case_id in sorted(cases):
        fill_case(case_id)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="Tylko pobierz dane z rnapuzzles.org")
    parser.add_argument("--fill",     action="store_true", help="Tylko uzupełnij pobrane dane")
    args = parser.parse_args()

    if args.download:
        download_all()
    elif args.fill:
        fill_all()
    else:
        download_all()
        fill_all()

    print("\nGotowe.")

if __name__ == "__main__":
    main()
