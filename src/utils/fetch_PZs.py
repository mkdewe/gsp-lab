import os
import requests
import pandas as pd
from bs4 import BeautifulSoup

BASE_URL = "https://www.rnapuzzles.org/table/2000/01/01/"
SUFFIXES = ["", "a", "b", "tRNA", "tBox", "Free", "Bound", "NMR", "noend"]


# Output directory: results/other_metrics/rnapuzzles.github.io
PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
OUTPUT_BASE_DIR = os.path.join(PROJECT_DIR, "results", "other-metrics", "rnapuzzles.github.io")

def clean_column_name(name):
    return name.replace(" ", "_").replace("/", "_")

def extract_table_data(puzzle_id):
    url = f"{BASE_URL}{puzzle_id}-3d.html"
    print(f"Trying to fetch: {url}")
    try:
        response = requests.get(url)
        if response.status_code != 200:
            print(f"  HTTP {response.status_code} for {url}")
            return None
        soup = BeautifulSoup(response.content, "html.parser")
        table = soup.find("table", class_="sortable")
        if not table:
            print(f"  No table found for {puzzle_id}")
            return None

        headers = [th.get_text(strip=True) for th in table.find_all("tr")[0].find_all("th")]
        rows = []
        for row in table.find_all("tr")[1:]:
            cols = [td.get_text(strip=True) for td in row.find_all("td")]
            if cols:
                rows.append(cols)

        if not rows:
            print(f"  No data rows for {puzzle_id}")
            return None

        return pd.DataFrame(rows, columns=headers)
    except Exception as e:
        print(f"  Error fetching {puzzle_id}: {e}")
        return None

def save_columns_to_csv(df, puzzle_id):
    puzzle_dir = os.path.join(OUTPUT_BASE_DIR, puzzle_id)
    os.makedirs(puzzle_dir, exist_ok=True)

    for col in df.columns:
        if col.lower() in ["lab", "num", "detail", "best sol."]:
            continue

        cleaned_col = clean_column_name(col)
        filename = f"{cleaned_col}.csv"  
        filepath = os.path.join(puzzle_dir, filename)

        df_subset = df[["Lab", "Num", col]]
        df_subset.to_csv(filepath, index=False)
        print(f"  Saved: {filepath}")

    index_file = os.path.join(puzzle_dir, "index.csv")
    pd.DataFrame({"Column": [c for c in df.columns if c.lower() not in ["lab", "num", "detail", "best sol."]]}).to_csv(index_file, index=False)
    print(f"  Index saved: {index_file}")

def process_pz_series():
    """Process PZ series with suffixes"""
    i = 1
    while True:
        base_id = f"PZ{i}"
        print(f"\nChecking PZ group: {base_id}")
        found_data_in_group = False

        for suffix in SUFFIXES:
            full_id = base_id + suffix
            df = extract_table_data(full_id)
            if df is not None and not df.empty:
                print(f"Found data for {full_id}, saving...")
                save_columns_to_csv(df, full_id)
                found_data_in_group = True
            else:
                print(f"No data for {full_id}")

        if not found_data_in_group:
            print(f"No data found for group {base_id}. Stopping PZ series search.")
            break

        i += 1

def process_cr_list():
    """Process specific CR IDs without suffixes"""
    cr_ids = [
        "CR1107", "CR1108", "CR1116", "CR1117", "CR1126", 
        "CR1128", "CR1136", "CR1138", "CR1149", "CR1156", 
        "CR1189", "CR1190"
    ]
    
    print(f"\nProcessing CR series...")
    for cr_id in cr_ids:
        print(f"\nChecking CR ID: {cr_id}")
        df = extract_table_data(cr_id)
        if df is not None and not df.empty:
            print(f"Found data for {cr_id}, saving...")
            save_columns_to_csv(df, cr_id)
        else:
            print(f"No data for {cr_id}")

def process_all():
    """Process both PZ series and specific CR IDs"""
    print("Starting data download from RNA Puzzles...")
    
    # Process PZ series
    process_pz_series()
    
    # Process CR IDs
    process_cr_list()
    
    print("\nAll data processing completed!")

def main():
    # You can choose to run either specific function or all
    process_all()
    
    # Or run them separately:
    # process_pz_series()  # For PZ series only
    # process_cr_list()    # For CR IDs only

if __name__ == "__main__":
    main()