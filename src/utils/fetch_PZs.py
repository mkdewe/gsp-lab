import os
import requests
import pandas as pd
from bs4 import BeautifulSoup

BASE_URL = "https://www.rnapuzzles.org/table/2000/01/01/"
SUFFIXES = ["", "a", "b", "tRNA", "tBox", "Free", "Bound"]

# Output directory: results/other_metrics/rnapuzzles.github.io
PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
OUTPUT_BASE_DIR = os.path.join(PROJECT_DIR, "results", "other-metrics", "rnapuzzles.github.io")

def clean_column_name(name):
    return name.replace(" ", "_").replace("/", "_")

def extract_table_data(puzzle_id):
    url = f"{BASE_URL}{puzzle_id}-3d.html"
    print(f"Trying to fetch: {url}")
    response = requests.get(url)
    if response.status_code != 200:
        return None
    soup = BeautifulSoup(response.content, "html.parser")
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
        print(f"Saved: {filepath}")

    index_file = os.path.join(puzzle_dir, "index.csv")
    pd.DataFrame({"Column": [c for c in df.columns if c.lower() not in ["lab", "num", "detail", "best sol."]]}).to_csv(index_file, index=False)
    print(f"Index saved: {index_file}")


def main():
    i = 1
    while True:
        base_id = f"PZ{i}"
        print(f"\nChecking group: {base_id}")
        found_data_in_group = False

        for suffix in SUFFIXES:
            full_id = base_id + suffix
            try:
                df = extract_table_data(full_id)
                if df is not None and not df.empty:
                    print(f"Found data for {full_id}, saving...")
                    save_columns_to_csv(df, full_id)
                    found_data_in_group = True
                else:
                    print(f"No data for {full_id}")
            except Exception as e:
                print(f"Error with {full_id}: {e}")

        if not found_data_in_group:
            print(f"No data found for group {base_id}. Stopping further search.")
            break

        i += 1

if __name__ == "__main__":
    main()
