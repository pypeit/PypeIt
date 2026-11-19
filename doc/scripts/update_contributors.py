
"""
Update the contributors list.  This should *not* be executed every time the
documentation is updated.  Instead, it should be run periodically with intention
so that the contributors.csv file is updated correctly.

This script was written almost entirely by Claude.  K. Westfall reviewed it,
applied some modest edits, and tested it.
"""
import requests
import time
import csv
from pathlib import Path

from IPython import embed

def load_known_contributors(csv_path):
    """Load the existing contributors CSV into a dict keyed by GitHub handle."""
    known = {}
    _csv_path = Path(csv_path).absolute()
    if not _csv_path.is_file():
        raise FileNotFoundError(f'{_csv_path} does not exist!')
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or len(row) < 3:
                continue
            handle, name, valid = row[0], row[1], row[2]
            known[handle] = {"name": name, "valid": bool(valid)}
    return known

def get_contributors(owner, repo, csv_path="contributors.csv",
                      outfile="contributors.txt", token=None):
    headers = {"Authorization": f"token {token}"} if token else {}

    known = load_known_contributors(csv_path)

    new_entries = []  # rows to append to the CSV

    # Step 1: get all contributors (login + contribution count) from the repo
    contributors = []
    page = 1
    while True:
        url = f"https://api.github.com/repos/{owner}/{repo}/contributors"
        resp = requests.get(url, params={"per_page": 100, "page": page}, headers=headers)
        resp.raise_for_status()
        batch = resp.json()
        if not batch:
            break
        contributors.extend(batch)
        page += 1

    # Step 2: resolve names, using the CSV cache wherever possible
    resolved = []
    for c in contributors:
        login = c["login"]

        if login in known:
            name = known[login]["name"]
            valid = known[login]["valid"]
        else:
            # Not in our records yet -- look it up via the API
            user_resp = requests.get(f"https://api.github.com/users/{login}", headers=headers)
            user_resp.raise_for_status()
            name = user_resp.json().get("name") or ""
            valid = True
            new_entries.append([login, name, valid])
            known[login] = {"name": name, "valid": valid}  # avoid dup lookups this run
            time.sleep(0.1)  # be polite to the API

        resolved.append({
            "login": login,
            "name": name,
            "valid": valid,
            "contributions": c["contributions"],
        })

    # Step 3: append any newly discovered contributors to the CSV
    if new_entries:
        with open(csv_path, "a", newline="") as f:
            writer = csv.writer(f, lineterminator="\n")
            writer.writerows(new_entries)
        print(f"Added {len(new_entries)} new contributor(s) to {csv_path}")

    # Step 4: write the full contributor report (skip anyone flagged invalid)
    with open(outfile, "w") as f:
        f.write("login\tname\tcontributions\n")
        for r in resolved:
            if r["valid"] == "invalid":
                continue
            f.write(f"{r['login']}    {r['name']}    {r['contributions']}\n")

    print(f"Wrote report to {outfile} ({len(new_entries)} new lookups, "
          f"{len(resolved) - len(new_entries)} from cache)")


if __name__ == '__main__':
    get_contributors("pypeit", "PypeIt")
