"""Recorded Rotterdam endpoint construction and separate-censoring analysis."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from estimation import evaluate

def endpoints(data):
    """Retain recurrence censoring for late deaths and death follow-up separately."""
    n = len(data)
    if n != 2982 or len(np.unique(data["pid"])) != n:
        raise ValueError("Expected 2,982 unique recorded persons")
    rec = data["recur"] == 1
    early = ~rec & (data["death"] == 1) & (data["dtime"] <= data["rtime"])
    late = ~rec & (data["death"] == 1) & (data["dtime"] > data["rtime"])
    pt = np.where(rec, data["rtime"], np.where(early, data["dtime"], data["rtime"]))
    pe = (rec | early).astype(int)
    p = np.c_[pt, pe]
    d = np.c_[data["dtime"], data["death"]]
    if np.sum(late) != 43 or np.any(pt[late] != data["rtime"][late]) or np.any(pe[late] != 0):
        raise ValueError("Recorded endpoint check failed")
    if np.sum(rec) != 1518 or np.sum(early) != 152 or np.sum(pe == 0) != 1312:
        raise ValueError("Endpoint branch totals differ")
    return p, d

def analyze(p, d):
    """Use the declared day horizon and separate censoring; no bootstrap draws."""
    return evaluate(p, d, 1826.25, 1.0, common=False, check_optimality=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_csv", type=Path, help="Locally supplied Rotterdam source fields")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = np.genfromtxt(args.source_csv, names=True, delimiter=",")
    p, d = endpoints(data)
    result = analyze(p, d)
    from encoding import clean
    args.output.write_text(json.dumps(clean(result), indent=2, allow_nan=False) + "\n")

if __name__ == "__main__":
    main()
