"""Check stable inputs and small exact examples without simulation draws."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time


def main():
    root = Path(__file__).resolve().parents[1]
    start = time.perf_counter()
    entries = (root / "checksums.sha256").read_text().splitlines()
    for entry in entries:
        digest, name = entry.split("  ", 1)
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"Input checksum does not match: {name}")
    subprocess.run([sys.executable, str(root / "tests/test_method.py")], check=True)
    subprocess.run([sys.executable, str(root / "tests/test_analysis.py")], check=True)
    subprocess.run([sys.executable, str(root / "tests/test_release_reporting.py")], check=True)
    subprocess.run([sys.executable, str(root / "examples/minimal_example.py")], check=True)
    definition = json.loads((root / "data/scenario_definitions.json").read_text())
    if len(definition["cells"]) != 27 or sum(c["replicates"] for c in definition["cells"]) != 39000:
        raise ValueError("The scenario definitions do not match the reported study design.")
    if definition["base_seed"] != 20261004:
        raise ValueError("The scenario seed does not match the reported study design.")
    print(f"Verified {len(entries)} stable scientific checksums.")
    print(f"Deterministic checks completed in {time.perf_counter() - start:.2f} seconds")


if __name__ == "__main__":
    main()
