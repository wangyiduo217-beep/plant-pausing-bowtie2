#!/usr/bin/env python3
"""Download and verify publication weights (no model dependencies required)."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.pretrained import main

if __name__ == "__main__":
    raise SystemExit(main(["download", *sys.argv[1:]]))
