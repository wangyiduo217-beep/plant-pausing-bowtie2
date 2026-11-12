#!/usr/bin/env python3
"""Predict plus/minus window scores for 1,024-bp forward-reference DNA."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.pretrained import main

if __name__ == "__main__":
    raise SystemExit(main(["predict", *sys.argv[1:]]))
