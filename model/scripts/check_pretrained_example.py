#!/usr/bin/env python3
"""Run actual model inference and compare every example score to its reference."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.pretrained import main

if __name__ == "__main__":
    raise SystemExit(main(["check-example", *sys.argv[1:]]))
