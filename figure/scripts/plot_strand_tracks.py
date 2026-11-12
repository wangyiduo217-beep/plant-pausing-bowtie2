#!/usr/bin/env python3
"""Plot observed/predicted plus-above, minus-below signal tracks."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "model/src"))
from plant_pausing_bowtie2.strand_model import main

if __name__ == "__main__":
    raise SystemExit(main(["plot", *sys.argv[1:]]))
