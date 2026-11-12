#!/usr/bin/env python3
"""Generate raw plus/minus genome tracks from a trained model."""
from plant_pausing_bowtie2.strand_model import main

if __name__ == "__main__":
    raise SystemExit(main(["predict", *__import__("sys").argv[1:]]))
