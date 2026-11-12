#!/usr/bin/env python3
"""Prepare chromosome-separated strand-model manifests."""
from plant_pausing_bowtie2.strand_model import main

if __name__ == "__main__":
    raise SystemExit(main(["prepare", *__import__("sys").argv[1:]]))
