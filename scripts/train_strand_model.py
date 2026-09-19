#!/usr/bin/env python3
"""Train one species-specific strand-resolved model."""
from plant_pausing_bowtie2.strand_model import main

if __name__ == "__main__":
    raise SystemExit(main(["train", *__import__("sys").argv[1:]]))
