#!/usr/bin/env python3
"""Plot observed/predicted plus-above, minus-below signal tracks."""
from plant_pausing_bowtie2.strand_model import main

if __name__ == "__main__":
    raise SystemExit(main(["plot", *__import__("sys").argv[1:]]))
