#!/usr/bin/env python3
"""Call per-SRX stranded GRO-seq transcript intervals from accepted BAM files."""

import sys

from plant_pausing_bowtie2.y_labels import main


if __name__ == "__main__":
    raise SystemExit(main(["call", *sys.argv[1:]]))
