#!/usr/bin/env python3
"""Run interval calling and y-table generation in one command."""

import sys

from plant_pausing_bowtie2.y_labels import main


if __name__ == "__main__":
    raise SystemExit(main(["run", *sys.argv[1:]]))
