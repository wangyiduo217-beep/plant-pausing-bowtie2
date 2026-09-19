#!/usr/bin/env python3
"""Build strand-resolved consensus intervals and sliding-window y tables."""

import sys

from plant_pausing_bowtie2.y_labels import main


if __name__ == "__main__":
    raise SystemExit(main(["build", *sys.argv[1:]]))
