#!/usr/bin/env python3
"""Script entry point for the configuration-driven SRA/FASTQ-to-BAM workflow."""

from plant_pausing_bowtie2.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
