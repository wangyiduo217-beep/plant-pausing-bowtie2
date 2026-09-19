"""Reproducible Bowtie2 processing for plant nascent-transcription reads."""

__version__ = "0.3.0"

from .workflow import run_sample

__all__ = ["run_sample", "__version__"]
