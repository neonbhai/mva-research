"""Generators for the Track 2 documents that must not drift from the pipeline.

A rejection list maintained by hand agrees with the pipeline on the day it is
written and never again. Everything in here renders from a live triage run, so the
document cannot say something the code does not.
"""

from __future__ import annotations
