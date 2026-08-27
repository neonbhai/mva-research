"""Scale-measurement tooling.

Nothing in this package is part of the pipeline. It exists to answer one
question before the real data arrives: what does ``mva.ingestion`` actually do
when it is handed a whole-genome callset instead of a twelve-line fixture?

Two modules:

* :mod:`tools.scale.generate_wgs_vcf` fabricates a WGS-shaped VCF. Every base,
  coordinate, depth and genotype in it is invented. It is a **throughput
  phantom**: it is structurally like a real callset so that parser and model
  costs are measured against realistic work, and it is biologically meaningless.
  No conclusion about variants, genes or disease may be drawn from it.
* :mod:`tools.scale.harness` runs the measurements, each in its own subprocess
  so that peak RSS is attributable to one stage rather than to whatever the
  previous stage forgot to free.
"""
