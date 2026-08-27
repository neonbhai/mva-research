# Tech debt

Deferred work with the cost of deferring it and the condition that should trigger
paying it down. Ordered by that cost, not by effort.

| ID | Item | Cost of deferral | Trigger |
|---|---|---|---|
| TD-01 | Annotation is a local TSV, not VEP/gnomAD/ClinVar | Demo consequence and frequency values are fabricated. No real-data claim is possible until this is replaced. | Before the first real-data run. |
| TD-02 | No known-answer validation set | Sensitivity on real cases is unmeasured; weights are uncalibrated. | Before quoting any performance number. |
| TD-03 | No structural/copy-number variant calling | In a chromosomal-instability disorder, an SV or CNV is a plausible causal class we cannot see at all. | Before treating a negative result as informative. |
| TD-04 | Phenotype scoring has no ontology-graph propagation | An observed child term does not credit its associated parent term, so true matches are under-scored. | When a real HPO release is wired in. |
| TD-05 | Cross-platform determinism untested | Byte-identity is verified on macOS/arm64 only; the container exists but is unexercised. | Before a third party reruns the submission. |
| TD-06 | Network denial is Python-level only | C extensions and subprocesses can still reach the network; the guard is a tripwire, not a boundary. | Before the first real-data run — pair with an OS-level control. |
| TD-07 | No trio/segregation support | Phase stays UNKNOWN for every pair, capping the inheritance component. This is the single largest score-limiting gap. | If parental data becomes available. |
| TD-08 | Drug catalogue is synthetic | Track 2 output demonstrates the reasoning machinery, not a real repurposing candidate. | Before any external presentation of a specific compound. |
| TD-09 | No human domain-expert review | Model-based adversarial review is a filter, not a substitute. | Before submission. |
| TD-11 | Snakemake rules re-run the stage prefix | `mva run <stage>` executes everything up to that stage, so a full DAG run recomputes early stages once per rule (~20s vs ~3s for `mva run all`). Correctness is unaffected. | If DAG runtime becomes a constraint; fix is a resume-from-artifacts CLI mode. |
| TD-12 | Container base image is not digest-pinned | `workflow/containers/Dockerfile` carries a `<PASTE HERE>` placeholder rather than a fabricated digest. The image is reproducible by tag only. | Before using the container for a reproducibility claim. |
| TD-10 | `epcr` is a rank-ordering confidence, not a calibrated probability | The challenge's F-max metric sweeps EPCR thresholds, so miscalibration costs points even when the ranking is right. | If F-max matters more than rank points. |
