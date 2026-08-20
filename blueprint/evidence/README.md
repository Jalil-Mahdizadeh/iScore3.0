# Evidence Tables

These tab-separated inventories are the audit trail behind the committee prose.

- publications.tsv classifies each reviewed work by task, protein/ligand representation, ligand/pose information, public source, strict relevance, and review depth.
- software.tsv records the effective interface found in source/documentation, ligand-coordinate status, observed licence, reproducibility issues, and reuse decision.
- datasets.tsv records the reviewed release/snapshot, intended role, labels, structure/site coverage, scale, access/licence, caveat, and source.
- leakage_threats.tsv maps each leakage or shortcut mechanism to its control, required audit artifact, and release rule.
- search_log.md documents search sources, example queries, source-code audit procedure, saturation rule, unresolved items, and evidence-quality labels.

Every material claim should trace to an official URL, DOI, repository, or deposit. Counts and licence statements are time-stamped observations and must be reverified when a dependency is pinned.

TSV files use one header row, UTF-8 text, literal tabs, and no embedded newlines. They are validated for a constant field count after edits.
