# Third-party source audit

External source trees are cloned into the ignored `source_cache/` directory. The repository versions only audit records containing the canonical URL, exact revision, retrieval date, licence result, files inspected, interfaces relied on, and cryptographic hashes for checkpoints or other executable artifacts.

Absence of a detected licence is treated as **no permission to redistribute or incorporate**, not as an open-source licence.

[`reviewed_sources.tsv`](reviewed_sources.tsv) is the tracked ledger for Gate-0/1 source inspections. It records exact commits, effective interfaces, visible licences, inspected paths, and representative file hashes. Source-cache contents are evidence inputs only and are never imported into iScore3.0 unless a separate reuse decision and licence review permits it.
