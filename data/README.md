# Data lifecycle

There is no active iScore3.0 experimental gate. Gate-4A data contracts, ledgers,
manifests, splits, and evidence are frozen at Git tag
`gate4a-terminal-2026-08-22`; the preceding programme is frozen at
`gate3-terminal-2026-08-21`.

Ignored local raw, interim, feature, external, or checkpoint files may remain for
recovery, but they are inactive and must not be interpreted as admitted data for
a future hypothesis. Any successor project requires a new, preregistered data
governance contract, fresh development/test separation, immutable source hashes,
and an explicit rule for whether historical gate data are excluded or used only
as labelled retrospective sensitivity evidence.

No process may overwrite historical source snapshots. New acquisitions require
new versioned paths, manifests, hashes, timestamps, terms, and processing
commands.
