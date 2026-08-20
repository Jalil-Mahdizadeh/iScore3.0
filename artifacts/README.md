# Artifact contract

Each experiment directory is immutable and content-addressed by its input manifests, resolved configuration, code revision, feature schemas, and any upstream checkpoint hashes. Versioned Gate-0/1 artifacts are limited to small metrics, predictions, and provenance summaries. Large checkpoints and logs remain outside Git and are referenced by checksums.
