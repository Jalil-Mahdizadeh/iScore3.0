# Gate-4A split artifacts

This directory contains immutable, outcome-blind row-to-fold assignments and
leakage diagnostics.

The primary evaluation is double-cold: every test pair contains both an unseen
ligand scaffold component and an unseen receptor binding-site structural/family
component. All stereochemical/state variants of a ligand remain together. All
mutants, constructs, close homologues, and structurally similar pockets remain
together.

Random-pair splits are forbidden. Fold construction, exclusion reasons, source
publication edges, ligand similarity edges, sequence edges, and pocket-structure
edges must be frozen before outcome-model fitting.
