# Gate-2 frozen configuration

`feasibility.yaml` is the original pre-fit protocol. Amendments 01 and 02 record
only dataset mapping/minimum-count decisions made before any Gate-2 outcome fit.
`feasibility-effective-v3.yaml` is the canonical frozen configuration used for the
completed run; it records the v3 strict data, structural leakage graph, fixed
encoders, nested component/OOD evaluation, low-capacity tensor interaction,
negative controls, uncertainty rule, and progression criterion.

The Gate-2 progression check is `FAIL` / no-go. These files must not be altered
to rationalize the outcome; any future study requires a new phase identifier and
a new pre-fit freeze.
