from iscore3.gate4a.source_audit import audit_matrix_values


def test_matrix_audit_preserves_blank_censoring() -> None:
    header = ["Accession Number", "Entrez Gene Symbol", "Kinase", "Drug A", "Drug B"]
    rows = [
        ["NP_1", "GENE1", "KIN1", "10", ""],
        ["NP_2", "GENE2", "KIN2", "10000", "500"],
    ]
    audit = audit_matrix_values(
        header,
        rows,
        metadata_columns=3,
        source_sha256="fixture",
        sheet_name="fixture",
        blank_is_censored=True,
        censor_limit_nm=10_000.0,
    )
    assert audit.pair_count == 4
    assert audit.exact_count == 3
    assert audit.right_censored_count == 1
    assert audit.exact_numeric_10000_count == 1
    assert audit.invalid_count == 0
