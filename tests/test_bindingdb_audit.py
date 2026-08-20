from iscore3.data.bindingdb_audit import _candidate_assessment, _pdb_ids, reconcile


def _candidate(**updates):
    row = {
        "BindingDB Reactant_set_id": "123",
        "Ligand InChI Key": "AAAA-BBBB-C",
        "Kd (nM)": "25",
        "Curation/DataSource": "BindingDB",
        "Article DOI": "10.1000/example",
        "PMID": "123456",
        "Date of publication": "2020",
        "Date in BindingDB": "2021-01-01",
        "Ligand HET ID in PDB": "LIG",
        "PDB ID(s) for Ligand-Target Complex": "1ABC, 2DEF",
        "UniProt (SwissProt) Primary ID of Target Chain 1": "P12345",
    }
    row.update(updates)
    return row


def _pilot(**updates):
    row = {
        "pdb_id": "1ABC",
        "ligand_comp_id": "LIG",
        "value_nm": "25.0",
        "inchikey": "AAAA-BBBB-C",
        "uniprot_accession": "P12345",
        "citation_doi": "10.1000/EXAMPLE",
        "citation_pubmed": "123456",
        "observation_id": "obs-1",
        "role": "supervised_s0",
        "label_quarantined": "False",
    }
    row.update(updates)
    return row


def test_pdb_id_parser_is_token_aware():
    assert _pdb_ids("1abc, 2DEF") == {"1ABC", "2DEF"}
    assert "1ABC" not in _pdb_ids("X1ABCY")


def test_candidate_assessment_requires_all_provenance_axes():
    result = _candidate_assessment(_pilot(), _candidate())
    assert result["high_confidence"] is True
    assert result["kd_match"] is True

    result = _candidate_assessment(_pilot(), _candidate(**{"Article DOI": "10.1000/wrong", "PMID": "999"}))
    assert result["same_structure_measurement_publication"] is False
    assert result["measurement_publication_present"] is True
    assert result["high_confidence"] is True


def test_censored_kd_fails_high_confidence():
    result = _candidate_assessment(_pilot(), _candidate(**{"Kd (nM)": ">25"}))
    assert result["kd_match"] is False
    assert result["high_confidence"] is False


def test_reference_label_is_not_exported_or_required():
    pilot = _pilot(
        role="site_reference_only",
        value_nm="",
        label_quarantined="True",
    )
    rows, summary = reconcile([pilot], [_candidate(**{"Kd (nM)": "999"})])
    assert rows[0]["kd_exact_match"] == ""
    assert rows[0]["label_quarantine_preserved"] is True
    assert summary["reference_label_quarantine_preserved"] == 1
