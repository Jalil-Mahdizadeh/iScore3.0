from iscore3.gate4a.compound_mapping import (
    CompoundIdentity,
    candidate_mapping_state,
    candidate_queries,
    extract_pubchem_properties,
)


def test_alternative_name_precedes_source_codes() -> None:
    identity = CompoundIdentity(1, "AZD-6244/ARRY-886", "Selumetinib")
    assert candidate_queries(identity) == (
        "Selumetinib",
        "AZD-6244/ARRY-886",
        "AZD-6244",
        "ARRY-886",
    )


def test_derivative_identity_is_quarantined_without_base_compound_query() -> None:
    identity = CompoundIdentity(13, "BIBF-1120 (derivative)", "")
    assert candidate_queries(identity) == ()
    assert candidate_mapping_state(identity, ()) == "quarantined_ambiguous_source_identity"


def test_pubchem_response_extraction_requires_structure_and_cid() -> None:
    response = {
        "PropertyTable": {
            "Properties": [
                {"CID": 1, "SMILES": "CCO", "InChIKey": "X"},
                {"CID": 2, "Title": "missing structure"},
            ]
        }
    }
    properties = extract_pubchem_properties(response)
    assert properties == ({"CID": 1, "SMILES": "CCO", "InChIKey": "X"},)
    identity = CompoundIdentity(1, "ethanol", "")
    assert candidate_mapping_state(identity, properties) == "candidate_requires_manual_verification"
