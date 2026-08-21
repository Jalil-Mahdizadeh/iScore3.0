import csv
from pathlib import Path

from iscore3.gate03.curation import (
    curate_bindingdb_series,
    exact_positive_number,
    normalize_condition,
    publication_id,
)


def test_publication_and_exact_value_contract():
    assert publication_id("https://doi.org/10.1/ABC.", "123") == "doi:10.1/abc"
    assert publication_id("", "123") == "pmid:123"
    assert publication_id("", "") == ""
    assert exact_positive_number(" 25 ") == 25
    assert exact_positive_number(">25") is None
    assert exact_positive_number("0") is None
    assert normalize_condition(" 7.00 ") == "7"
    assert normalize_condition("") == "<missing>"


def _row(index: int, *, kd: str, smiles: str, doi: str = "10.1/example"):
    return {
        "BindingDB Reactant_set_id": str(index),
        "Ligand SMILES": smiles,
        "Ligand InChI Key": f"KEY-{index}",
        "Target Name": "Target",
        "Kd (nM)": kd,
        "pH": "7.0",
        "Temp (C)": "25 C",
        "Curation/DataSource": "BindingDB",
        "Article DOI": doi,
        "PMID": "123",
        "Authors": "A",
        "Date of publication": "2020",
        "PDB ID(s) for Ligand-Target Complex": "1ABC",
        "Number of Protein Chains in Target (>1 implies a multichain complex)": "1",
        "BindingDB Target Chain Sequence 1": "ACDEFGHIK",
        "PDB ID(s) of Target Chain 1": "1ABC",
        "UniProt (SwissProt) Primary ID of Target Chain 1": "P12345",
        "UniProt (TrEMBL) Primary ID of Target Chain 1": "",
    }


def test_series_curation_requires_publication_and_preserves_rows(tmp_path: Path):
    rows = [
        _row(1, kd="10", smiles="CC"),
        _row(2, kd="100", smiles="CCC"),
        _row(3, kd="20", smiles="CCCC", doi=""),
    ]
    rows[-1]["PMID"] = ""
    path = tmp_path / "bindingdb.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summaries, observations, audit = curate_bindingdb_series(
        bindingdb_tsv=path, minimum_ligands=2
    )
    assert len(summaries) == 1
    assert len(observations) == 2
    assert summaries[0]["ligand_count"] == "2"
    assert audit["row_census"]["excluded_missing_stable_publication"] == 1
    assert {row["source_reactant_set_ids"] for row in observations} == {"1", "2"}
