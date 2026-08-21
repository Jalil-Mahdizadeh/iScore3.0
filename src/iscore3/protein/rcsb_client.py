"""Small gate-neutral RCSB API client used by receptor-view acquisition."""

from __future__ import annotations

from http.client import RemoteDisconnected
import json
import time
from typing import Any, Iterator, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/graphql"
COORDINATE_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
USER_AGENT = "iScore3.0/0.1 (scientific-reproducibility; RCSB-cached-requests)"

ENTRY_QUERY = r"""
query IScore3Entries($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_accession_info { deposit_date initial_release_date revision_date }
    rcsb_binding_affinity { comp_id type value unit }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    rcsb_primary_citation { title pdbx_database_id_DOI pdbx_database_id_PubMed }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity { pdbx_description pdbx_mutation pdbx_fragment }
      rcsb_polymer_entity_container_identifiers {
        entity_id asym_ids auth_asym_ids
        reference_sequence_identifiers { database_accession database_name }
      }
      rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    }
    nonpolymer_entities {
      rcsb_id
      rcsb_nonpolymer_entity_container_identifiers {
        entity_id asym_ids auth_asym_ids nonpolymer_comp_id
      }
      nonpolymer_comp {
        chem_comp { id type formula_weight formula }
        rcsb_chem_comp_descriptor { SMILES SMILES_stereo InChI InChIKey }
      }
    }
  }
}
""".strip()


class RcsbClientError(RuntimeError):
    """Raised when the RCSB service response violates the client contract."""


def chunks(values: Sequence[str], size: int) -> Iterator[list[str]]:
    if size <= 0:
        raise ValueError("chunk size must be positive")
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def post_json(url: str, payload: Mapping[str, Any], *, attempts: int = 4) -> dict[str, Any]:
    body = json.dumps(payload, separators=(",", ":"), allow_nan=False).encode("utf-8")
    request = Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            with urlopen(request, timeout=120) as response:
                result = json.loads(response.read().decode("utf-8"))
            if not isinstance(result, dict):
                raise RcsbClientError(f"non-object JSON response from {url}")
            if result.get("errors"):
                raise RcsbClientError(f"Graph/API errors from {url}: {result['errors']}")
            return result
        except (
            HTTPError,
            URLError,
            TimeoutError,
            RemoteDisconnected,
            ConnectionResetError,
            json.JSONDecodeError,
        ) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RcsbClientError(f"request failed after {attempts} attempts: {last_error}")
