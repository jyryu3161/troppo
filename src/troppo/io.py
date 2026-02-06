"""
Data I/O utilities for the troppo package.

Provides functions for loading expression data, metabolic models,
detecting gene ID nomenclatures, and converting between gene ID types.
"""

import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import pandas as pd

# Canonical gene ID type constants
ENTREZ_ID = "entrez_id"
ENSEMBL_GENE_ID = "ensembl_gene_id"
HGNC_ID = "hgnc_id"
SYMBOL = "symbol"

_ID_PATTERNS = [
    (ENTREZ_ID, re.compile(r"^\d+$")),
    (ENSEMBL_GENE_ID, re.compile(r"^ENSG\d+", re.IGNORECASE)),
    (HGNC_ID, re.compile(r"^HGNC:", re.IGNORECASE)),
]


# ---------------------------------------------------------------------------
# Gene ID detection
# ---------------------------------------------------------------------------

def detect_gene_id_type(ids: Sequence[str]) -> str:
    """Detect gene identifier nomenclature via pattern-based majority vote.

    Samples up to the first 20 IDs and classifies each as one of:
    ``entrez_id``, ``ensembl_gene_id``, ``hgnc_id``, or ``symbol``.

    Parameters
    ----------
    ids : Sequence[str]
        Gene identifiers to classify.

    Returns
    -------
    str
        Detected ID type string.
    """
    sample = [str(i).strip() for i in list(ids)[:20] if pd.notna(i)]
    votes: List[str] = []
    for gid in sample:
        matched = False
        for id_type, pattern in _ID_PATTERNS:
            if pattern.match(gid):
                votes.append(id_type)
                matched = True
                break
        if not matched:
            votes.append(SYMBOL)
    if not votes:
        return SYMBOL
    return Counter(votes).most_common(1)[0][0]


# ---------------------------------------------------------------------------
# Expression data loading
# ---------------------------------------------------------------------------

def load_expression_csv(
    path: Union[str, Path],
    gene_col: Optional[str] = None,
    sample_cols: Optional[List[str]] = None,
    sep: str = "auto",
    gene_id_type: str = "auto",
) -> pd.DataFrame:
    """Load a CSV/TSV expression matrix and return a tidy DataFrame.

    Parameters
    ----------
    path : str or Path
        Path to the expression file.
    gene_col : str, optional
        Column name containing gene identifiers.  If *None*, the first
        column is used (or the existing index if it is non-numeric).
    sample_cols : list of str, optional
        Columns to keep as sample data.  If *None*, all columns except
        *gene_col* are used.
    sep : str
        Column separator.  ``"auto"`` sniffs for tab or comma.
    gene_id_type : str
        Gene ID nomenclature (``"auto"`` to detect automatically).

    Returns
    -------
    pd.DataFrame
        Genes as index (``index.name`` set to detected/specified ID type),
        samples as columns.
    """
    path = Path(path)

    # Auto-detect separator
    if sep == "auto":
        with open(path, "r") as fh:
            first_line = fh.readline()
        sep = "\t" if "\t" in first_line else ","

    df = pd.read_csv(path, sep=sep)

    # Resolve gene column
    if gene_col is None:
        gene_col = df.columns[0]
    df = df.set_index(gene_col)

    # Subset to requested samples
    if sample_cols is not None:
        df = df[sample_cols]

    # Coerce expression values to numeric where possible
    df = df.apply(pd.to_numeric, errors="coerce")

    # Detect / assign gene ID type
    if gene_id_type == "auto":
        gene_id_type = detect_gene_id_type(df.index)
    df.index.name = gene_id_type

    return df


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

_BIOMASS_PATTERNS = re.compile(r"biomass|growth", re.IGNORECASE)


def load_model(path: Union[str, Path], biomass_reaction: str = "auto"):
    """Load a metabolic model from SBML or JSON.

    Parameters
    ----------
    path : str or Path
        Path to an SBML (``.xml``) or JSON (``.json``) model file.
    biomass_reaction : str
        Reaction ID to set as the objective.  ``"auto"`` searches
        reaction IDs and names for *biomass* or *growth*.

    Returns
    -------
    cobra.Model
        The loaded COBRApy model with objective set (if found).
    """
    import cobra  # heavyweight optional dependency

    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".json":
        model = cobra.io.load_json_model(str(path))
    elif suffix in (".xml", ".sbml"):
        model = cobra.io.read_sbml_model(str(path))
    else:
        raise ValueError(f"Unsupported model format: {suffix}")

    if biomass_reaction == "auto":
        for rxn in model.reactions:
            if _BIOMASS_PATTERNS.search(rxn.id) or _BIOMASS_PATTERNS.search(rxn.name):
                model.objective = rxn.id
                break
    elif biomass_reaction is not None:
        model.objective = biomass_reaction

    return model


# ---------------------------------------------------------------------------
# Gene ID conversion
# ---------------------------------------------------------------------------

def convert_gene_ids(
    ids: Sequence[str],
    from_type: str,
    to_type: str,
    report: bool = True,
) -> Dict[str, str]:
    """Convert gene identifiers between nomenclatures via HGNC.

    Wraps :func:`troppo.omics.id_converter.idConverter`.

    Parameters
    ----------
    ids : Sequence[str]
        Gene identifiers to convert.
    from_type : str
        Source nomenclature (e.g. ``"symbol"``, ``"entrez_id"``).
    to_type : str
        Target nomenclature.
    report : bool
        If *True*, print conversion success/failure statistics.

    Returns
    -------
    dict
        Mapping ``{original_id: converted_id}``.
    """
    from troppo.omics.id_converter import idConverter

    id_list = list(ids)
    mapping = idConverter(id_list, from_type, to_type) or {}

    if report:
        total = len(id_list)
        converted = len(mapping)
        failed = total - converted
        print(
            f"Gene ID conversion ({from_type} -> {to_type}): "
            f"{converted}/{total} succeeded, {failed} failed"
        )

    return mapping
