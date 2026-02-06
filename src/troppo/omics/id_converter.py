"""
Gene ID converter using the HGNC complete set database.

Provides functions for:
- idConverter: Convert gene IDs between nomenclatures
- searchNomenclature: Detect which nomenclature a list of IDs belongs to

Uses file-based caching to avoid redundant downloads.
"""

import os
import pandas as pd
import urllib.request
import urllib.error
import shutil
from pathlib import Path
from datetime import date
from functools import lru_cache

# Module-level cache for the HGNC DataFrame
_hgnc_cache = None


def _get_cache_dir() -> Path:
    """Get cache directory for HGNC data."""
    cache_dir = Path.home() / '.troppo' / 'cache'
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _get_HGNC() -> str:
    """
    Download the HGNC complete set file (cached daily).

    Returns
    -------
    str
        Path to the downloaded HGNC TSV file.
    """
    now = date.today()
    cache_dir = _get_cache_dir()
    path = cache_dir / f"hgnc_complete_set_{now}.tsv"

    if path.is_file():
        return str(path)

    url = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"

    try:
        with urllib.request.urlopen(url) as response, open(path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
    except (urllib.error.URLError, OSError) as e:
        # Try to find any existing cached file
        existing = sorted(cache_dir.glob("hgnc_complete_set_*.tsv"))
        if existing:
            print(f'Network error, using cached HGNC file: {existing[-1].name}')
            return str(existing[-1])
        # Fallback to local file in working directory
        local = Path("hgnc_complete_set_2018-09-13.tsv")
        if local.is_file():
            print('Please check internet connection, using locally available HGNC file')
            return str(local)
        raise RuntimeError(f"Cannot download HGNC data and no cache available: {e}")

    return str(path)


def _load_hgnc_df() -> pd.DataFrame:
    """Load HGNC data as a cached DataFrame."""
    global _hgnc_cache
    if _hgnc_cache is not None:
        return _hgnc_cache
    file = _get_HGNC()
    _hgnc_cache = pd.read_csv(file, sep='\t', low_memory=False)
    return _hgnc_cache


def idConverter(ids: list, old: str, new: str) -> dict:
    """
    Convert gene IDs from one nomenclature to another using HGNC.

    Parameters
    ----------
    ids : list or set
        IDs to convert
    old : str
        Source nomenclature column name in HGNC
    new : str
        Target nomenclature column name in HGNC

    Returns
    -------
    dict or None
        {old_id: new_id} mapping, or None on error
    """
    ds = _load_hgnc_df()
    try:
        d = dict(zip(ds[old.lower()], ds[new.lower()]))
    except KeyError:
        print('The ID designation is incorrect, must be one of:\n',
              ds.columns.values.tolist())
        return None

    res = {x: str(d[x]).split('.')[0] for x in ids if x in d}
    return res


def searchNomenclature(ids: list) -> str:
    """
    Detect gene ID nomenclature using pandas-based matching (efficient).

    Parameters
    ----------
    ids : list
        List of gene IDs (all using the same nomenclature)

    Returns
    -------
    str or None
        The nomenclature designation according to HGNC.
    """
    ds = _load_hgnc_df()
    ids_to_check = list(ids)  # copy to avoid mutation

    # Sample up to 20 IDs for efficiency
    sample_ids = ids_to_check[:min(20, len(ids_to_check))]

    matches = {}

    for test_id in sample_ids:
        test_str = str(test_id)
        for col_idx, col in enumerate(ds.columns):
            col_values = ds[col].astype(str).values
            if test_str in col_values:
                matches[col] = matches.get(col, 0) + 1
                break  # found in first matching column

    if not matches:
        print('No match was found for the provided ids')
        return None

    # Return the column with the most matches
    best_col = max(matches, key=matches.get)
    return best_col


if __name__ == "__main__":
    a = ['A1BG', 'HGNC:1']
    print(idConverter(a, 'symbol', 'entrez_id'))
