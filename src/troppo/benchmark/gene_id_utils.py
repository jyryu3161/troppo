"""
Gene ID Utilities for Benchmark Framework

This module provides utilities for handling different gene ID nomenclatures
in benchmark validations.

다양한 gene ID nomenclature를 처리하는 유틸리티
"""

from typing import List, Set, Dict, Optional
import warnings


def detect_id_type(gene_ids: List[str]) -> str:
    """
    Detect the type of gene IDs based on patterns

    Parameters
    ----------
    gene_ids : List[str]
        List of gene IDs

    Returns
    -------
    str
        Detected ID type ('entrez', 'ensembl', 'symbol', 'hgnc', 'unknown')
    """
    if not gene_ids:
        return 'unknown'

    # Sample a few IDs
    sample_size = min(10, len(gene_ids))
    sample = gene_ids[:sample_size]

    # Count patterns
    entrez_count = sum(1 for g in sample if isinstance(g, (int, float)) or (isinstance(g, str) and g.isdigit()))
    ensembl_count = sum(1 for g in sample if isinstance(g, str) and g.startswith('ENSG'))
    hgnc_count = sum(1 for g in sample if isinstance(g, str) and g.startswith('HGNC:'))

    # Determine type
    if entrez_count >= sample_size * 0.7:
        return 'entrez_id'
    elif ensembl_count >= sample_size * 0.7:
        return 'ensembl_gene_id'
    elif hgnc_count >= sample_size * 0.7:
        return 'hgnc_id'
    else:
        return 'symbol'  # Default to gene symbol


def convert_gene_ids(
    gene_ids: List[str],
    from_type: Optional[str] = None,
    to_type: str = 'symbol',
    use_hgnc: bool = True
) -> Dict[str, str]:
    """
    Convert gene IDs from one nomenclature to another

    Parameters
    ----------
    gene_ids : List[str]
        List of gene IDs to convert
    from_type : str, optional
        Source ID type (auto-detect if None)
    to_type : str
        Target ID type
    use_hgnc : bool
        Whether to use HGNC database for conversion

    Returns
    -------
    Dict[str, str]
        Mapping from original IDs to converted IDs
    """
    if from_type is None:
        from_type = detect_id_type(gene_ids)
        print(f"Auto-detected gene ID type: {from_type}")

    if from_type == to_type:
        return {g: g for g in gene_ids}

    if use_hgnc:
        try:
            from troppo.omics.id_converter import idConverter

            # Convert IDs
            converted = idConverter(gene_ids, from_type, to_type)

            if converted:
                print(f"Converted {len(converted)}/{len(gene_ids)} genes from {from_type} to {to_type}")
                return converted
            else:
                warnings.warn(f"Failed to convert IDs from {from_type} to {to_type}")
                return {}
        except Exception as e:
            warnings.warn(f"ID conversion failed: {str(e)}")
            return {}
    else:
        warnings.warn("HGNC-based conversion disabled, returning empty mapping")
        return {}


def standardize_gene_ids(
    gene_ids: List[str],
    target_type: str = 'symbol'
) -> Dict[str, str]:
    """
    Standardize a list of gene IDs to a common nomenclature

    Parameters
    ----------
    gene_ids : List[str]
        List of gene IDs (possibly mixed types)
    target_type : str
        Target nomenclature

    Returns
    -------
    Dict[str, str]
        Mapping from original to standardized IDs
    """
    # Clean up IDs
    cleaned_ids = [str(g).strip() for g in gene_ids if g]

    # Detect and convert
    return convert_gene_ids(cleaned_ids, to_type=target_type)


def match_gene_lists(
    list1: List[str],
    list2: List[str],
    auto_convert: bool = True,
    target_type: str = 'symbol'
) -> tuple:
    """
    Match two gene lists, handling different ID nomenclatures

    Parameters
    ----------
    list1 : List[str]
        First gene list
    list2 : List[str]
        Second gene list
    auto_convert : bool
        Whether to auto-convert IDs
    target_type : str
        Target nomenclature for matching

    Returns
    -------
    tuple
        (matched_genes, list1_converted, list2_converted)
    """
    if auto_convert:
        # Detect types
        type1 = detect_id_type(list1)
        type2 = detect_id_type(list2)

        print(f"List 1 detected as: {type1}")
        print(f"List 2 detected as: {type2}")

        # Convert both to target type
        if type1 != target_type:
            conv1 = convert_gene_ids(list1, from_type=type1, to_type=target_type)
            list1_standardized = list(conv1.values())
        else:
            list1_standardized = list1

        if type2 != target_type:
            conv2 = convert_gene_ids(list2, from_type=type2, to_type=target_type)
            list2_standardized = list(conv2.values())
        else:
            list2_standardized = list2

    else:
        list1_standardized = list1
        list2_standardized = list2

    # Find matches
    set1 = set(list1_standardized)
    set2 = set(list2_standardized)
    matched = set1 & set2

    return matched, list1_standardized, list2_standardized


class FlexibleGeneIDHandler:
    """
    Handler for flexible gene ID processing in benchmarks

    Examples
    --------
    >>> handler = FlexibleGeneIDHandler()
    >>> handler.add_essential_genes(['2597', '3156'], id_type='entrez_id')
    >>> handler.add_model_genes(['GAPDH', 'HMGCR'], id_type='symbol')
    >>> matched = handler.match_genes()
    """

    def __init__(self, target_nomenclature: str = 'symbol'):
        """
        Initialize handler

        Parameters
        ----------
        target_nomenclature : str
            Target nomenclature for standardization
        """
        self.target_nomenclature = target_nomenclature
        self.gene_lists = {}
        self.converted_lists = {}

    def add_gene_list(
        self,
        name: str,
        genes: List[str],
        id_type: Optional[str] = None
    ):
        """
        Add a gene list

        Parameters
        ----------
        name : str
            Name of the gene list
        genes : List[str]
            Gene IDs
        id_type : str, optional
            ID type (auto-detect if None)
        """
        self.gene_lists[name] = {
            'genes': genes,
            'id_type': id_type or detect_id_type(genes)
        }

        # Convert to target nomenclature
        if self.gene_lists[name]['id_type'] != self.target_nomenclature:
            converted = convert_gene_ids(
                genes,
                from_type=self.gene_lists[name]['id_type'],
                to_type=self.target_nomenclature
            )
            self.converted_lists[name] = list(converted.values())
        else:
            self.converted_lists[name] = genes

    def get_converted_list(self, name: str) -> List[str]:
        """Get converted gene list"""
        return self.converted_lists.get(name, [])

    def match_lists(self, list1_name: str, list2_name: str) -> Set[str]:
        """Match two gene lists"""
        set1 = set(self.get_converted_list(list1_name))
        set2 = set(self.get_converted_list(list2_name))
        return set1 & set2


def prepare_validation_genes(
    essential_genes: Optional[List[str]] = None,
    non_essential_genes: Optional[List[str]] = None,
    essential_id_type: Optional[str] = None,
    non_essential_id_type: Optional[str] = None,
    model_genes: Optional[List[str]] = None,
    model_id_type: Optional[str] = None,
    target_type: str = 'symbol'
) -> Dict[str, List[str]]:
    """
    Prepare gene lists for validation, handling ID conversions

    Parameters
    ----------
    essential_genes : List[str], optional
        Essential genes
    non_essential_genes : List[str], optional
        Non-essential genes
    essential_id_type : str, optional
        Essential genes ID type
    non_essential_id_type : str, optional
        Non-essential genes ID type
    model_genes : List[str], optional
        Genes in the model
    model_id_type : str, optional
        Model genes ID type
    target_type : str
        Target nomenclature

    Returns
    -------
    Dict[str, List[str]]
        Standardized gene lists
    """
    result = {}

    if essential_genes:
        if essential_id_type and essential_id_type != target_type:
            conv = convert_gene_ids(essential_genes, from_type=essential_id_type, to_type=target_type)
            result['essential'] = list(conv.values())
        else:
            result['essential'] = essential_genes

    if non_essential_genes:
        if non_essential_id_type and non_essential_id_type != target_type:
            conv = convert_gene_ids(non_essential_genes, from_type=non_essential_id_type, to_type=target_type)
            result['non_essential'] = list(conv.values())
        else:
            result['non_essential'] = non_essential_genes

    if model_genes:
        if model_id_type and model_id_type != target_type:
            conv = convert_gene_ids(model_genes, from_type=model_id_type, to_type=target_type)
            result['model'] = list(conv.values())
        else:
            result['model'] = model_genes

    return result


# Preset mappings for common ID types
ID_TYPE_ALIASES = {
    'entrez': 'entrez_id',
    'ncbi': 'entrez_id',
    'ensembl': 'ensembl_gene_id',
    'ensg': 'ensembl_gene_id',
    'hgnc': 'hgnc_id',
    'gene_symbol': 'symbol',
    'gene_name': 'symbol'
}


def normalize_id_type(id_type: str) -> str:
    """Normalize ID type string"""
    id_type_lower = id_type.lower().strip()
    return ID_TYPE_ALIASES.get(id_type_lower, id_type_lower)
