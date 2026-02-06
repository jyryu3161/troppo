"""
Expression Data Utilities for Gene ID Handling

This module provides utilities to handle expression data with different gene ID
nomenclatures, ensuring compatibility with metabolic models.
"""

from typing import Dict, Optional, List
import warnings


def detect_expression_data_id_type(expression_data: Dict[str, float]) -> str:
    """
    Detect the gene ID type used in expression data.

    Parameters
    ----------
    expression_data : Dict[str, float]
        Expression data as {gene_id: expression_value}

    Returns
    -------
    str
        Detected ID type ('entrez_id', 'ensembl_gene_id', 'symbol', 'hgnc_id')
    """
    from troppo.io import detect_gene_id_type
    gene_ids = list(expression_data.keys())
    return detect_gene_id_type(gene_ids)


# ID type normalization aliases
_ID_TYPE_ALIASES = {
    'entrez': 'entrez_id',
    'ncbi': 'entrez_id',
    'ensembl': 'ensembl_gene_id',
    'ensg': 'ensembl_gene_id',
    'hgnc': 'hgnc_id',
    'gene_symbol': 'symbol',
    'gene_name': 'symbol',
}


def normalize_id_type(id_type: str) -> str:
    """Normalize ID type string to canonical form."""
    id_type_lower = id_type.lower().strip()
    return _ID_TYPE_ALIASES.get(id_type_lower, id_type_lower)


def prepare_omics_container_for_model(
    omics_container,
    model_wrapper,
    auto_convert: bool = True,
    target_id_type: Optional[str] = None,
    verbose: bool = True
):
    """
    Prepare OmicsContainer by converting gene IDs to match the model.

    Parameters
    ----------
    omics_container : OmicsContainer
        The omics container with expression data
    model_wrapper : ModelBasedWrapper
        The model wrapper containing the metabolic model
    auto_convert : bool, optional
        Automatically convert IDs if they don't match (default: True)
    target_id_type : str, optional
        Target ID type to convert to. If None, uses model's ID type
    verbose : bool, optional
        Print conversion information (default: True)

    Returns
    -------
    OmicsContainer
        The prepared omics container (potentially with converted IDs)
    """
    import copy
    from troppo.io import detect_gene_id_type

    container = copy.deepcopy(omics_container)

    # Detect expression data ID type
    expr_data_id_type = detect_expression_data_id_type(container.get_Data())

    if verbose:
        print(f"Expression data ID type detected: {expr_data_id_type}")

    # Determine target ID type
    if target_id_type is None:
        try:
            model_genes = [g.id for g in model_wrapper.model_reader.model.genes]
            target_id_type = detect_gene_id_type(model_genes[:100])
            if verbose:
                print(f"Model gene ID type detected: {target_id_type}")
        except Exception as e:
            if verbose:
                print(f"Could not detect model gene ID type: {e}")
                print("Using 'symbol' as default target type")
            target_id_type = 'symbol'
    else:
        target_id_type = normalize_id_type(target_id_type)
        if verbose:
            print(f"Target gene ID type specified: {target_id_type}")

    # Check if conversion is needed
    if expr_data_id_type == target_id_type:
        if verbose:
            print(f"Expression data IDs already match model IDs ({target_id_type})")
        return container

    # Convert if auto_convert is enabled
    if auto_convert:
        if verbose:
            print(f"Converting expression data IDs from {expr_data_id_type} to {target_id_type}...")

        try:
            if container.nomenclature is None:
                container.nomenclature = expr_data_id_type
            container.convertIds(target_id_type)
            if verbose:
                print(f"Conversion complete: {len(container.get_Data())} genes retained")
            return container
        except Exception as e:
            warnings.warn(
                f"Failed to convert gene IDs from {expr_data_id_type} to {target_id_type}: {e}\n"
                "Returning original container."
            )
            return omics_container
    else:
        warnings.warn(
            f"Expression data uses {expr_data_id_type} but model uses {target_id_type}.\n"
            "ID mismatch may cause issues. Set auto_convert=True to fix."
        )
        return container


def create_compatible_data_map(
    omics_container,
    model_wrapper,
    auto_convert: bool = True,
    and_func=min,
    or_func=max,
    verbose: bool = True
):
    """
    Create OmicsDataMap with automatic gene ID conversion if needed.

    Parameters
    ----------
    omics_container : OmicsContainer
        The omics container with expression data
    model_wrapper : ModelBasedWrapper
        The model wrapper containing the metabolic model
    auto_convert : bool, optional
        Automatically convert IDs if needed (default: True)
    and_func : callable, optional
        Function to apply for AND in GPR rules (default: min)
    or_func : callable, optional
        Function to apply for OR in GPR rules (default: max)
    verbose : bool, optional
        Print information (default: True)

    Returns
    -------
    OmicsDataMap
        Data map with reaction scores
    """
    prepared_container = prepare_omics_container_for_model(
        omics_container,
        model_wrapper,
        auto_convert=auto_convert,
        verbose=verbose
    )

    if verbose:
        print("Creating data map with GPR-based gene-to-reaction mapping...")

    data_map = prepared_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=and_func,
        or_func=or_func
    )

    if verbose:
        print(f"Data map created: {len(data_map.get_scores())} reactions scored")

    return data_map


def create_data_map_from_dict(
    expression_data: Dict[str, float],
    model_wrapper,
    omicstype: str = 'transcriptomics',
    condition: str = 'sample',
    gene_id_type: Optional[str] = None,
    auto_convert: bool = True,
    and_func=min,
    or_func=max,
    verbose: bool = True
):
    """
    Create OmicsDataMap directly from expression data dictionary.

    Parameters
    ----------
    expression_data : Dict[str, float]
        Expression data as {gene_id: expression_value}
    model_wrapper : ModelBasedWrapper
        The model wrapper containing the metabolic model
    omicstype : str, optional
        Type of omics data (default: 'transcriptomics')
    condition : str, optional
        Sample condition name (default: 'sample')
    gene_id_type : str, optional
        Gene ID type. If None, auto-detected.
    auto_convert : bool, optional
        Automatically convert IDs if needed (default: True)
    and_func : callable, optional
        Function for AND in GPR rules (default: min)
    or_func : callable, optional
        Function for OR in GPR rules (default: max)
    verbose : bool, optional
        Print information (default: True)

    Returns
    -------
    OmicsDataMap
        Data map with reaction scores
    """
    from troppo.omics import OmicsContainer

    if gene_id_type is None:
        gene_id_type = detect_expression_data_id_type(expression_data)
        if verbose:
            print(f"Gene ID type auto-detected: {gene_id_type}")
    else:
        gene_id_type = normalize_id_type(gene_id_type)
        if verbose:
            print(f"Gene ID type specified: {gene_id_type}")

    if verbose:
        print(f"Creating OmicsContainer with {len(expression_data)} genes...")

    omics_container = OmicsContainer(
        omicstype=omicstype,
        condition=condition,
        data=expression_data,
        nomenclature=gene_id_type
    )

    return create_compatible_data_map(
        omics_container,
        model_wrapper,
        auto_convert=auto_convert,
        and_func=and_func,
        or_func=or_func,
        verbose=verbose
    )


__all__ = [
    'detect_expression_data_id_type',
    'normalize_id_type',
    'prepare_omics_container_for_model',
    'create_compatible_data_map',
    'create_data_map_from_dict',
]
