"""
Expression Data Utilities for Gene ID Handling

발현 데이터의 유전자 ID 처리를 위한 유틸리티 함수들

This module provides utilities to handle expression data with different gene ID
nomenclatures, ensuring compatibility with metabolic models.

이 모듈은 다양한 유전자 ID 명명법을 가진 발현 데이터를 처리하여
대사 모델과의 호환성을 보장하는 유틸리티를 제공합니다.
"""

from typing import Dict, Optional, List
import warnings


def detect_expression_data_id_type(expression_data: Dict[str, float]) -> str:
    """
    Detect the gene ID type used in expression data

    발현 데이터에서 사용된 유전자 ID 타입을 감지합니다

    Parameters
    ----------
    expression_data : Dict[str, float]
        Expression data as {gene_id: expression_value}

    Returns
    -------
    str
        Detected ID type ('entrez_id', 'ensembl_gene_id', 'symbol', 'hgnc_id', 'uniprot_ids')

    Examples
    --------
    >>> data = {'2597': 10.5, '3156': 8.2, '5230': 12.1}  # Entrez IDs
    >>> detect_expression_data_id_type(data)
    'entrez_id'

    >>> data = {'GAPDH': 10.5, 'HMGCR': 8.2, 'PGK1': 12.1}  # Symbols
    >>> detect_expression_data_id_type(data)
    'symbol'
    """
    from troppo.benchmark.gene_id_utils import detect_id_type

    gene_ids = list(expression_data.keys())
    return detect_id_type(gene_ids)


def prepare_omics_container_for_model(
    omics_container,
    model_wrapper,
    auto_convert: bool = True,
    target_id_type: Optional[str] = None,
    verbose: bool = True
):
    """
    Prepare OmicsContainer by converting gene IDs to match the model

    모델에 맞춰 OmicsContainer의 유전자 ID를 변환하여 준비합니다

    This function ensures that the gene IDs in the expression data match
    the gene IDs in the metabolic model, enabling proper mapping through
    Gene-Protein-Reaction (GPR) rules.

    이 함수는 발현 데이터의 유전자 ID가 대사 모델의 유전자 ID와 일치하도록
    보장하여 Gene-Protein-Reaction (GPR) 규칙을 통한 올바른 매핑을 가능하게 합니다.

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

    Examples
    --------
    >>> # Expression data with Entrez IDs, model with Symbols
    >>> omics_container = OmicsContainer(
    ...     omicstype='transcriptomics',
    ...     condition='sample1',
    ...     data={'2597': 10.5, '3156': 8.2, '5230': 12.1}  # Entrez IDs
    ... )
    >>> prepared_container = prepare_omics_container_for_model(
    ...     omics_container,
    ...     model_wrapper
    ... )
    >>> # Now container has gene IDs matching the model
    """
    import copy
    from troppo.benchmark.gene_id_utils import detect_id_type, normalize_id_type

    # Make a copy to avoid modifying the original
    container = copy.deepcopy(omics_container)

    # Detect expression data ID type
    expr_data_id_type = detect_expression_data_id_type(container.get_Data())

    if verbose:
        print(f"Expression data ID type detected: {expr_data_id_type}")

    # Determine target ID type
    if target_id_type is None:
        # Try to detect model gene ID type
        try:
            model_genes = [g.id for g in model_wrapper.model_reader.model.genes]
            target_id_type = detect_id_type(model_genes[:100])
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
            # Set nomenclature if not already set
            if container.nomenclature is None:
                container.nomenclature = expr_data_id_type

            # Use OmicsContainer's built-in conversion
            container.convertIds(target_id_type)

            if verbose:
                print(f"Conversion complete: {len(container.get_Data())} genes retained")

            return container

        except Exception as e:
            warnings.warn(
                f"Failed to convert gene IDs from {expr_data_id_type} to {target_id_type}: {e}\n"
                "Returning original container. This may cause issues with GPR mapping."
            )
            return omics_container
    else:
        warnings.warn(
            f"Expression data uses {expr_data_id_type} but model uses {target_id_type}.\n"
            "ID mismatch may cause issues with GPR mapping. Set auto_convert=True to fix."
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
    Create OmicsDataMap with automatic gene ID conversion if needed

    필요시 자동 유전자 ID 변환을 수행하여 OmicsDataMap을 생성합니다

    This is a convenience function that combines ID conversion and data map creation.
    It ensures that expression data gene IDs are compatible with the model before
    creating the reaction score mapping.

    이 함수는 ID 변환과 데이터 맵 생성을 결합한 편의 함수입니다.
    반응 점수 매핑을 생성하기 전에 발현 데이터 유전자 ID가 모델과 호환되는지 확인합니다.

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

    Examples
    --------
    >>> # Expression data with Entrez IDs
    >>> omics_container = OmicsContainer(
    ...     omicstype='transcriptomics',
    ...     condition='sample1',
    ...     data={'2597': 10.5, '3156': 8.2}  # Entrez IDs: GAPDH, HMGCR
    ... )
    >>>
    >>> # Create data map with automatic ID conversion
    >>> data_map = create_compatible_data_map(
    ...     omics_container,
    ...     model_wrapper,
    ...     auto_convert=True,
    ...     verbose=True
    ... )
    >>> # data_map now contains reaction scores based on converted gene IDs
    """
    # Prepare container with proper gene IDs
    prepared_container = prepare_omics_container_for_model(
        omics_container,
        model_wrapper,
        auto_convert=auto_convert,
        verbose=verbose
    )

    # Create data map using the prepared container
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
    Create OmicsDataMap directly from expression data dictionary

    발현 데이터 딕셔너리로부터 직접 OmicsDataMap을 생성합니다

    This is the most convenient function for creating a data map from raw
    expression data. It handles OmicsContainer creation and ID conversion
    automatically.

    이 함수는 원시 발현 데이터로부터 데이터 맵을 생성하는 가장 편리한 함수입니다.
    OmicsContainer 생성과 ID 변환을 자동으로 처리합니다.

    Parameters
    ----------
    expression_data : Dict[str, float]
        Expression data as {gene_id: expression_value}
    model_wrapper : ModelBasedWrapper
        The model wrapper containing the metabolic model
    omicstype : str, optional
        Type of omics data ('transcriptomics' or 'proteomics', default: 'transcriptomics')
    condition : str, optional
        Sample condition name (default: 'sample')
    gene_id_type : str, optional
        Gene ID type used in expression_data. If None, auto-detected.
        Options: 'entrez_id', 'ensembl_gene_id', 'symbol', 'hgnc_id', 'uniprot_ids'
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

    Examples
    --------
    >>> # Expression data with Entrez IDs
    >>> expression_data = {
    ...     '2597': 10.5,   # GAPDH
    ...     '3156': 8.2,    # HMGCR
    ...     '5230': 12.1,   # PGK1
    ...     '5315': 9.8     # PKM
    ... }
    >>>
    >>> # Create data map directly
    >>> data_map = create_data_map_from_dict(
    ...     expression_data,
    ...     model_wrapper,
    ...     gene_id_type='entrez_id',  # Optional: auto-detected if not specified
    ...     verbose=True
    ... )
    """
    from troppo.omics import OmicsContainer
    from troppo.benchmark.gene_id_utils import normalize_id_type

    # Detect gene ID type if not specified
    if gene_id_type is None:
        gene_id_type = detect_expression_data_id_type(expression_data)
        if verbose:
            print(f"Gene ID type auto-detected: {gene_id_type}")
    else:
        gene_id_type = normalize_id_type(gene_id_type)
        if verbose:
            print(f"Gene ID type specified: {gene_id_type}")

    # Create OmicsContainer
    if verbose:
        print(f"Creating OmicsContainer with {len(expression_data)} genes...")

    omics_container = OmicsContainer(
        omicstype=omicstype,
        condition=condition,
        data=expression_data,
        nomenclature=gene_id_type
    )

    # Create compatible data map
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
    'prepare_omics_container_for_model',
    'create_compatible_data_map',
    'create_data_map_from_dict'
]
