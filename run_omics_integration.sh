#!/bin/bash

###############################################################################
# Troppo 오믹스 통합 실행 스크립트
# Troppo Omics Integration Runner Script
#
# 이 스크립트는 Troppo 패키지를 사용하여 오믹스 데이터를 대사 모델에 통합하는
# 여러 방법론(GIMME, tINIT, iMAT, FastCORE)을 자동으로 실행합니다.
#
# This script automatically runs multiple omics integration methods
# (GIMME, tINIT, iMAT, FastCORE) using the Troppo package.
#
# Usage:
#   ./run_omics_integration.sh [method] [model_path] [omics_data_path]
#
# Arguments:
#   method          - Integration method: gimme, tinit, imat, fastcore, or all (default: all)
#   model_path      - Path to SBML model file (default: tests/data/HumanGEM_Consistent_COVID19_HAM.xml)
#   omics_data_path - Path to omics data CSV file (default: tests/data/Desai-GTEx_ensembl.csv)
#
# Examples:
#   ./run_omics_integration.sh
#   ./run_omics_integration.sh gimme
#   ./run_omics_integration.sh all custom_model.xml custom_data.csv
###############################################################################

# 색상 정의 (Color definitions)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 기본 설정 (Default settings)
METHOD=${1:-all}
MODEL_PATH=${2:-"tests/data/HumanGEM_Consistent_COVID19_HAM.xml"}
OMICS_PATH=${3:-"tests/data/Desai-GTEx_ensembl.csv"}
OUTPUT_DIR="omics_integration_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 로그 함수 (Logging functions)
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 배너 출력 (Print banner)
print_banner() {
    echo ""
    echo "=========================================================================="
    echo "  Troppo Omics Integration Pipeline"
    echo "  오믹스 통합 파이프라인"
    echo "=========================================================================="
    echo ""
}

# 도움말 출력 (Print help)
print_help() {
    echo "Usage: $0 [method] [model_path] [omics_data_path]"
    echo ""
    echo "Arguments:"
    echo "  method          - Integration method: gimme, tinit, imat, fastcore, or all (default: all)"
    echo "  model_path      - Path to SBML model file"
    echo "  omics_data_path - Path to omics data CSV file"
    echo ""
    echo "Examples:"
    echo "  $0"
    echo "  $0 gimme"
    echo "  $0 all custom_model.xml custom_data.csv"
    echo ""
}

# Python 스크립트 생성 (Generate Python script)
generate_python_script() {
    local method=$1
    local output_file="run_${method}_${TIMESTAMP}.py"

    cat > "${OUTPUT_DIR}/${output_file}" << 'PYTHON_SCRIPT'
import sys
import pandas as pd
import cobra
import re
import numpy as np
import json
from datetime import datetime

# Troppo imports
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper
from troppo.omics.integration import (
    ContinuousScoreIntegrationStrategy,
    ThresholdSelectionIntegrationStrategy,
    DefaultCoreIntegrationStrategy
)

# Reconstruction methods
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
from troppo.methods.reconstruction.imat import IMAT, IMATProperties
from troppo.methods.reconstruction.fastcore import FASTCORE, FASTCOREProperties

def main():
    # Command line arguments
    method = sys.argv[1] if len(sys.argv) > 1 else 'all'
    model_path = sys.argv[2] if len(sys.argv) > 2 else 'tests/data/HumanGEM_Consistent_COVID19_HAM.xml'
    omics_path = sys.argv[3] if len(sys.argv) > 3 else 'tests/data/Desai-GTEx_ensembl.csv'
    output_dir = sys.argv[4] if len(sys.argv) > 4 else 'omics_integration_results'

    print("=" * 80)
    print("Troppo Omics Integration")
    print("=" * 80)
    print(f"Method: {method}")
    print(f"Model: {model_path}")
    print(f"Omics Data: {omics_path}")
    print(f"Output Directory: {output_dir}")
    print("=" * 80)

    # GPR parsing function
    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    # Load model
    print("\n[1/6] Loading metabolic model...")
    model = cobra.io.read_sbml_model(model_path)
    print(f"  ✓ Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    # Load omics data
    print("\n[2/6] Loading omics data...")
    omics_data = pd.read_csv(omics_path, index_col=0)
    print(f"  ✓ Omics data loaded: {omics_data.shape[0]} samples, {omics_data.shape[1]} genes")

    # Create omics container
    print("\n[3/6] Creating omics container...")
    reader = TabularReader(
        path_or_df=omics_data,
        nomenclature='ensemble_gene_id',
        omics_type='transcriptomics'
    )
    omics_container = reader.to_containers()[0]
    print(f"  ✓ Container created: {len(omics_container.get_Data())} genes")

    # Create model wrapper
    print("\n[4/6] Creating model wrapper...")
    model_wrapper = ModelBasedWrapper(
        model=model,
        ttg_ratio=9999,
        gpr_gene_parse_function=replace_alt_transcripts
    )
    print(f"  ✓ Model wrapper created")

    # Map genes to reactions
    print("\n[5/6] Mapping genes to reactions...")
    data_map = omics_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=min,
        or_func=sum
    )
    print(f"  ✓ Mapping complete: {len(data_map.get_scores())} reaction scores")

    # Run integration methods
    print("\n[6/6] Running integration method(s)...")
    results = {}

    if method == 'all' or method == 'gimme':
        print("\n  Running GIMME...")
        results['GIMME'] = run_gimme(data_map, model_wrapper)
        print(f"    ✓ GIMME complete: {len(results['GIMME'])} reactions selected")

    if method == 'all' or method == 'tinit':
        print("\n  Running tINIT...")
        results['tINIT'] = run_tinit(data_map, model_wrapper)
        print(f"    ✓ tINIT complete: {len(results['tINIT'])} reactions selected")

    if method == 'all' or method == 'imat':
        print("\n  Running iMAT...")
        results['iMAT'] = run_imat(data_map, model_wrapper)
        print(f"    ✓ iMAT complete: {len(results['iMAT'])} reactions selected")

    if method == 'all' or method == 'fastcore':
        print("\n  Running FastCORE...")
        results['FastCORE'] = run_fastcore(data_map, model_wrapper)
        print(f"    ✓ FastCORE complete: {len(results['FastCORE'])} reactions selected")

    # Save results
    print("\n[7/7] Saving results...")
    save_results(results, model_wrapper, output_dir)
    print(f"  ✓ Results saved to {output_dir}/")

    # Print summary
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"{'Method':<15} {'Reactions':<15} {'% of Original':<20}")
    print("-" * 80)
    for method_name, result in results.items():
        num_rxns = len(result)
        percentage = num_rxns / len(model.reactions) * 100
        print(f"{method_name:<15} {num_rxns:<15} {percentage:<20.2f}%")
    print("-" * 80)
    print(f"{'Original Model':<15} {len(model.reactions):<15} {100.0:<20}%")
    print("=" * 80)
    print("\n✓ Pipeline completed successfully!")

def run_gimme(data_map, model_wrapper):
    """Run GIMME algorithm"""
    def score_apply(reaction_map_scores):
        return {k: 0 if v is None else v for k, v in reaction_map_scores.items()}

    continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
    scores = continuous_integration.integrate(data_map=data_map)

    idx_biomass = model_wrapper.model_reader.r_ids.index('biomass_human')

    properties = GIMMEProperties(
        exp_vector=[v for k, v in scores.items()],
        obj_frac=0.8,
        objectives=[{idx_biomass: 1}],
        preprocess=True,
        flux_threshold=0.8,
        solver='CPLEX',
        reaction_ids=model_wrapper.model_reader.r_ids,
        metabolite_ids=model_wrapper.model_reader.m_ids
    )

    gimme = GIMME(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    return gimme.run()

def run_tinit(data_map, model_wrapper):
    """Run tINIT algorithm"""
    threshold_value = np.percentile(
        [v for v in data_map.get_scores().values() if v is not None], 75
    )

    threshold_integration = ThresholdSelectionIntegrationStrategy(thresholds=threshold_value)
    core_reactions = threshold_integration.integrate(data_map=data_map)

    properties = tINITProperties(
        core=core_reactions,
        solver='CPLEX',
        reaction_ids=model_wrapper.model_reader.r_ids,
        metabolite_ids=model_wrapper.model_reader.m_ids
    )

    tinit = tINIT(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    return tinit.run()

def run_imat(data_map, model_wrapper):
    """Run iMAT algorithm"""
    score_values = [v for v in data_map.get_scores().values() if v is not None]
    high_threshold = np.percentile(score_values, 67)
    low_threshold = np.percentile(score_values, 33)

    highly_expressed = set([
        k for k, v in data_map.get_scores().items()
        if v is not None and v >= high_threshold
    ])
    lowly_expressed = set([
        k for k, v in data_map.get_scores().items()
        if v is not None and v <= low_threshold
    ])

    properties = IMATProperties(
        rh=highly_expressed,
        rl=lowly_expressed,
        epsilon=1.0,
        threshold=1e-3,
        solver='CPLEX',
        reaction_ids=model_wrapper.model_reader.r_ids,
        metabolite_ids=model_wrapper.model_reader.m_ids
    )

    imat = IMAT(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    return imat.run()

def run_fastcore(data_map, model_wrapper):
    """Run FastCORE algorithm"""
    core_threshold = np.percentile(
        [v for v in data_map.get_scores().values() if v is not None], 50
    )

    core = set([
        k for k, v in data_map.get_scores().items()
        if v is not None and v >= core_threshold
    ])

    properties = FASTCOREProperties(
        core=core,
        solver='CPLEX',
        reaction_ids=model_wrapper.model_reader.r_ids,
        metabolite_ids=model_wrapper.model_reader.m_ids
    )

    fastcore = FASTCORE(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    return fastcore.run()

def save_results(results, model_wrapper, output_dir):
    """Save results to files"""
    import os
    os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save reaction IDs for each method
    for method_name, result in results.items():
        reaction_ids = [model_wrapper.model_reader.r_ids[i] for i in result]

        output_file = f"{output_dir}/{method_name}_reactions_{timestamp}.txt"
        with open(output_file, 'w') as f:
            f.write(f"# {method_name} Selected Reactions\n")
            f.write(f"# Total: {len(reaction_ids)} reactions\n")
            f.write(f"# Timestamp: {timestamp}\n\n")
            for rxn_id in reaction_ids:
                f.write(f"{rxn_id}\n")

    # Save summary
    summary_file = f"{output_dir}/summary_{timestamp}.json"
    summary = {
        'timestamp': timestamp,
        'methods': {
            method_name: {
                'num_reactions': len(result),
                'reaction_indices': result
            }
            for method_name, result in results.items()
        }
    }

    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

if __name__ == '__main__':
    main()
PYTHON_SCRIPT

    echo "${output_file}"
}

# 메인 실행 함수 (Main execution function)
run_integration() {
    print_banner

    # 입력 검증 (Validate inputs)
    if [[ "$METHOD" != "all" && "$METHOD" != "gimme" && "$METHOD" != "tinit" && "$METHOD" != "imat" && "$METHOD" != "fastcore" ]]; then
        log_error "Invalid method: $METHOD"
        log_info "Valid methods: gimme, tinit, imat, fastcore, all"
        print_help
        exit 1
    fi

    if [[ ! -f "$MODEL_PATH" ]]; then
        log_error "Model file not found: $MODEL_PATH"
        exit 1
    fi

    if [[ ! -f "$OMICS_PATH" ]]; then
        log_error "Omics data file not found: $OMICS_PATH"
        exit 1
    fi

    # 출력 디렉토리 생성 (Create output directory)
    mkdir -p "$OUTPUT_DIR"
    log_info "Output directory created: $OUTPUT_DIR"

    # 설정 출력 (Print configuration)
    log_info "Configuration:"
    echo "  Method: $METHOD"
    echo "  Model: $MODEL_PATH"
    echo "  Omics Data: $OMICS_PATH"
    echo "  Output: $OUTPUT_DIR"
    echo ""

    # Python 스크립트 생성 (Generate Python script)
    log_info "Generating Python execution script..."
    python_script=$(generate_python_script "$METHOD")
    log_success "Python script generated: $OUTPUT_DIR/$python_script"

    # Python 스크립트 실행 (Run Python script)
    log_info "Executing omics integration pipeline..."
    echo ""

    python "$OUTPUT_DIR/$python_script" "$METHOD" "$MODEL_PATH" "$OMICS_PATH" "$OUTPUT_DIR"

    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo ""
        log_success "Pipeline completed successfully!"
        log_info "Results saved in: $OUTPUT_DIR/"
    else
        echo ""
        log_error "Pipeline failed with exit code: $exit_code"
        exit $exit_code
    fi
}

# 도움말 옵션 처리 (Handle help option)
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_help
    exit 0
fi

# 실행 (Execute)
run_integration
