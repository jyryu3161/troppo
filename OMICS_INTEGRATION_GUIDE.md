# Troppo 오믹스 통합 가이드 (Omics Integration Guide)

이 가이드는 Troppo 패키지를 사용하여 오믹스 데이터를 대사 모델에 통합하는 방법을 설명합니다.

This guide explains how to integrate omics data into metabolic models using the Troppo package.

## 📋 목차 (Table of Contents)

- [소개 (Introduction)](#소개-introduction)
- [설치 (Installation)](#설치-installation)
- [사용 방법 (Usage)](#사용-방법-usage)
  - [Jupyter 노트북 튜토리얼](#jupyter-노트북-튜토리얼)
  - [Shell Script 실행](#shell-script-실행)
- [통합 방법론 (Integration Methods)](#통합-방법론-integration-methods)
- [예제 (Examples)](#예제-examples)
- [유전자 ID 명명법 처리 (Handling Gene ID Nomenclatures)](#유전자-id-명명법-처리-handling-gene-id-nomenclatures)
- [문제 해결 (Troubleshooting)](#문제-해결-troubleshooting)

---

## 소개 (Introduction)

Troppo는 오믹스 데이터(transcriptomics, proteomics 등)를 대사 모델에 통합하여 조직 특이적(tissue-specific) 또는 조건 특이적(condition-specific) 대사 네트워크를 재구성하는 파이썬 패키지입니다.

Troppo is a Python package for integrating omics data (transcriptomics, proteomics, etc.) into metabolic models to reconstruct tissue-specific or condition-specific metabolic networks.

### 구현된 통합 방법론 (Implemented Integration Methods)

1. **GIMME** (Gene Inactivity Moderated by Metabolism and Expression)
2. **tINIT** (Task-driven Integrative Network Inference for Tissues)
3. **iMAT** (Integrative Metabolic Analysis Tool)
4. **FastCORE** (Fast Consistency-based Reconstruction)
5. **CORDA** (Cost Optimization Reaction Dependency Assessment)

---

## 설치 (Installation)

### 필수 요구사항 (Prerequisites)

- Python 3.7 이상
- CPLEX, Gurobi, 또는 GLPK 같은 LP/MILP 솔버

### Troppo 설치

```bash
# PyPI로부터 설치 (Install from PyPI)
pip install troppo

# 또는 GitHub에서 최신 버전 설치 (Or install latest version from GitHub)
pip install git+https://github.com/BioSystemsUM/troppo.git
```

### 의존성 패키지 (Dependencies)

```bash
pip install cobra framed pandas numpy
```

---

## 사용 방법 (Usage)

### Jupyter 노트북 튜토리얼

포괄적인 튜토리얼 노트북이 제공됩니다:

```bash
jupyter notebook tests/Troppo_tutorial_omics_integration.ipynb
```

이 노트북에는 다음이 포함되어 있습니다:
- 데이터 로딩 및 전처리
- 각 통합 방법론의 상세한 설명
- 단계별 실행 예제
- 결과 비교 및 분석

This notebook includes:
- Data loading and preprocessing
- Detailed explanation of each integration method
- Step-by-step execution examples
- Result comparison and analysis

### Shell Script 실행

자동화된 파이프라인 실행:

```bash
# 모든 방법론 실행 (Run all methods)
./run_omics_integration.sh

# 특정 방법론만 실행 (Run specific method)
./run_omics_integration.sh gimme

# 커스텀 데이터로 실행 (Run with custom data)
./run_omics_integration.sh all path/to/model.xml path/to/omics_data.csv
```

#### Shell Script 옵션

```
Usage: ./run_omics_integration.sh [method] [model_path] [omics_data_path]

Arguments:
  method          - Integration method: gimme, tinit, imat, fastcore, or all (default: all)
  model_path      - Path to SBML model file (default: tests/data/HumanGEM_Consistent_COVID19_HAM.xml)
  omics_data_path - Path to omics data CSV file (default: tests/data/Desai-GTEx_ensembl.csv)
```

---

## 통합 방법론 (Integration Methods)

### 1. GIMME

**특징:**
- 연속적인 유전자 발현 값 사용
- 발현이 낮은 반응의 플럭스를 최소화하여 모델 재구성
- 빠른 실행 속도

**적용 사례:**
- Quantitative gene expression data (RNA-seq, microarray)
- 발현 값이 연속적인 경우

**Python 예제:**
```python
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

properties = GIMMEProperties(
    exp_vector=expression_scores,
    obj_frac=0.8,
    objectives=[{biomass_idx: 1}],
    solver='CPLEX'
)

gimme = GIMME(S=S_matrix, lb=lower_bounds, ub=upper_bounds, properties=properties)
result = gimme.run()
```

### 2. tINIT

**특징:**
- Task 기반 접근 방식
- 조직 특이적 대사 기능(tasks) 보존
- 높은 발현을 가진 반응 우선 선택

**적용 사례:**
- Tissue-specific model reconstruction
- 특정 대사 기능을 보존해야 하는 경우

**Python 예제:**
```python
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

properties = tINITProperties(
    core=high_expression_reactions,
    solver='CPLEX'
)

tinit = tINIT(S=S_matrix, lb=lower_bounds, ub=upper_bounds, properties=properties)
result = tinit.run()
```

### 3. iMAT

**특징:**
- 발현 데이터를 high/moderate/low로 분류
- 발현 패턴과 플럭스 분포의 일치성 최대화
- 질적(qualitative) 발현 데이터에 적합

**적용 사례:**
- Binary or categorical expression data
- 발현 수준을 명확히 구분할 수 있는 경우

**Python 예제:**
```python
from troppo.methods.reconstruction.imat import IMAT, IMATProperties

properties = IMATProperties(
    rh=highly_expressed_reactions,
    rl=lowly_expressed_reactions,
    epsilon=1.0,
    solver='CPLEX'
)

imat = IMAT(S=S_matrix, lb=lower_bounds, ub=upper_bounds, properties=properties)
result = imat.run()
```

### 4. FastCORE

**특징:**
- Core 반응 세트를 보존하는 최소 일관성 네트워크 생성
- 매우 빠른 실행 속도
- 단순하고 효율적인 알고리즘

**적용 사례:**
- 빠른 프로토타이핑이 필요한 경우
- Core 반응 세트가 명확한 경우

**Python 예제:**
```python
from troppo.methods.reconstruction.fastcore import FASTCORE, FASTCOREProperties

properties = FASTCOREProperties(
    core=core_reactions,
    solver='CPLEX'
)

fastcore = FASTCORE(S=S_matrix, lb=lower_bounds, ub=upper_bounds, properties=properties)
result = fastcore.run()
```

---

## 예제 (Examples)

### 기본 워크플로우 (Basic Workflow)

```python
import pandas as pd
import cobra
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper

# 1. 모델 및 데이터 로드
model = cobra.io.read_sbml_model('model.xml')
omics_data = pd.read_csv('expression_data.csv', index_col=0)

# 2. 오믹스 컨테이너 생성
reader = TabularReader(
    path_or_df=omics_data,
    nomenclature='ensemble_gene_id',
    omics_type='transcriptomics'
)
omics_container = reader.to_containers()[0]

# 3. 모델 래퍼 생성
model_wrapper = ModelBasedWrapper(model=model, ttg_ratio=9999)

# 4. 유전자-반응 매핑
data_map = omics_container.get_integrated_data_map(
    model_reader=model_wrapper.model_reader,
    and_func=min,
    or_func=sum
)

# 5. 통합 방법론 적용 (예: GIMME)
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

continuous_integration = ContinuousScoreIntegrationStrategy(
    score_apply=lambda x: {k: 0 if v is None else v for k, v in x.items()}
)
scores = continuous_integration.integrate(data_map=data_map)

idx_biomass = model_wrapper.model_reader.r_ids.index('biomass_reaction')
properties = GIMMEProperties(
    exp_vector=[v for k, v in scores.items()],
    obj_frac=0.8,
    objectives=[{idx_biomass: 1}],
    solver='CPLEX'
)

gimme = GIMME(
    S=model_wrapper.S,
    lb=model_wrapper.lb,
    ub=model_wrapper.ub,
    properties=properties
)

result = gimme.run()
print(f"Selected {len(result)} reactions")
```

### 여러 방법론 비교 (Comparing Multiple Methods)

```python
# 각 방법론 실행
methods_results = {
    'GIMME': run_gimme(data_map, model_wrapper),
    'tINIT': run_tinit(data_map, model_wrapper),
    'iMAT': run_imat(data_map, model_wrapper),
    'FastCORE': run_fastcore(data_map, model_wrapper)
}

# 결과 비교
for method_name, result in methods_results.items():
    num_reactions = len(result)
    percentage = num_reactions / len(model.reactions) * 100
    print(f"{method_name}: {num_reactions} reactions ({percentage:.2f}%)")

# 공통 반응 찾기
common_reactions = set.intersection(*[set(r) for r in methods_results.values()])
print(f"Common reactions across all methods: {len(common_reactions)}")
```

---

## 워크플로우 다이어그램 (Workflow Diagram)

```
┌─────────────────────────────────────────────────────────────────┐
│                    Omics Integration Pipeline                    │
└─────────────────────────────────────────────────────────────────┘

1. Data Loading
   ├── SBML Model (.xml)
   └── Omics Data (.csv, .tsv)
         ↓
2. Preprocessing
   ├── OmicsContainer creation
   ├── ModelBasedWrapper creation
   └── ID matching & validation
         ↓
3. Gene-Reaction Mapping
   ├── Apply GPR rules (AND/OR logic)
   ├── Gene scores → Reaction scores
   └── Create OmicsDataMap
         ↓
4. Integration Strategy
   ├── Continuous (GIMME)
   ├── Threshold (tINIT, FastCORE)
   └── Classification (iMAT)
         ↓
5. Reconstruction Algorithm
   ├── GIMME: Minimize low-expression flux
   ├── tINIT: Task-based selection
   ├── iMAT: Expression pattern matching
   └── FastCORE: Core preservation
         ↓
6. Results
   ├── Selected reaction indices
   ├── Tissue-specific model
   └── Analysis & validation
```

---

## 데이터 형식 (Data Formats)

### 오믹스 데이터 CSV 형식

샘플이 행(rows), 유전자가 열(columns)인 경우:

```
,ENSG00000000001,ENSG00000000002,ENSG00000000003,...
Sample1,10.5,20.3,5.7,...
Sample2,12.1,18.9,6.2,...
Sample3,9.8,22.1,4.9,...
```

유전자가 행(rows), 샘플이 열(columns)인 경우:

```
Gene,Sample1,Sample2,Sample3,...
ENSG00000000001,10.5,12.1,9.8,...
ENSG00000000002,20.3,18.9,22.1,...
ENSG00000000003,5.7,6.2,4.9,...
```

### 대사 모델

- **형식**: SBML (.xml)
- **요구사항**:
  - GPR (Gene-Protein-Reaction) rules 포함
  - 유전자 ID가 오믹스 데이터와 일치하는 nomenclature 사용

---

## 유전자 ID 명명법 처리 (Handling Gene ID Nomenclatures)

### 다양한 ID 타입 지원 (Support for Multiple ID Types)

Troppo는 발현 데이터와 대사 모델 간 유전자 ID 불일치를 자동으로 처리합니다.

Troppo automatically handles gene ID mismatches between expression data and metabolic models.

#### 지원되는 ID 타입 (Supported ID Types)

- **Entrez ID**: NCBI Entrez Gene IDs (예: '2597', '3156')
- **Ensembl Gene ID**: Ensembl Gene IDs (예: 'ENSG00000111640')
- **Symbol**: HGNC Gene Symbols (예: 'GAPDH', 'HMGCR')
- **HGNC ID**: HGNC IDs (예: 'HGNC:4141')
- **UniProt ID**: UniProt IDs

### 방법 1: 딕셔너리에서 직접 데이터 맵 생성 (Direct Dictionary to Data Map)

가장 간단한 방법:

```python
from troppo.omics import create_data_map_from_dict

# 발현 데이터 (Entrez ID 사용)
expression_data = {
    '2597': 100.5,   # GAPDH (Entrez ID)
    '3156': 85.2,    # HMGCR (Entrez ID)
    '5230': 95.8,    # PGK1 (Entrez ID)
}

# 자동 ID 변환을 통한 데이터 맵 생성
data_map = create_data_map_from_dict(
    expression_data=expression_data,
    model_wrapper=model_wrapper,
    gene_id_type='entrez_id',  # ID 타입 명시 (또는 자동 감지)
    auto_convert=True,          # 자동 변환 활성화
    verbose=True
)
```

### 방법 2: OmicsContainer 사용 (Using OmicsContainer)

더 세밀한 제어가 필요한 경우:

```python
from troppo.omics import OmicsContainer, create_compatible_data_map

# Entrez ID로 OmicsContainer 생성
omics_container = OmicsContainer(
    omicstype='transcriptomics',
    condition='sample1',
    data=expression_data,
    nomenclature='entrez_id'
)

# 자동 ID 변환을 통한 호환 데이터 맵 생성
data_map = create_compatible_data_map(
    omics_container=omics_container,
    model_wrapper=model_wrapper,
    auto_convert=True,
    verbose=True
)
```

### 방법 3: 벤치마크에서 혼합 ID 타입 사용 (Mixed ID Types in Benchmark)

발현 데이터와 검증 데이터에 서로 다른 ID 타입 사용:

```python
from troppo.omics import create_data_map_from_dict
from troppo.benchmark import BenchmarkRunner

# 발현 데이터: Entrez IDs
expression_data = {'2597': 100.5, '3156': 85.2}

# 데이터 맵 생성 (Entrez ID 자동 변환)
data_map = create_data_map_from_dict(
    expression_data,
    model_wrapper,
    gene_id_type='entrez_id',
    auto_convert=True
)

# 벤치마크: 필수 유전자는 Ensembl ID, 비필수 유전자는 Symbol
runner = BenchmarkRunner(
    model_wrapper=model_wrapper,
    data_map=data_map,
    methods=['gimme', 'fastcore'],
    # 필수 유전자 (Ensembl IDs)
    essential_genes=['ENSG00000111640', 'ENSG00000102144'],
    essential_genes_id_type='ensembl_gene_id',
    # 비필수 유전자 (Symbols)
    non_essential_genes=['HMGCR', 'TP53'],
    non_essential_genes_id_type='symbol'
)

comparison = runner.run_benchmark(validate_essentiality=True)
```

### 자동 ID 감지 (Automatic ID Detection)

ID 타입을 지정하지 않으면 자동으로 감지됩니다:

```python
from troppo.omics import detect_expression_data_id_type

# ID 타입 자동 감지
id_type = detect_expression_data_id_type(expression_data)
print(f"Detected ID type: {id_type}")

# 자동 감지를 사용한 데이터 맵 생성 (gene_id_type 생략)
data_map = create_data_map_from_dict(
    expression_data,
    model_wrapper,
    auto_convert=True  # ID 타입은 자동 감지됨
)
```

### 완전한 예제 (Complete Example)

전체 워크플로우 예제는 다음 파일을 참조하세요:

```bash
python examples/expression_data_with_entrez_ids.py
```

이 예제는 다음을 포함합니다:
- Entrez ID 발현 데이터 사용
- 자동 ID 변환
- 혼합 ID 타입을 사용한 벤치마킹
- 여러 명명법 통합 방법

---

## 문제 해결 (Troubleshooting)

### 일반적인 문제 및 해결 방법

#### 1. Solver 관련 오류

**문제:**
```
SolverNotFound: No solver found
```

**해결:**
```bash
# CPLEX 설치 확인
python -c "import cplex; print('CPLEX installed')"

# 또는 GLPK 사용
pip install swiglpk

# Properties에서 solver 변경
properties = GIMMEProperties(..., solver='GLPK')
```

#### 2. ID 매칭 문제

**문제:**
```
Warning: Many genes could not be mapped to reactions
```

**해결:**
1. Nomenclature 확인:
   ```python
   # 오믹스 데이터의 ID 형식 확인
   print(omics_container.get_Nomenclature())

   # 필요시 ID 변환
   omics_container.convertIds('symbol')  # Ensembl → Gene Symbol
   ```

2. 모델의 유전자 ID 확인:
   ```python
   model_genes = [g.id for g in model.genes]
   print(f"Model uses: {model_genes[:5]}")
   ```

#### 3. 메모리 부족

**문제:**
대용량 모델 처리 시 메모리 부족

**해결:**
1. Preprocessing 활성화:
   ```python
   properties = GIMMEProperties(..., preprocess=True)
   ```

2. Core 반응 수 줄이기:
   ```python
   # 더 높은 임계값 사용
   threshold = np.percentile(scores, 90)  # 상위 10%만 선택
   ```

#### 4. 실행 시간이 너무 긴 경우

**해결:**
1. FastCORE 사용 (가장 빠름)
2. Solver 파라미터 조정
3. 모델 단순화 (불필요한 반응 제거)

---

## 성능 비교 (Performance Comparison)

| Method    | Speed  | Memory | Accuracy | Use Case                    |
|-----------|--------|--------|----------|-----------------------------|
| GIMME     | Fast   | Low    | Good     | Continuous expression data  |
| tINIT     | Medium | Medium | High     | Task-based reconstruction   |
| iMAT      | Medium | Medium | Good     | Categorical expression data |
| FastCORE  | Very Fast | Low | Medium   | Quick prototyping          |
| CORDA     | Slow   | High   | High     | Comprehensive analysis      |

---

## 추가 자료 (Additional Resources)

### 공식 문서
- Troppo Documentation: http://troppo-bisbi.readthedocs.io/
- GitHub Repository: https://github.com/BioSystemsUM/troppo

### 논문 (Publications)

1. **GIMME**: Becker SA, Palsson BO. (2008). *Genome Biology*
2. **tINIT**: Agren et al. (2014). *Molecular Systems Biology*
3. **iMAT**: Shlomi et al. (2008). *Nature Biotechnology*
4. **FastCORE**: Vlassis et al. (2014). *PLoS Computational Biology*

### 관련 패키지
- COBRApy: https://opencobra.github.io/cobrapy/
- FRAMED: https://github.com/cdanielmachado/framed

---

## 라이센스 (License)

Troppo is released under the GNU General Public License v3.0

---

## 기여 (Contributing)

버그 리포트, 기능 요청, 풀 리퀘스트를 환영합니다!

Bug reports, feature requests, and pull requests are welcome!

GitHub Issues: https://github.com/BioSystemsUM/troppo/issues

---

## 인용 (Citation)

이 패키지를 사용하는 경우 다음과 같이 인용해주세요:

If you use this package, please cite:

```
Troppo: A Python package for tissue-specific reconstruction and phenotype prediction
Centre of Biological Engineering, University of Minho
```

---

**Last Updated:** 2025-10-25
**Version:** 1.0
**Authors:** Troppo Development Team
