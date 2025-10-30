# Troppo 확장성 가이드 (Extensibility Guide)

이 문서는 Troppo 패키지를 확장하는 방법을 설명합니다.

This document explains how to extend the Troppo package.

---

## 📋 목차 (Table of Contents)

- [개요](#개요-overview)
- [방법론 레지스트리 시스템](#방법론-레지스트리-시스템)
- [새로운 방법론 추가하기](#새로운-방법론-추가하기)
- [벤치마크 프레임워크](#벤치마크-프레임워크)
- [성능 비교 및 평가](#성능-비교-및-평가)
- [예제](#예제)

---

## 개요 (Overview)

Troppo는 새로운 오믹스 통합 방법론을 쉽게 추가하고 기존 방법론과 비교할 수 있도록 설계되었습니다.

Troppo is designed to make it easy to add new omics integration methods and compare them with existing methods.

### 주요 기능 (Key Features)

1. **플러그인 시스템 (Plugin System)**
   - 새로운 방법론을 등록하고 관리하는 레지스트리
   - 데코레이터 기반 간편 등록
   - 메타데이터 관리 (속도, 메모리, 복잡도 등)

2. **벤치마크 프레임워크 (Benchmark Framework)**
   - 여러 방법론 동시 실행 및 비교
   - 자동 성능 지표 수집
   - 생물학적 검증 지표
   - 결과 시각화 및 레포트 생성

3. **확장성 (Extensibility)**
   - 표준 인터페이스를 따르는 모든 방법론 지원
   - Hook 시스템으로 커스텀 동작 추가
   - 유연한 설정 관리

---

## 방법론 레지스트리 시스템

### 레지스트리란?

`MethodRegistry`는 모든 사용 가능한 오믹스 통합 방법론을 중앙에서 관리하는 시스템입니다.

The `MethodRegistry` is a centralized system for managing all available omics integration methods.

### 레지스트리 사용하기

```python
from troppo.methods.registry import MethodRegistry

# 등록된 모든 방법론 보기
methods = MethodRegistry.list_methods()
print(methods)  # ['gimme', 'tinit', 'imat', 'fastcore', ...]

# 특정 방법론 정보 가져오기
metadata = MethodRegistry.get_method('gimme')
print(metadata.description)
print(metadata.speed)
print(metadata.complexity)

# 모든 메타데이터 출력
MethodRegistry.print_registry()
```

### 필터링

```python
# Continuous scores를 요구하는 방법론만
continuous_methods = MethodRegistry.list_methods(requires_continuous_scores=True)

# Reconstruction 카테고리만
reconstruction_methods = MethodRegistry.list_methods(category='reconstruction')

# 모든 메타데이터 가져오기
all_metadata = MethodRegistry.get_all_metadata()
```

---

## 새로운 방법론 추가하기

### Step 1: Properties 클래스 정의

```python
from cobamp.utilities.property_management import PropertyDictionary

class MyMethodProperties(PropertyDictionary):
    """
    Properties for MyMethod
    """

    def __init__(
        self,
        expression_scores: dict,
        threshold: float = 0.5,
        solver: str = 'CPLEX',
        reaction_ids: list = None,
        metabolite_ids: list = None,
        **kwargs
    ):
        # Define mandatory properties
        mandatory = {
            'expression_scores': dict,
            'solver': str
        }

        # Define optional properties
        optional = {
            'threshold': float,
            'reaction_ids': list,
            'metabolite_ids': list
        }

        super().__init__(mandatory, optional)

        # Set properties
        self.expression_scores = expression_scores
        self.threshold = threshold
        self.solver = solver
        self.reaction_ids = reaction_ids or []
        self.metabolite_ids = metabolite_ids or []

        # Store additional kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)
```

### Step 2: Method 클래스 구현

```python
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm
import numpy as np

class MyMethod(ContextSpecificModelReconstructionAlgorithm):
    """
    My custom omics integration method

    This method does [describe your algorithm here]...
    """

    def __init__(self, S, lb, ub, properties: MyMethodProperties):
        """
        Initialize the method

        Parameters
        ----------
        S : array-like
            Stoichiometric matrix
        lb : array-like
            Lower bounds
        ub : array-like
            Upper bounds
        properties : MyMethodProperties
            Method properties
        """
        self.S = S
        self.lb = lb
        self.ub = ub
        self.properties = properties
        self.solver = properties.solver

    def run(self):
        """
        Execute the method

        Returns
        -------
        list
            List of selected reaction indices
        """
        # Your algorithm implementation here

        # Example: select reactions with scores above threshold
        scores = self.properties.expression_scores
        threshold = self.properties.threshold

        selected = []
        for idx, rxn_id in enumerate(self.properties.reaction_ids):
            if rxn_id in scores and scores[rxn_id] >= threshold:
                selected.append(idx)

        return selected

    @property
    def properties_class(self):
        """Return the properties class"""
        return MyMethodProperties
```

### Step 3: 레지스트리에 등록

#### 방법 1: 직접 등록

```python
from troppo.methods.registry import MethodRegistry

MethodRegistry.register(
    name='my_method',
    method_class=MyMethod,
    properties_class=MyMethodProperties,
    description='My custom integration method',
    category='reconstruction',
    requires_continuous_scores=True,
    speed='fast',  # fast, medium, slow
    memory='low',  # low, medium, high
    complexity='low',  # low, medium, high
    reference='Your Name et al. (2024)'
)
```

#### 방법 2: 데코레이터 사용

```python
from troppo.methods.registry import register_method

@register_method(
    name='my_method',
    description='My custom integration method',
    requires_continuous_scores=True,
    speed='fast'
)
class MyMethod(ContextSpecificModelReconstructionAlgorithm):
    properties_class = MyMethodProperties  # Required!

    def __init__(self, S, lb, ub, properties):
        # ...

    def run(self):
        # ...
```

### Step 4: 사용하기

```python
# Create properties
properties = MyMethodProperties(
    expression_scores={'rxn1': 0.8, 'rxn2': 0.6},
    threshold=0.5,
    solver='CPLEX',
    reaction_ids=['rxn1', 'rxn2', 'rxn3']
)

# Create method instance using registry
method = MethodRegistry.create_method(
    'my_method',
    S=S_matrix,
    lb=lower_bounds,
    ub=upper_bounds,
    properties=properties
)

# Run
result = method.run()
print(f"Selected {len(result)} reactions")
```

---

## 벤치마크 프레임워크

### 기본 사용법

```python
from troppo.benchmark import BenchmarkRunner

# Create runner
runner = BenchmarkRunner(
    model_wrapper=model_wrapper,
    data_map=data_map,
    methods=['gimme', 'tinit', 'imat', 'my_method'],  # Include your method!
    biomass_reaction='biomass_reaction_id',
    verbose=True
)

# Run benchmark
comparison = runner.run_benchmark(
    validate_biology=True,
    validate_consistency=True
)

# View results
print(comparison.get_summary_dataframe())
```

### 빠른 벤치마크

```python
from troppo.benchmark import quick_benchmark

# One-line benchmark
summary = quick_benchmark(
    model_wrapper=model_wrapper,
    data_map=data_map,
    methods=['gimme', 'my_method']
)

print(summary)
```

### 커스텀 설정으로 벤치마크

```python
# Define method-specific configurations
method_configs = {
    'my_method': {
        'threshold': 0.7,
        'custom_param': 'value'
    },
    'gimme': {
        'obj_frac': 0.9
    }
}

# Run with custom configs
comparison = runner.run_benchmark(
    method_configs=method_configs,
    validate_biology=True
)
```

---

## 성능 비교 및 평가

### 수집되는 지표 (Collected Metrics)

#### 1. 성능 지표 (Performance Metrics)
- **execution_time**: 실행 시간 (초)
- **peak_memory**: 최대 메모리 사용량 (MB)

#### 2. 모델 크기 지표 (Model Size Metrics)
- **num_reactions_selected**: 선택된 반응 수
- **percentage_retained**: 유지된 반응 비율

#### 3. 생물학적 검증 지표 (Biological Validation Metrics)
- **biomass_flux**: Biomass 생성 플럭스
- **growth_rate**: 성장 속도
- **task_completion**: 대사 task 완료 여부
- **task_completion_rate**: Task 완료율

#### 4. 네트워크 품질 지표 (Network Quality Metrics)
- **network_consistency**: 네트워크 일관성
- **blocked_reactions**: 막힌 반응 수
- **num_blocked_reactions**: 막힌 반응 개수

### 결과 분석

```python
# Get summary table
summary_df = comparison.get_summary_dataframe()

# Get overlap matrix (Jaccard similarity)
overlap_matrix = comparison.get_overlap_matrix()

# Find common reactions
common_reactions = comparison.get_common_reactions()

# Find unique reactions for a method
unique = comparison.get_unique_reactions('my_method')

# Rank methods by metric
rankings = comparison.rank_methods('execution_time')
for rank, (method, time) in enumerate(rankings, 1):
    print(f"{rank}. {method}: {time:.2f}s")
```

### 시각화

```python
from troppo.benchmark.visualization import (
    plot_performance_comparison,
    plot_overlap_heatmap,
    plot_pareto_front,
    plot_radar_chart,
    create_comparison_report
)

# Performance comparison
plot_performance_comparison(comparison)

# Overlap heatmap
plot_overlap_heatmap(comparison)

# Pareto front (trade-off analysis)
plot_pareto_front(
    comparison,
    x_metric='execution_time',
    y_metric='percentage_retained'
)

# Radar chart
plot_radar_chart(comparison, normalize=True)

# Generate comprehensive report
create_comparison_report(
    comparison,
    output_dir='benchmark_results',
    include_plots=True
)
```

### 결과 저장 및 로드

```python
# Save to JSON
comparison.save_to_json('results.json')

# Load from JSON
from troppo.benchmark import BenchmarkComparison
loaded = BenchmarkComparison.load_from_json('results.json')
```

---

## 예제

### 완전한 예제: 새로운 방법론 추가 및 벤치마크

```python
# 1. Define Properties
class SimpleThresholdProperties(PropertyDictionary):
    def __init__(self, scores: dict, threshold: float, **kwargs):
        mandatory = {'scores': dict, 'threshold': float}
        optional = {'reaction_ids': list}
        super().__init__(mandatory, optional)

        self.scores = scores
        self.threshold = threshold
        self.reaction_ids = kwargs.get('reaction_ids', [])
        self.solver = kwargs.get('solver', 'CPLEX')

# 2. Implement Method
class SimpleThreshold(ContextSpecificModelReconstructionAlgorithm):
    """Simple threshold-based selection"""

    properties_class = SimpleThresholdProperties

    def __init__(self, S, lb, ub, properties):
        self.S = S
        self.lb = lb
        self.ub = ub
        self.properties = properties

    def run(self):
        selected = []
        for idx, rxn_id in enumerate(self.properties.reaction_ids):
            score = self.properties.scores.get(rxn_id, 0)
            if score >= self.properties.threshold:
                selected.append(idx)
        return selected

# 3. Register
MethodRegistry.register(
    'simple_threshold',
    SimpleThreshold,
    SimpleThresholdProperties,
    description='Simple threshold-based method',
    requires_continuous_scores=True,
    speed='fast',
    memory='low'
)

# 4. Benchmark
from troppo.benchmark import quick_benchmark

summary = quick_benchmark(
    model_wrapper=model_wrapper,
    data_map=data_map,
    methods=['gimme', 'tinit', 'simple_threshold']
)

print(summary)
```

---

## 고급 기능 (Advanced Features)

### Hook 시스템

```python
from troppo.methods.registry import MethodRegistry

# Add a pre-registration hook
def my_pre_register_hook(name, method_class, properties_class, kwargs):
    print(f"Registering method: {name}")

MethodRegistry.add_hook('pre_register', my_pre_register_hook)

# Add a post-creation hook
def my_post_create_hook(name, instance):
    print(f"Created instance of: {name}")

MethodRegistry.add_hook('post_create', my_post_create_hook)
```

### 커스텀 검증 지표

```python
class MyBenchmarkRunner(BenchmarkRunner):
    def _validate_biology(self, result, selected_indices):
        super()._validate_biology(result, selected_indices)

        # Add custom validation
        # ... your custom logic here ...

        result.custom_metric = my_custom_calculation()
```

---

## 모범 사례 (Best Practices)

### 1. 명확한 문서화

```python
class MyMethod(ContextSpecificModelReconstructionAlgorithm):
    """
    My Method Name

    Brief description of what the method does.

    Algorithm:
    1. Step 1 description
    2. Step 2 description
    3. ...

    Parameters
    ----------
    S : array-like
        Stoichiometric matrix
    ...

    References
    ----------
    .. [1] Author et al. (Year). Title. Journal.

    Examples
    --------
    >>> properties = MyMethodProperties(...)
    >>> method = MyMethod(S, lb, ub, properties)
    >>> result = method.run()
    """
```

### 2. 에러 처리

```python
def run(self):
    try:
        # Your algorithm
        result = self._do_something()
        return result
    except Exception as e:
        raise RuntimeError(f"MyMethod failed: {str(e)}")
```

### 3. 유효성 검사

```python
def __init__(self, S, lb, ub, properties):
    # Validate inputs
    if S.shape[0] != len(lb):
        raise ValueError("Dimension mismatch between S and lb")

    if properties.threshold < 0 or properties.threshold > 1:
        raise ValueError("Threshold must be between 0 and 1")

    self.S = S
    self.lb = lb
    self.ub = ub
    self.properties = properties
```

### 4. 테스트

```python
import unittest

class TestMyMethod(unittest.TestCase):
    def test_basic_functionality(self):
        # Create test data
        S = np.array([[1, -1, 0], [0, 1, -1]])
        lb = np.array([0, 0, 0])
        ub = np.array([10, 10, 10])

        properties = MyMethodProperties(
            scores={'rxn1': 0.8, 'rxn2': 0.6},
            threshold=0.5,
            reaction_ids=['rxn1', 'rxn2', 'rxn3']
        )

        method = MyMethod(S, lb, ub, properties)
        result = method.run()

        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
```

---

## 문제 해결 (Troubleshooting)

### 방법론이 레지스트리에 나타나지 않음

```python
# Check if registered
if 'my_method' in MethodRegistry.list_methods():
    print("Method is registered")
else:
    print("Method NOT registered")
    # Re-register
    MethodRegistry.register(...)
```

### Properties 클래스를 찾을 수 없음

데코레이터를 사용할 때는 반드시 `properties_class` 속성을 설정하세요:

```python
@register_method('my_method', ...)
class MyMethod(...):
    properties_class = MyMethodProperties  # Required!
```

### 벤치마크 실행 오류

```python
# Enable verbose mode for debugging
runner = BenchmarkRunner(..., verbose=True)

# Catch errors
try:
    comparison = runner.run_benchmark()
except Exception as e:
    print(f"Benchmark failed: {str(e)}")
    import traceback
    traceback.print_exc()
```

---

## 리소스 (Resources)

### 파일 구조

```
troppo/
├── src/troppo/
│   ├── methods/
│   │   ├── base.py                 # Base classes
│   │   ├── registry.py             # Method registry
│   │   └── reconstruction/
│   │       ├── gimme.py
│   │       ├── tinit.py
│   │       └── ...
│   └── benchmark/
│       ├── __init__.py
│       ├── framework.py            # Benchmark framework
│       └── visualization.py        # Visualization tools
├── examples/
│   └── custom_method_example.py    # Complete example
├── tests/
│   ├── Troppo_tutorial_benchmark.ipynb
│   └── ...
└── docs/
    └── EXTENSIBILITY_GUIDE.md      # This file
```

### 관련 문서

- **OMICS_INTEGRATION_GUIDE.md**: 오믹스 통합 사용 가이드
- **examples/custom_method_example.py**: 완전한 커스텀 방법론 예제
- **tests/Troppo_tutorial_benchmark.ipynb**: 벤치마크 튜토리얼

### 온라인 리소스

- GitHub: https://github.com/BioSystemsUM/troppo
- Documentation: http://troppo-bisbi.readthedocs.io/
- Issues: https://github.com/BioSystemsUM/troppo/issues

---

## 기여하기 (Contributing)

새로운 방법론을 개발했다면 Troppo 프로젝트에 기여해주세요!

If you've developed a new method, consider contributing to the Troppo project!

### Pull Request 절차

1. Fork the repository
2. Create a new branch (`git checkout -b feature/my-method`)
3. Implement your method following the guidelines
4. Add tests
5. Update documentation
6. Submit a pull request

### 요구사항

- [ ] Properties 클래스 구현
- [ ] Method 클래스 구현
- [ ] 레지스트리에 등록
- [ ] 문서화 (docstrings)
- [ ] 단위 테스트
- [ ] 사용 예제
- [ ] 벤치마크 테스트

---

## 라이센스 (License)

Troppo is released under the GNU General Public License v3.0

---

**Last Updated:** 2025-10-25
**Version:** 1.0
**Authors:** Troppo Development Team
