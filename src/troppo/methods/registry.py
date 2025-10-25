"""
Method Registry System for Troppo

This module provides a plugin-like system for registering and managing
omics integration methods, making it easy to add new methods and compare them.

새로운 방법론을 쉽게 등록하고 관리할 수 있는 플러그인 시스템입니다.
"""

from typing import Dict, Type, Callable, Optional, Any, List
import inspect
from abc import ABC, abstractmethod


class MethodMetadata:
    """
    Metadata for a reconstruction method
    재구성 방법에 대한 메타데이터
    """

    def __init__(
        self,
        name: str,
        method_class: Type,
        properties_class: Type,
        description: str = "",
        category: str = "reconstruction",
        requires_continuous_scores: bool = False,
        requires_threshold_scores: bool = False,
        supports_tasks: bool = False,
        complexity: str = "medium",  # low, medium, high
        speed: str = "medium",  # fast, medium, slow
        memory: str = "medium",  # low, medium, high
        reference: str = ""
    ):
        """
        Initialize method metadata

        Parameters
        ----------
        name : str
            Method name (e.g., 'gimme', 'tinit')
        method_class : Type
            The method class (e.g., GIMME, tINIT)
        properties_class : Type
            The properties class (e.g., GIMMEProperties)
        description : str
            Method description
        category : str
            Method category (reconstruction, gapfill, etc.)
        requires_continuous_scores : bool
            Whether method requires continuous scores
        requires_threshold_scores : bool
            Whether method requires threshold-based scores
        supports_tasks : bool
            Whether method supports metabolic tasks
        complexity : str
            Implementation complexity (low, medium, high)
        speed : str
            Execution speed (fast, medium, slow)
        memory : str
            Memory usage (low, medium, high)
        reference : str
            Citation or reference
        """
        self.name = name
        self.method_class = method_class
        self.properties_class = properties_class
        self.description = description
        self.category = category
        self.requires_continuous_scores = requires_continuous_scores
        self.requires_threshold_scores = requires_threshold_scores
        self.supports_tasks = supports_tasks
        self.complexity = complexity
        self.speed = speed
        self.memory = memory
        self.reference = reference

    def to_dict(self) -> Dict[str, Any]:
        """Convert metadata to dictionary"""
        return {
            'name': self.name,
            'description': self.description,
            'category': self.category,
            'requires_continuous_scores': self.requires_continuous_scores,
            'requires_threshold_scores': self.requires_threshold_scores,
            'supports_tasks': self.supports_tasks,
            'complexity': self.complexity,
            'speed': self.speed,
            'memory': self.memory,
            'reference': self.reference
        }


class MethodRegistry:
    """
    Registry for omics integration methods
    오믹스 통합 방법론 레지스트리

    This class manages all available reconstruction methods and provides
    utilities for method discovery, instantiation, and comparison.

    사용 가능한 모든 재구성 방법을 관리하고 메서드 검색, 인스턴스화 및 비교를 위한
    유틸리티를 제공합니다.
    """

    _methods: Dict[str, MethodMetadata] = {}
    _hooks: Dict[str, List[Callable]] = {
        'pre_register': [],
        'post_register': [],
        'pre_create': [],
        'post_create': []
    }

    @classmethod
    def register(
        cls,
        name: str,
        method_class: Type,
        properties_class: Type,
        **kwargs
    ) -> None:
        """
        Register a new reconstruction method

        Parameters
        ----------
        name : str
            Unique method name
        method_class : Type
            Method class
        properties_class : Type
            Properties class
        **kwargs
            Additional metadata (description, category, etc.)

        Examples
        --------
        >>> MethodRegistry.register(
        ...     'gimme',
        ...     GIMME,
        ...     GIMMEProperties,
        ...     description='Gene Inactivity Moderated by Metabolism and Expression',
        ...     requires_continuous_scores=True,
        ...     speed='fast'
        ... )
        """
        # Run pre-register hooks
        for hook in cls._hooks['pre_register']:
            hook(name, method_class, properties_class, kwargs)

        # Create metadata
        metadata = MethodMetadata(
            name=name,
            method_class=method_class,
            properties_class=properties_class,
            **kwargs
        )

        # Register method
        cls._methods[name.lower()] = metadata

        # Run post-register hooks
        for hook in cls._hooks['post_register']:
            hook(name, metadata)

    @classmethod
    def unregister(cls, name: str) -> None:
        """Unregister a method"""
        if name.lower() in cls._methods:
            del cls._methods[name.lower()]

    @classmethod
    def get_method(cls, name: str) -> Optional[MethodMetadata]:
        """Get method metadata by name"""
        return cls._methods.get(name.lower())

    @classmethod
    def list_methods(
        cls,
        category: Optional[str] = None,
        requires_continuous_scores: Optional[bool] = None,
        requires_threshold_scores: Optional[bool] = None
    ) -> List[str]:
        """
        List available methods with optional filtering

        Parameters
        ----------
        category : str, optional
            Filter by category
        requires_continuous_scores : bool, optional
            Filter by continuous score requirement
        requires_threshold_scores : bool, optional
            Filter by threshold score requirement

        Returns
        -------
        List[str]
            List of method names
        """
        methods = cls._methods.values()

        if category is not None:
            methods = [m for m in methods if m.category == category]

        if requires_continuous_scores is not None:
            methods = [m for m in methods if m.requires_continuous_scores == requires_continuous_scores]

        if requires_threshold_scores is not None:
            methods = [m for m in methods if m.requires_threshold_scores == requires_threshold_scores]

        return [m.name for m in methods]

    @classmethod
    def get_all_metadata(cls) -> Dict[str, Dict[str, Any]]:
        """Get metadata for all registered methods"""
        return {name: meta.to_dict() for name, meta in cls._methods.items()}

    @classmethod
    def create_method(
        cls,
        name: str,
        S,
        lb,
        ub,
        properties,
        **kwargs
    ):
        """
        Create a method instance

        Parameters
        ----------
        name : str
            Method name
        S : array
            Stoichiometric matrix
        lb : array
            Lower bounds
        ub : array
            Upper bounds
        properties
            Properties object for the method
        **kwargs
            Additional arguments

        Returns
        -------
        Method instance

        Examples
        --------
        >>> gimme = MethodRegistry.create_method(
        ...     'gimme',
        ...     S=S_matrix,
        ...     lb=lower_bounds,
        ...     ub=upper_bounds,
        ...     properties=gimme_props
        ... )
        """
        # Run pre-create hooks
        for hook in cls._hooks['pre_create']:
            hook(name, S, lb, ub, properties, kwargs)

        metadata = cls.get_method(name)
        if metadata is None:
            raise ValueError(f"Method '{name}' not found in registry. "
                           f"Available methods: {', '.join(cls.list_methods())}")

        # Create instance
        instance = metadata.method_class(S=S, lb=lb, ub=ub, properties=properties, **kwargs)

        # Run post-create hooks
        for hook in cls._hooks['post_create']:
            hook(name, instance)

        return instance

    @classmethod
    def create_properties(
        cls,
        name: str,
        **kwargs
    ):
        """
        Create a properties instance for a method

        Parameters
        ----------
        name : str
            Method name
        **kwargs
            Properties parameters

        Returns
        -------
        Properties instance
        """
        metadata = cls.get_method(name)
        if metadata is None:
            raise ValueError(f"Method '{name}' not found in registry")

        return metadata.properties_class(**kwargs)

    @classmethod
    def add_hook(cls, hook_type: str, func: Callable) -> None:
        """
        Add a hook function

        Parameters
        ----------
        hook_type : str
            Hook type (pre_register, post_register, pre_create, post_create)
        func : Callable
            Hook function
        """
        if hook_type not in cls._hooks:
            raise ValueError(f"Invalid hook type: {hook_type}")
        cls._hooks[hook_type].append(func)

    @classmethod
    def remove_hook(cls, hook_type: str, func: Callable) -> None:
        """Remove a hook function"""
        if hook_type in cls._hooks and func in cls._hooks[hook_type]:
            cls._hooks[hook_type].remove(func)

    @classmethod
    def print_registry(cls) -> None:
        """Print all registered methods with their metadata"""
        print("=" * 80)
        print("Registered Omics Integration Methods")
        print("=" * 80)

        for name, metadata in sorted(cls._methods.items()):
            print(f"\n{name.upper()}")
            print("-" * 40)
            print(f"  Description: {metadata.description}")
            print(f"  Category: {metadata.category}")
            print(f"  Continuous scores: {metadata.requires_continuous_scores}")
            print(f"  Threshold scores: {metadata.requires_threshold_scores}")
            print(f"  Supports tasks: {metadata.supports_tasks}")
            print(f"  Complexity: {metadata.complexity}")
            print(f"  Speed: {metadata.speed}")
            print(f"  Memory: {metadata.memory}")
            if metadata.reference:
                print(f"  Reference: {metadata.reference}")

        print("\n" + "=" * 80)


def register_method(
    name: str,
    description: str = "",
    category: str = "reconstruction",
    **metadata_kwargs
):
    """
    Decorator for registering a method

    Examples
    --------
    >>> @register_method(
    ...     'my_method',
    ...     description='My custom method',
    ...     requires_continuous_scores=True
    ... )
    ... class MyMethod(ContextSpecificModelReconstructionAlgorithm):
    ...     pass
    """
    def decorator(cls):
        # Try to infer properties class from class attributes
        properties_class = getattr(cls, 'properties_class', None)
        if properties_class is None:
            # Look for a nested Properties class
            properties_class = getattr(cls, 'Properties', None)

        if properties_class is None:
            raise ValueError(f"Could not find properties class for {cls.__name__}")

        MethodRegistry.register(
            name=name,
            method_class=cls,
            properties_class=properties_class,
            description=description,
            category=category,
            **metadata_kwargs
        )

        return cls

    return decorator


# Auto-register existing methods
def auto_register_methods():
    """
    Automatically register existing Troppo methods
    기존 Troppo 방법론을 자동으로 등록합니다
    """
    try:
        from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
        MethodRegistry.register(
            'gimme',
            GIMME,
            GIMMEProperties,
            description='Gene Inactivity Moderated by Metabolism and Expression',
            requires_continuous_scores=True,
            speed='fast',
            memory='low',
            complexity='low',
            reference='Becker SA, Palsson BO. (2008). Genome Biology'
        )
    except ImportError:
        pass

    try:
        from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
        MethodRegistry.register(
            'tinit',
            tINIT,
            tINITProperties,
            description='Task-driven Integrative Network Inference for Tissues',
            requires_threshold_scores=True,
            supports_tasks=True,
            speed='medium',
            memory='medium',
            complexity='medium',
            reference='Agren et al. (2014). Molecular Systems Biology'
        )
    except ImportError:
        pass

    try:
        from troppo.methods.reconstruction.imat import IMAT, IMATProperties
        MethodRegistry.register(
            'imat',
            IMAT,
            IMATProperties,
            description='Integrative Metabolic Analysis Tool',
            requires_threshold_scores=True,
            speed='medium',
            memory='medium',
            complexity='medium',
            reference='Shlomi et al. (2008). Nature Biotechnology'
        )
    except ImportError:
        pass

    try:
        from troppo.methods.reconstruction.fastcore import FASTCORE, FASTCOREProperties
        MethodRegistry.register(
            'fastcore',
            FASTCORE,
            FASTCOREProperties,
            description='Fast Consistency-based Reconstruction',
            requires_threshold_scores=True,
            speed='fast',
            memory='low',
            complexity='low',
            reference='Vlassis et al. (2014). PLoS Computational Biology'
        )
    except ImportError:
        pass

    try:
        from troppo.methods.reconstruction.corda import CORDA, CORDAProperties
        MethodRegistry.register(
            'corda',
            CORDA,
            CORDAProperties,
            description='Cost Optimization Reaction Dependency Assessment',
            requires_threshold_scores=True,
            speed='slow',
            memory='high',
            complexity='high',
            reference='Schultz & Qutub (2016). PLoS Computational Biology'
        )
    except ImportError:
        pass

    try:
        from troppo.methods.reconstruction.swiftcore import SWIFTCORE, SWIFTCOREProperties
        MethodRegistry.register(
            'swiftcore',
            SWIFTCORE,
            SWIFTCOREProperties,
            description='SWIFTCORE - Fast CORE reconstruction',
            requires_threshold_scores=True,
            speed='fast',
            memory='low',
            complexity='low'
        )
    except ImportError:
        pass


# Auto-register on import
auto_register_methods()
