import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch

from troppo.results import ReconstructionResults


class _PicklableReaction:
    """Minimal picklable reaction for testing save/load."""
    def __init__(self, id, subsystem=''):
        self.id = id
        self.subsystem = subsystem


class _PicklableModel:
    """Minimal picklable model for testing save/load."""
    def __init__(self):
        self.reactions = [
            _PicklableReaction('R1', 'Glycolysis'),
            _PicklableReaction('R2', 'TCA cycle'),
        ]
        self.metabolites = []
        self.genes = []


class TestReconstructionResults:
    """Test ReconstructionResults with mock data."""

    @pytest.fixture
    def mock_model(self):
        """Create a mock COBRA model."""
        model = MagicMock()
        model.reactions = [
            MagicMock(id='R1', subsystem='Glycolysis'),
            MagicMock(id='R2', subsystem='Glycolysis'),
            MagicMock(id='R3', subsystem='TCA cycle'),
            MagicMock(id='R4', subsystem='TCA cycle'),
            MagicMock(id='R5', subsystem='Pentose phosphate'),
        ]
        model.metabolites = [MagicMock(id=f'M{i}') for i in range(10)]
        model.genes = [MagicMock(id=f'G{i}') for i in range(5)]
        return model

    @pytest.fixture
    def sample_results(self, mock_model):
        """Create sample reconstruction results."""
        results_dict = {
            'sample1': {'R1': True, 'R2': True, 'R3': True, 'R4': False, 'R5': False},
            'sample2': {'R1': True, 'R2': False, 'R3': True, 'R4': True, 'R5': False},
            'sample3': {'R1': True, 'R2': True, 'R3': True, 'R4': True, 'R5': True},
        }
        return ReconstructionResults(results_dict, mock_model, metadata={'solver': 'CPLEX'})

    def test_n_samples(self, sample_results):
        assert sample_results.n_samples == 3

    def test_sample_names(self, sample_results):
        assert set(sample_results.sample_names) == {'sample1', 'sample2', 'sample3'}

    def test_reaction_matrix_shape(self, sample_results):
        assert sample_results.reaction_matrix.shape == (3, 5)

    def test_n_reactions_per_sample(self, sample_results):
        counts = sample_results.n_reactions_per_sample
        assert counts['sample1'] == 3
        assert counts['sample2'] == 3
        assert counts['sample3'] == 5

    def test_get_active_reactions(self, sample_results):
        active = sample_results.get_active_reactions('sample1')
        assert set(active) == {'R1', 'R2', 'R3'}

    def test_get_common_reactions(self, sample_results):
        common = sample_results.get_common_reactions()
        assert common == {'R1', 'R3'}

    def test_get_unique_reactions(self, sample_results):
        unique = sample_results.get_unique_reactions('sample3')
        assert unique == {'R5'}

    def test_summary(self, sample_results):
        summary = sample_results.summary()
        assert 'n_active_reactions' in summary.columns
        assert 'pct_retained' in summary.columns
        assert len(summary) == 3

    def test_subsystem_coverage(self, sample_results):
        coverage = sample_results.subsystem_coverage()
        assert len(coverage) == 3  # 3 samples
        assert 'Glycolysis' in coverage.columns

    def test_jaccard_similarity(self, sample_results):
        jac = sample_results.jaccard_similarity()
        assert jac.shape == (3, 3)
        # Diagonal should be 1.0
        for name in sample_results.sample_names:
            assert jac.loc[name, name] == pytest.approx(1.0)

    def test_to_csv(self, sample_results, tmp_path):
        csv_path = str(tmp_path / "results.csv")
        sample_results.to_csv(csv_path)
        df = pd.read_csv(csv_path, index_col=0)
        assert df.shape == (3, 5)

    def test_save_load(self, tmp_path):
        """Test save/load with a pickle-compatible model (no MagicMock)."""
        model = _PicklableModel()
        results_dict = {
            'sample1': {'R1': True, 'R2': False},
            'sample2': {'R1': True, 'R2': True},
        }
        results = ReconstructionResults(results_dict, model, metadata={'solver': 'CPLEX'})

        pkl_path = str(tmp_path / "results.pkl")
        results.save(pkl_path)
        loaded = ReconstructionResults.load(pkl_path)
        assert loaded.n_samples == 2
        assert set(loaded.sample_names) == {'sample1', 'sample2'}

    def test_repr(self, sample_results):
        r = repr(sample_results)
        assert 'n_samples=3' in r
        assert 'n_reactions=5' in r
