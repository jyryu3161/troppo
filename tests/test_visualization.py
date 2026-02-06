import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing

from troppo.results import ReconstructionResults


@pytest.fixture
def mock_results():
    """Create mock results for visualization testing."""
    model = MagicMock()
    model.reactions = [
        MagicMock(id=f'R{i}', subsystem=f'Subsystem_{i % 5}') for i in range(20)
    ]
    model.metabolites = [MagicMock(id=f'M{i}') for i in range(10)]
    model.genes = [MagicMock(id=f'G{i}') for i in range(5)]

    np.random.seed(42)
    results_dict = {}
    for i in range(5):
        sample_name = f'sample_{i}'
        results_dict[sample_name] = {
            f'R{j}': bool(np.random.random() > 0.3) for j in range(20)
        }
    return ReconstructionResults(results_dict, model)


class TestVisualization:
    def test_plot_model_sizes(self, mock_results):
        from troppo.visualization import plot_model_sizes
        fig = plot_model_sizes(mock_results)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_plot_reaction_overlap_heatmap(self, mock_results):
        from troppo.visualization import plot_reaction_overlap_heatmap
        fig = plot_reaction_overlap_heatmap(mock_results)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_plot_subsystem_coverage(self, mock_results):
        from troppo.visualization import plot_subsystem_coverage
        fig = plot_subsystem_coverage(mock_results, top_n=5)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_plot_reaction_frequency(self, mock_results):
        from troppo.visualization import plot_reaction_frequency
        fig = plot_reaction_frequency(mock_results, top_n=10)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_plot_overview(self, mock_results):
        from troppo.visualization import plot_overview
        fig = plot_overview(mock_results)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_plot_overview_save(self, mock_results, tmp_path):
        from troppo.visualization import plot_overview
        save_dir = str(tmp_path / "figures")
        fig = plot_overview(mock_results, save_dir=save_dir)
        assert fig is not None
        import os
        assert os.path.exists(os.path.join(save_dir, 'overview.pdf'))
        import matplotlib.pyplot as plt
        plt.close(fig)
