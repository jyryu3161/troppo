import numpy as np
import pandas as pd
import pickle
from typing import Optional, List, Set, Dict
import warnings


class ReconstructionResults:
    """Container for multi-sample tINIT reconstruction results.

    Stores the boolean reaction activity matrix across samples and provides
    methods for querying, comparing, validating, and exporting reconstructed
    context-specific metabolic models.

    Parameters
    ----------
    results_dict : Dict[str, Dict[str, bool]]
        Mapping of ``{sample_name: {reaction_id: True/False}}`` for each
        sample processed by the tINIT pipeline.
    model : cobra.Model
        The base COBRA model used for reconstruction.  A reference is kept
        so that sample-specific sub-models can be derived on demand.
    metadata : dict, optional
        Pipeline metadata (solver name, algorithm parameters, timestamps,
        etc.) carried along for provenance tracking.
    """

    def __init__(
        self,
        results_dict: Dict[str, Dict[str, bool]],
        model,
        metadata: dict = None,
    ):
        # Store model reference
        self._model = model
        self._metadata = metadata or {}

        # Build reaction matrix (samples x reactions boolean DataFrame)
        self._reaction_matrix = pd.DataFrame(results_dict).T
        self._reaction_matrix = self._reaction_matrix.fillna(False).astype(bool)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def reaction_matrix(self) -> pd.DataFrame:
        """Boolean DataFrame with shape ``(n_samples, n_reactions)``.

        Rows are sample names, columns are reaction identifiers.  A ``True``
        entry means the reaction was retained in the corresponding
        context-specific reconstruction.
        """
        return self._reaction_matrix

    @property
    def sample_names(self) -> list:
        """Ordered list of sample names."""
        return list(self._reaction_matrix.index)

    @property
    def n_samples(self) -> int:
        """Number of samples in the result set."""
        return len(self._reaction_matrix)

    @property
    def n_reactions_per_sample(self) -> pd.Series:
        """Number of active reactions per sample."""
        return self._reaction_matrix.sum(axis=1)

    # ------------------------------------------------------------------
    # Model retrieval
    # ------------------------------------------------------------------

    def get_model(self, sample: str):
        """Return a COBRA model for a specific sample.

        Inactive reactions are removed from a copy of the base model and
        orphan metabolites are cleaned up.

        Parameters
        ----------
        sample : str
            Sample name (must exist in :pyattr:`sample_names`).

        Returns
        -------
        cobra.Model
            A pruned copy of the base model containing only the reactions
            that were active in *sample*.
        """
        import cobra

        model_copy = self._model.copy()
        active = set(
            self._reaction_matrix.loc[sample][
                self._reaction_matrix.loc[sample]
            ].index
        )
        to_remove = [r for r in model_copy.reactions if r.id not in active]
        model_copy.remove_reactions(to_remove, remove_orphans=True)
        return model_copy

    # ------------------------------------------------------------------
    # Reaction queries
    # ------------------------------------------------------------------

    def get_active_reactions(self, sample: str) -> list:
        """List of reaction IDs active in *sample*.

        Parameters
        ----------
        sample : str
            Sample name.

        Returns
        -------
        list of str
            Reaction identifiers retained during reconstruction.
        """
        row = self._reaction_matrix.loc[sample]
        return list(row[row].index)

    def get_common_reactions(self) -> set:
        """Reactions active in **all** samples.

        Returns
        -------
        set of str
            Reaction identifiers that are universally retained.
        """
        return set(
            self._reaction_matrix.columns[self._reaction_matrix.all(axis=0)]
        )

    def get_unique_reactions(self, sample: str) -> set:
        """Reactions active **only** in the given sample.

        Parameters
        ----------
        sample : str
            Sample name.

        Returns
        -------
        set of str
            Reaction identifiers that are active in *sample* and inactive
            in every other sample.
        """
        others = self._reaction_matrix.drop(sample)
        active_here = set(
            self._reaction_matrix.loc[sample][
                self._reaction_matrix.loc[sample]
            ].index
        )
        active_elsewhere = set(
            self._reaction_matrix.columns[others.any(axis=0)]
        )
        return active_here - active_elsewhere

    # ------------------------------------------------------------------
    # Summary & comparison
    # ------------------------------------------------------------------

    def summary(self) -> pd.DataFrame:
        """Per-sample summary statistics.

        Returns
        -------
        pd.DataFrame
            Columns: ``n_active_reactions``, ``n_total_reactions``,
            ``pct_retained``.  Indexed by sample name.
        """
        data = []
        for sample in self.sample_names:
            n_active = self._reaction_matrix.loc[sample].sum()
            n_total = len(self._reaction_matrix.columns)
            data.append(
                {
                    "sample": sample,
                    "n_active_reactions": int(n_active),
                    "n_total_reactions": n_total,
                    "pct_retained": round(100.0 * n_active / n_total, 1),
                }
            )
        return pd.DataFrame(data).set_index("sample")

    def subsystem_coverage(self) -> pd.DataFrame:
        """Subsystem-level active reaction fraction per sample.

        For every metabolic subsystem defined in the base model, computes the
        fraction of reactions that are active in each sample.

        Returns
        -------
        pd.DataFrame
            Shape ``(n_samples, n_subsystems)`` with values in ``[0, 1]``.
        """
        rxn_subsystems = {}
        for rxn in self._model.reactions:
            sub = getattr(rxn, "subsystem", "") or "Unknown"
            rxn_subsystems[rxn.id] = sub

        # Build subsystem x sample matrix
        subsystems = sorted(set(rxn_subsystems.values()))
        data = {}
        for sample in self.sample_names:
            active = set(self.get_active_reactions(sample))
            sample_data = {}
            for sub in subsystems:
                sub_rxns = [r for r, s in rxn_subsystems.items() if s == sub]
                if sub_rxns:
                    sample_data[sub] = sum(
                        1 for r in sub_rxns if r in active
                    ) / len(sub_rxns)
                else:
                    sample_data[sub] = 0.0
            data[sample] = sample_data
        return pd.DataFrame(data).T

    def jaccard_similarity(self) -> pd.DataFrame:
        """Pairwise Jaccard similarity matrix between samples.

        Returns
        -------
        pd.DataFrame
            Symmetric ``(n_samples, n_samples)`` matrix with Jaccard indices
            in ``[0, 1]``.
        """
        n = self.n_samples
        names = self.sample_names
        matrix = np.zeros((n, n))
        for i in range(n):
            set_i = set(self.get_active_reactions(names[i]))
            for j in range(i, n):
                set_j = set(self.get_active_reactions(names[j]))
                if set_i or set_j:
                    jac = len(set_i & set_j) / len(set_i | set_j)
                else:
                    jac = 1.0
                matrix[i, j] = matrix[j, i] = jac
        return pd.DataFrame(matrix, index=names, columns=names)

    # ------------------------------------------------------------------
    # Flux analysis
    # ------------------------------------------------------------------

    def biomass_fluxes(self) -> pd.Series:
        """Compute the objective (biomass) flux for each sample model.

        Each sample-specific model is optimised with the default solver
        settings.  If optimisation fails the flux is recorded as ``0.0``.

        Returns
        -------
        pd.Series
            Indexed by sample name with the objective value.
        """
        fluxes = {}
        for sample in self.sample_names:
            try:
                m = self.get_model(sample)
                sol = m.optimize()
                fluxes[sample] = (
                    sol.objective_value if sol.status == "optimal" else 0.0
                )
            except Exception:
                fluxes[sample] = 0.0
        return pd.Series(fluxes, name="biomass_flux")

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------

    def validate_gene_essentiality(
        self,
        essential_genes: list,
        non_essential_genes: list = None,
    ) -> pd.DataFrame:
        """Validate gene essentiality predictions across samples.

        Uses :class:`troppo.benchmark.validation.GeneEssentialityValidator`
        to perform single-gene knockouts and compare predicted essentiality
        against known essential / non-essential gene lists.

        Parameters
        ----------
        essential_genes : list of str
            Gene identifiers known to be essential.
        non_essential_genes : list of str, optional
            Gene identifiers known to be non-essential.

        Returns
        -------
        pd.DataFrame
            Validation metrics per sample (accuracy, sensitivity,
            specificity, etc.).
        """
        from troppo.benchmark.validation import GeneEssentialityValidator

        validator = GeneEssentialityValidator(
            essential_genes=essential_genes,
            non_essential_genes=non_essential_genes or [],
        )
        results = []
        for sample in self.sample_names:
            m = self.get_model(sample)
            result = validator.validate(m)
            row = result.to_dict()
            row["sample"] = sample
            results.append(row)
        return pd.DataFrame(results).set_index("sample")

    def validate_theoretical_yields(
        self,
        substrate: str,
        product: str,
        conditions: dict = None,
    ) -> pd.DataFrame:
        """Validate theoretical product yields across samples.

        Uses :class:`troppo.benchmark.validation.TheoreticalYieldCalculator`
        to compute maximum theoretical yields for a given substrate/product
        pair in each reconstructed model.

        Parameters
        ----------
        substrate : str
            Exchange reaction or metabolite identifier for the carbon source.
        product : str
            Exchange reaction or metabolite identifier for the target product.
        conditions : dict, optional
            Additional medium or constraint overrides.

        Returns
        -------
        pd.DataFrame
            Yield results per sample.
        """
        from troppo.benchmark.validation import TheoreticalYieldCalculator

        calc = TheoreticalYieldCalculator()
        results = []
        for sample in self.sample_names:
            m = self.get_model(sample)
            yield_results = calc.calculate_yields(
                m, carbon_sources=[substrate], product=product
            )
            for r in yield_results:
                row = r.to_dict()
                row["sample"] = sample
                results.append(row)
        return pd.DataFrame(results).set_index("sample")

    # ------------------------------------------------------------------
    # Export / serialisation
    # ------------------------------------------------------------------

    def to_csv(self, path: str):
        """Export the reaction matrix to a CSV file.

        Parameters
        ----------
        path : str
            Destination file path.
        """
        self._reaction_matrix.to_csv(path)

    def to_sbml(self, output_dir: str):
        """Export each sample-specific model as an SBML file.

        Parameters
        ----------
        output_dir : str
            Directory in which to write ``<sample>.xml`` files.  Created if
            it does not already exist.
        """
        import os
        import cobra

        os.makedirs(output_dir, exist_ok=True)
        for sample in self.sample_names:
            m = self.get_model(sample)
            cobra.io.write_sbml_model(
                m, os.path.join(output_dir, f"{sample}.xml")
            )

    def save(self, path: str):
        """Persist the results object via :mod:`pickle`.

        Parameters
        ----------
        path : str
            Destination file path.
        """
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: str) -> "ReconstructionResults":
        """Load a previously saved :class:`ReconstructionResults`.

        Parameters
        ----------
        path : str
            Path to the pickled file.

        Returns
        -------
        ReconstructionResults
        """
        with open(path, "rb") as f:
            return pickle.load(f)

    # ------------------------------------------------------------------
    # Dunder methods
    # ------------------------------------------------------------------

    def __repr__(self):
        return (
            f"ReconstructionResults(n_samples={self.n_samples}, "
            f"n_reactions={len(self._reaction_matrix.columns)})"
        )
