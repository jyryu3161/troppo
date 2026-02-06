import abc

from . import OmicsDataMap

MINSUM = (min, sum)
MINMAX = (min, max)


class ScoreIntegrationStrategy:
    """Base class for score integration strategies."""
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def integrate(self, data_map: OmicsDataMap):
        pass


class ReactionProtectionMixin:
    """Mixin for protecting specific reactions from removal."""

    def __init__(self, protected_reactions: list):
        self.protected_reactions = protected_reactions


class ContinuousScoreIntegrationStrategy(ScoreIntegrationStrategy):
    """Integrate continuous reaction scores directly."""

    def __init__(self, score_apply=None):
        self.score_apply = score_apply

    def integrate(self, data_map: OmicsDataMap) -> dict:
        return data_map.get_scores() if self.score_apply is None else self.score_apply(data_map.get_scores())


class CustomSelectionIntegrationStrategy(ScoreIntegrationStrategy):
    """Integrate scores using custom group functions."""

    def __init__(self, group_functions: dict):
        self.group_functions = group_functions

    def integrate(self, data_map: OmicsDataMap) -> dict:
        tvals = [f(data_map) for f in self.group_functions]
        return tvals[0] if len(tvals) < 2 else tvals


class AdjustedScoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    """Integrate scores with asymmetric normalization and reaction protection.

    Normalizes negative scores by dividing by the absolute value of the minimum
    negative score (not the maximum positive), producing symmetric [-1, max_pos] range.
    Protected reactions receive the maximum score.
    """

    def __init__(self, protected_reactions: list):
        ReactionProtectionMixin.__init__(self, protected_reactions)

    def integrate(self, data_map: OmicsDataMap) -> dict:
        scores_raw = data_map.get_scores()
        non_none_values = [v for v in scores_raw.values() if v is not None]

        if not non_none_values:
            return {k: 0 for k in scores_raw}

        # Normalize negative scores by abs(min_negative) to get symmetric range
        min_val = min(non_none_values)
        if min_val < 0:
            abs_min = abs(min_val)
            scores = {
                k: (v / abs_min if v < 0 else v) if v is not None else 0
                for k, v in scores_raw.items()
            }
        else:
            scores = {k: v if v is not None else 0 for k, v in scores_raw.items()}

        # Protected reactions get max score
        if self.protected_reactions:
            max_score = max(scores.values()) if scores else 1.0
            scores.update({x: max_score for x in self.protected_reactions if x in scores})

        return scores


class DefaultCoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    """Select reactions above a threshold as core reactions."""

    def __init__(self, threshold, protected_reactions: list):
        ReactionProtectionMixin.__init__(self, protected_reactions)
        self.__threshold = threshold

    def integrate(self, data_map: OmicsDataMap) -> list:
        return [[k for k, v in data_map.get_scores().items() if
                 (v is not None and v > self.__threshold) or k in self.protected_reactions]]


class ThresholdSelectionIntegrationStrategy(ScoreIntegrationStrategy):
    """Select reactions above one or more thresholds."""

    def __init__(self, thresholds):
        if isinstance(thresholds, (int, float)):
            self.thresholds = [thresholds]
        else:
            self.thresholds = thresholds

    def integrate(self, data_map: OmicsDataMap) -> list:
        tvals = [data_map.select(op='above', threshold=float(t)) for t in self.thresholds]
        return tvals[0] if len(tvals) < 2 else tvals
