"""General exception and warning classes."""


class StructuralCodesWarning(Warning):
    """Base class for StructurlCodes warnings."""


class NoConvergenceWarning(StructuralCodesWarning):
    """A warning that indicates that no convergence was reached for an
    iterative solver.
    """
