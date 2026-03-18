"""General exception and warning classes."""


class StructuralCodesWarning(Warning):
    """Base class for StructurlCodes warnings."""


class NoConvergenceWarning(StructuralCodesWarning):
    """Warning for lack of convergence.

    A warning that indicates that no convergence was reached for an
    iterative solver. The solution should be inspected with care.
    """


class InformationWarning(StructuralCodesWarning):
    """Warning for simple information.

    A warning that indicates that some information that should be taken care
    by user, but the calculation is not invalid.
    """
