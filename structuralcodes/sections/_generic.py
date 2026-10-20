"""The deprecated GenericSection class. Please use the BeamSection class
instead.
"""

import warnings

from ._beam_section import BeamSection


class GenericSection(BeamSection):
    """This is the GenericSection class which was renamed to BeamSection in
    v0.7.0.
    """

    def __init__(self, *args, **kwargs) -> None:
        """Since the GenericSection class was renamed to BeamSection in v0.7.0,
        this will initialize a BeamSection.
        """
        warnings.warn(
            (
                'The GenericSection class was renamed to BeamSection in '
                'v0.7.0. The GenericSection is scheduled for removal in a '
                'future release. Please use the BeamSection class instead. An '
                'object of type BeamSection is initialized.'
            ),
            category=DeprecationWarning,
        )
        super().__init__(*args, **kwargs)
