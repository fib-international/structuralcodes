import typing as t

import numpy as np
from scipy.stats import lognorm

corrosion_type = t.Literal['carbonation', 'chloride']
exposure_class = t.Literal[
    'sheltered',
    'unsheltered',
    'wet',
    'cyclic_dry_wet',
    'airborn_seawater',
    'submerged',
    'tidal_zone',
]
exposure_class_corrosion_type = {
    'carbonation': ['sheltered', 'unsheltered'],
    'chloride': [
        'wet',
        'cyclic_dry_wet',
        'airborn_seawater',
        'submerged',
        'tidal_zone',
    ],
}


def _validate_corrosion_exposure(
    corrosion_type: corrosion_type, exposure_class: exposure_class
) -> bool:
    """Returns True if corrosion_type and exposure_class are valid."""
    corr_type = corrosion_type.lower()
    exposure = exposure_class.lower()
    if corr_type not in ['carbonation', 'chloride']:
        return False
    return exposure in exposure_class_corrosion_type[corr_type]


_table_30_1_6 = {
    'carbonation': {
        'sheltered': {'mean': 2, 'std': 3},
        # the table shows "t", but 7 is correct
        'unsheltered': {'mean': 5, 'std': 7},
    },
    'chloride': {
        'wet': {'mean': 4, 'std': 6},
        'cyclic_dry_wet': {'mean': 30, 'std': 40},
        'airborn_seawater': {'mean': 30, 'std': 40},
        'submerged': {'mean': 4, 'std': 7},
        'tidal_zone': {'mean': 50, 'std': 100},
    },
}


def pcorr_rep(
    corrosion_type: corrosion_type, exposure_class: exposure_class
) -> dict:
    """Returns corrosion rate properties in micrometer/year.

    This implements MC 2020 table 30.1-6a and 30.1-6b.
    The return is a dictionary containing mean and standard deviation of pcorr
    in micrometer/year.

    Note:
        There is a typo in MC2020 Table 30.1-6a. In unshelterd conditions, the
        pcorr standard deviation is indicated as "t" while it should be equal
        to 7 micrometer/year according to fib Bulletin 111.

    Args:
        corrosion_type (str): The type of corrosion. Can be either
            "carbonation" or "chloride".
        exposure_class (str): The exposure class. Can be one of the following
            for "carbonation": "sheltered", "unsheltered".
            Can be one of the following for "chloride":
            "wet", "cyclic_dry_wet", "airborn_seawater", "submerged",
            "tidal_zone".

    Returns:
        dict: A dictionary containing the mean and standard deviation of the
            corrosion rate in micrometer/year.

    Raises:
        ValueError: if the corrosion type is not valid.
        ValueError: if the exposure class is not valid for the relative
            corrosion type.
    """
    if not _validate_corrosion_exposure(corrosion_type, exposure_class):
        raise ValueError(
            'Invalid corrosion_type or exposure_class.\n'
            'corrosion_type must be either "carbonation" or "chloride".\n'
            'exposure_class must be one of the following for "carbonation": '
            '"sheltered", "unsheltered".\n'
            'exposure_class must be one of the following for "chloride": '
            '"wet", "cyclic_dry_wet", "airborn_seawater", "submerged", '
            '"tidal_zone".'
        )
    return _table_30_1_6[corrosion_type.lower()][exposure_class.lower()]


def pcorr_fractile(
    corrosion_type: corrosion_type,
    exposure_class: exposure_class,
    fractile: float = 0.5,
) -> float:
    """Returns the fractile of the corrosion rate in micrometer/year.

    This functions uses data from MC 2020 table 30.1-6a and 30.1-6b to
    calculate the corrosion rate at a given fractile. The return is a float
    representing the corrosion rate in micrometer/year.

    Note: according to fib Bulletin 111, the corrosion rate can be assumed to
        follow a lognormal distribution.

    Args:
        corrosion_type (str): The type of corrosion. Can be either
            "carbonation" or "chloride".
        exposure_class (str): The exposure class. Can be one of the following
            for "carbonation": "sheltered", "unsheltered".
            Can be one of the following for "chloride":
            "wet", "cyclic_dry_wet", "airborn_seawater", "submerged",
            "tidal_zone".
        fractile (Optional[float]): The fractile for which to calculate the
            corrosion rate. The default value is 0.5, which corresponds to the
            median corrosion rate. The value should be in the range ]0,1[

    Returns:
        float: The corrosion rate at the specified fractile in micrometer/year.

    Raises:
        ValueError: if the corrosion type is not valid.
        ValueError: if the exposure class is not valid for the relative
            corrosion type.
        ValueError: if the fractile value is not valid.
    """
    if fractile <= 0 or fractile >= 1:
        raise ValueError(
            'The value of fractile is not valid, use a value between 0 and 1.'
        )
    pcorr_dic = pcorr_rep(corrosion_type, exposure_class)

    # Lognormal distribution
    s = np.log(1 + pcorr_dic['std'] ** 2 / pcorr_dic['mean'] ** 2) ** 0.5
    scale = (
        pcorr_dic['mean'] ** 2
        / (pcorr_dic['mean'] ** 2 + pcorr_dic['std'] ** 2) ** 0.5
    )
    lognorm_dist = lognorm(s=s, scale=scale)
    return lognorm_dist.ppf(fractile)


def kc(fc: float) -> float:
    """Compute kc reducing compressive strength factor.

    MC2020 30.1.10.3.2.1

    It does depend only on fc and not on corrosion level, even if it is valid
    for low to moderate corrosion values.

    Moderate and low corrosion levels have been defined to apply approximately
    for medium bar diameters with 5% weight loss and 0.25 depth of corrosion.

    Args:
        fc (float): Compressive strength of concrete in MPa.

    Returns:
        float: The reducing compressive strength factor kc.

    Raises:
        ValueError: if fc is not positive.
    """
    if fc <= 0:
        raise ValueError('fc must be a positive value.')
    eta_fc = min((30 / fc) ** (1 / 3.0), 1)
    return 0.75 * eta_fc
