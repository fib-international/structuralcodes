"""Functions from Section 6 of EN 1992-1-1:2023."""

import math
from typing import Dict, List, Literal, Optional

from scipy.interpolate import interp1d


def get_exposure_classes(
    env_level: Optional[int] = None,
    env_name: Optional[
        Literal[
            'none',
            'carbonation',
            'chlorides',
            'sea',
            'freeze',
            'chemical',
            'abrasion',
        ]
    ] = None,
) -> List[str]:
    """Returns a list with valid exposure classes.

    EN1992-1-1:2023 Table (6.1).

    Args:
        env_level (int, optional): Filters the result by only a determined env
            level.
        env_name (str, optional): The name of the environment agent.

    Returns:
        (List[str]): List with valid exposure classes.
    """
    if env_level is not None and env_level not in range(1, 8):
        raise ValueError(
            f'env_level should be between 1 and 7. Got {env_level}'
        )

    exposure_classes = {
        1: ['X0'],
        2: ['XC1', 'XC2', 'XC3', 'XC4'],
        3: ['XD1', 'XD2', 'XD3'],
        4: ['XS1', 'XS2', 'XS3'],
        5: ['XF1', 'XF2', 'XF3', 'XF4'],
        6: ['XA1', 'XA2', 'XA3'],
        7: ['XM1', 'XM2', 'XM3'],
    }

    env_name_map = {
        'none': 1,
        'carbonation': 2,
        'chlorides': 3,
        'sea': 4,
        'freeze': 5,
        'chemical': 6,
        'abrasion': 7,
    }

    # Determine the corresponding environment level if env_name is provided
    if env_name:
        env_level = env_name_map.get(env_name)

    if env_level:
        return exposure_classes.get(env_level, [])

    # Return all exposure classes if no specific level or name is provided
    return sum(exposure_classes.values(), [])


def get_exposure_class_description(exp_class: str) -> Dict[str, str]:
    """Returns the exposure class description.

    EN1992-1-1:2023 Table (6.1).

    Args:
        exp_class (str): Exposure class code.

    Returns:
        (Dict[str,str]): A dictinary containing the environment, the
        description and the examples.

    Raise:
        ValueError: If exposure class does not exist.
    """
    data = {
        'X0': {
            'main_description': 'No risk of corrosion or attack.',
            'environment': 'All exposures except where there is freeze/thaw, abrasion or chemical attack.',  # noqa: E501
            'examples': 'Plain concrete members without any reinforcement.',
        },
        'XC1': {
            'main_description': 'Corrosion of embedded metal induced by carbonation.',  # noqa: E501
            'environment': 'Dry.',
            'examples': 'Concrete inside buildings with low air humidity, where the corrosion rate will be insignificant.',  # noqa: E501
        },
        'XC2': {
            'main_description': 'Corrosion of embedded metal induced by carbonation.',  # noqa: E501
            'environment': 'Wet or permanent high humidity, rarely dry.',
            'examples': 'Concrete surfaces subject to long-term water contact or permanently submerged in water or permanently exposed to high humidity; many foundations; water containments (not external).',  # noqa: E501
        },
        'XC3': {
            'main_description': 'Corrosion of embedded metal induced by carbonation.',  # noqa: E501
            'environment': 'Moderate humidity.',
            'examples': 'Concrete inside buildings with moderate humidity and not permanent high humidity; External concrete sheltered from rain.',  # noqa: E501
        },
        'XC4': {
            'main_description': 'Corrosion of embedded metal induced by carbonation.',  # noqa: E501
            'environment': 'Cyclic wet and dry.',
            'examples': 'Concrete surfaces subject to cyclic water contact (e.g. external concrete not sheltered from rain as walls and facades).',  # noqa: E501
        },
        'XD1': {
            'main_description': 'Corrosion of embedded metal induced by chlorides, excluding sea water.',  # noqa: E501
            'environment': 'Moderate humidity.',
            'examples': 'Concrete surfaces exposed to airborne chlorides.',
        },
        'XD2': {
            'main_description': 'Corrosion of embedded metal induced by chlorides, excluding sea water.',  # noqa: E501
            'environment': 'Wet, rarely dry.',
            'examples': 'Swimming pools; Concrete components exposed to industrial waters containing chlorides.',  # noqa: E501
        },
        'XD3': {
            'main_description': 'Corrosion of embedded metal induced by chlorides, excluding sea water.',  # noqa: E501
            'environment': 'Cyclic wet and dry.',
            'examples': 'Parts of bridges exposed to water containing chlorides; Concrete roads, pavements and car park slabs in areas where de-icing agents are frequently used.',  # noqa: E501
        },
        'XS1': {
            'main_description': 'Corrosion of embedded metal induced by chlorides from sea water.',  # noqa: E501
            'environment': 'Exposed to airborne salt but not in direct contact with sea water.',  # noqa: E501
            'examples': 'Structures near to or on the coast.',
        },
        'XS2': {
            'main_description': 'Corrosion of embedded metal induced by chlorides from sea water.',  # noqa: E501
            'environment': 'Permanently submerged.',
            'examples': 'Parts of marine structures and structures in seawater.',  # noqa: E501
        },
        'XS3': {
            'main_description': 'Corrosion of embedded metal induced by chlorides from sea water.',  # noqa: E501
            'environment': 'Tidal, splash and spray zones.',
            'examples': 'Parts of marine structures and structures temporarily or permanently directly over sea water.',  # noqa: E501
        },
        'XF1': {
            'main_description': 'Freeze/Thaw Attack.',
            'environment': 'Moderate water saturation, without de-icing agent.',  # noqa: E501
            'examples': 'Vertical concrete surfaces exposed to rain and freezing.',  # noqa: E501
        },
        'XF2': {
            'main_description': 'Freeze/Thaw Attack.',
            'environment': 'Moderate water saturation, with de-icing agent.',
            'examples': 'Vertical concrete surfaces of road structures exposed to freezing and airborne de-icing agents.',  # noqa: E501
        },
        'XF3': {
            'main_description': 'Freeze/Thaw Attack.',
            'environment': 'High water saturation, without de-icing agents.',
            'examples': 'Horizontal concrete surfaces exposed to rain and freezing.',  # noqa: E501
        },
        'XF4': {
            'main_description': 'Freeze/Thaw Attack.',
            'environment': 'High water saturation with de-icing agents or sea water.',  # noqa: E501
            'examples': 'Road and bridge decks exposed to de-icing agents and freezing; concrete surfaces exposed to direct spray containing de-icing agents and freezing; splash zone of marine structures exposed to freezing.',  # noqa: E501
        },
        'XA1': {
            'main_description': 'Chemical attack.',
            'environment': 'Slightly aggressive chemical environment.',
            'examples': 'Natural soils and ground water according to Table 6.2.',  # noqa: E501
        },
        'XA2': {
            'main_description': 'Chemical attack.',
            'environment': 'Moderately aggressive chemical environment.',
            'examples': 'Natural soils and ground water according to Table 6.2.',  # noqa: E501
        },
        'XA3': {
            'main_description': 'Chemical attack.',
            'environment': 'Highly aggressive chemical environment.',
            'examples': 'Natural soils and ground water according to Table 6.2.',  # noqa: E501
        },
        'XM1': {
            'main_description': 'Mechanical attack of concrete by abrasion.',
            'environment': 'Moderate abrasion.',
            'examples': 'Members of industrial sites frequented by vehicles with pneumatic tyres.',  # noqa: E501
        },
        'XM2': {
            'main_description': 'Mechanical attack of concrete by abrasion.',
            'environment': 'Heavy abrasion.',
            'examples': 'Members of industrial sites frequented by forklifts with pneumatic or solid rubber tyres.',  # noqa: E501
        },
        'XM3': {
            'main_description': 'Mechanical attack of concrete by abrasion.',
            'environment': 'Extreme abrasion.',
            'examples': 'Members of industrial sites frequented by forklifts with elastomer or steel tyres or track vehicles.',  # noqa: E501
        },
    }

    exp_class_upper = exp_class.upper()
    if exp_class_upper not in data:
        raise ValueError(f"Exposure class '{exp_class}' does not exist.")

    return data[exp_class_upper]


def c_nom(c_min: float, delta_c_dev: float) -> float:
    """Calculate the nominal cover.

    EN1992-1-1:2023 Eq. (6.1).

    Args:
        c_min (float): Minimum cover in mm.
        delta_c_dev (float): Allowance in design for deviation in mm.

    Returns:
        float: Nominal cover in mm.

    Raises:
        ValueError: If any of the inputs are negative.
    """
    if c_min < 0:
        raise ValueError(f'cmin must not be negative. Got {c_min}')
    if delta_c_dev < 0:
        raise ValueError(f'delta_cdev must not be negative. Got {delta_c_dev}')

    return c_min + delta_c_dev


def c_min(
    c_min_dur: float,
    sum_c: float,
    c_min_b: float,
    additional_cover: float = 0,
) -> float:
    """Calculate the minimum cover.

    EN1992-1-1:2023 Eq. (6.2).

    Args:
        cmin_dur (float): Minimum cover required for environmental conditions
            in mm.
        sum_c (float): Sum of applicable reductions and additions in mm.
        c_min_b (float): Minimum cover for bond requirement in mm.
        additional_cover (float, optional): Additional cover for casting
            against soil or other special requirements in mm. Default is 0.

    Returns:
        float: Minimum cover in mm.

    Raises:
        ValueError: If any of the inputs are negative (except sum_c).
    """
    if c_min_dur < 0:
        raise ValueError(f'cmin_dur must not be negative. Got {c_min_dur}')
    if c_min_b < 0:
        raise ValueError(f'cmin_b must not be negative. Got {c_min_b}')
    if additional_cover < 0:
        raise ValueError(
            f'additional_cover must not be negative. Got {additional_cover}'
        )

    return max(c_min_dur + sum_c, c_min_b, 10) + additional_cover


def c_min_dur_carb(
    exposure_class: Literal['XC1', 'XC2', 'XC3', 'XC4'],
    design_service_life: int,
    exposure_resistance_class: Literal[
        'XRC 0.5',
        'XRC 1',
        'XRC 2',
        'XRC 3',
        'XRC 4',
        'XRC 5',
        'XRC 6',
        'XRC 7',
    ],
    delta_c_dur_gamma: float = 0,
) -> float:
    """Calculate the minimum concrete cover for durability against carbonation.

    EN1992-1-1:2023 Table. (6.3).

    Args:
        exposure_class (Literal): Exposure class for carbonation (XC1, XC2,
            XC3, XC4).
        design_service_life (Literal): Design service life in years (50 or
            100).
        exposure_resistance_class (Literal): Exposure resistance class (XRC
            0.5, XRC 1, XRC 2, XRC 3, XRC 4, XRC 5, XRC 6, XRC 7).
        delta_c_dur_gamma (float): Increment in the concrete cover in mm
            considering special requirements.

    Returns:
        float: Minimum concrete cover in mm.

    Raise:
        ValueError: If delta_c_dur_gamma is less than zero.
    """
    if delta_c_dur_gamma < 0:
        raise ValueError(
            f'delta_c_dur_gamma must be positive. Got {delta_c_dur_gamma}'
        )

    # Clamp the design service life value
    design_service_life = max(50, min(design_service_life, 100))

    # Minimum covers from Table 6.3
    covers = {
        'XRC 0.5': {
            'XC1': [10, 10],
            'XC2': [10, 10],
            'XC3': [10, 10],
            'XC4': [10, 10],
        },
        'XRC 1': {
            'XC1': [10, 10],
            'XC2': [10, 10],
            'XC3': [10, 15],
            'XC4': [10, 15],
        },
        'XRC 2': {
            'XC1': [10, 15],
            'XC2': [10, 15],
            'XC3': [15, 25],
            'XC4': [15, 25],
        },
        'XRC 3': {
            'XC1': [10, 15],
            'XC2': [15, 20],
            'XC3': [20, 30],
            'XC4': [20, 30],
        },
        'XRC 4': {
            'XC1': [10, 20],
            'XC2': [15, 25],
            'XC3': [25, 35],
            'XC4': [25, 40],
        },
        'XRC 5': {
            'XC1': [15, 25],
            'XC2': [20, 30],
            'XC3': [25, 45],
            'XC4': [30, 45],
        },
        'XRC 6': {
            'XC1': [15, 25],
            'XC2': [25, 35],
            'XC3': [35, 55],
            'XC4': [40, 55],
        },
        'XRC 7': {
            'XC1': [15, 30],
            'XC2': [25, 40],
            'XC3': [40, 60],
            'XC4': [45, 60],
        },
    }

    # Determine base cover
    base_cover_v = covers[exposure_resistance_class][exposure_class]
    base_cover = interp1d([50, 100], base_cover_v)(design_service_life)

    # Calculate minimum cover
    c_min_dur = base_cover + delta_c_dur_gamma
    return max(c_min_dur, 0)


def c_min_dur_chlo(
    exposure_class: Literal['XS1', 'XS2', 'XS3', 'XD1', 'XD2', 'XD3'],
    design_service_life: int,
    exposure_resistance_class: Literal[
        'XRDS 0.5',
        'XRDS 1',
        'XRDS 1.5',
        'XRDS 2',
        'XRDS 3',
        'XRDS 4',
        'XRDS 5',
        'XRDS 6',
        'XRDS 8',
        'XRDS 10',
    ],
    delta_c_dur_gamma: float = 0.0,
) -> float:
    """Calculate the minimum concrete cover for durability against chlorides.

    EN1992-1-1:2023 Table. (6.4).

    Args:
        exposure_class (Literal): Exposure class for chlorides (XS1, XS2, XS3,
            XD1, XD2, XD3).
        design_service_life (Literal): Design service life in years (50 or
            100).
        exposure_resistance_class (Literal): Exposure resistance class (XRDS
            0.5, XRDS 1, XRDS 1.5, XRDS 2, XRDS 3, XRDS 4, XRDS 5, XRDS 6, XRDS
            8, XRDS 10).
        delta_c_dur_gamma (float): Increment in the concrete cover in mm
            considering special requirements.

    Returns:
        float: Minimum concrete cover in mm.

    Raises:
        ValueError: If input values are invalid or the combination of exposure
            class, XDRS class and design service life is incompatible.
    """
    if delta_c_dur_gamma < 0:
        raise ValueError(
            f'delta_c_dur_gamma must be positive. Got {delta_c_dur_gamma}'
        )

    # Validate input
    design_service_life = max(50, min(design_service_life, 100))

    # Minimum covers from Table 6.4
    covers = {
        'XRDS 0.5': {
            'XS1': [20, 20],
            'XS2': [20, 30],
            'XS3': [30, 40],
            'XD1': [20, 20],
            'XD2': [20, 30],
            'XD3': [30, 40],
        },
        'XRDS 1': {
            'XS1': [20, 25],
            'XS2': [25, 35],
            'XS3': [35, 45],
            'XD1': [20, 25],
            'XD2': [25, 35],
            'XD3': [35, 45],
        },
        'XRDS 1.5': {
            'XS1': [25, 30],
            'XS2': [30, 40],
            'XS3': [40, 50],
            'XD1': [25, 30],
            'XD2': [30, 40],
            'XD3': [40, 50],
        },
        'XRDS 2': {
            'XS1': [25, 30],
            'XS2': [35, 45],
            'XS3': [45, 55],
            'XD1': [25, 30],
            'XD2': [35, 45],
            'XD3': [45, 55],
        },
        'XRDS 3': {
            'XS1': [30, 35],
            'XS2': [40, 50],
            'XS3': [55, 65],
            'XD1': [30, 35],
            'XD2': [40, 50],
            'XD3': [55, 65],
        },
        'XRDS 4': {
            'XS1': [30, 40],
            'XS2': [50, 60],
            'XS3': [60, 80],
            'XD1': [30, 40],
            'XD2': [50, 60],
            'XD3': [60, 80],
        },
        'XRDS 5': {
            'XS1': [35, 45],
            'XS2': [60, 70],
            'XS3': [70, float('nan')],
            'XD1': [35, 45],
            'XD2': [60, 70],
            'XD3': [70, float('nan')],
        },
        'XRDS 6': {
            'XS1': [40, 50],
            'XS2': [65, 80],
            'XS3': [float('nan'), float('nan')],
            'XD1': [40, 50],
            'XD2': [65, 80],
            'XD3': [float('nan'), float('nan')],
        },
        'XRDS 8': {
            'XS1': [45, 55],
            'XS2': [75, float('nan')],
            'XS3': [float('nan'), float('nan')],
            'XD1': [45, 55],
            'XD2': [75, float('nan')],
            'XD3': [float('nan'), float('nan')],
        },
        'XRDS 10': {
            'XS1': [50, 65],
            'XS2': [80, float('nan')],
            'XS3': [float('nan'), float('nan')],
            'XD1': [50, 65],
            'XD2': [80, float('nan')],
            'XD3': [float('nan'), float('nan')],
        },
    }

    # Determine base cover
    base_cover_v = covers[exposure_resistance_class][exposure_class]

    if not base_cover_v:  # Check for unavailable cover
        raise ValueError(
            'No cover available for the selected exposure'
            + ' resistance class and exposure class combination.'
        )
    base_cover = interp1d(
        [50, 100],
        base_cover_v,
    )(design_service_life)
    if math.isnan(base_cover):
        raise ValueError('Not a valid cover value found.')

    # Calculate minimum cover
    c_min_dur = base_cover + delta_c_dur_gamma
    return max(c_min_dur, 0)


def delta_c_min_30() -> float:
    """Default reduction of minimum cover for structures with design life of 30
    years or less unless specified in National Annex.

    EN1992-1-1:2023 (6.5.2.2(2)).

    Returns:
        float: Reduction distance in mm.
    """
    return -5.0


def delta_c_min_exc() -> float:
    """Default reduction of minimum cover for superior compaction or improved
    curing unless specified in National Annex.

    EN1992-1-1:2023 (6.5.2.2(3)).

    Returns:
        float: Reduction distance in mm.
    """
    return -5.0


def delta_c_min_p() -> float:
    """Default addition of minimum cover for prestressing tendons unless
    specified in National Annex.

    EN1992-1-1:2023 Table. (6.5.2.2(4)).

    Returns:
        float: Added distance in mm.
    """
    return 10.0


def delta_dur_red_1() -> float:
    """Default reduction of minimum cover for use of special measures of
    reinforcing steel unless specified in National Annex.

    EN1992-1-1:2023 (6.5.2.2(5)).

    Returns:
        float: Reduction distance in mm.
    """
    return -10.0


def delta_dur_red_2() -> float:
    """Default reduction of minimum cover for use of special measures of
    reinforcing steel unless specified in National Annex.

    EN1992-1-1:2023 (6.5.2.2(9)).

    Returns:
        float: Reduction distance in mm.
    """
    return -0.0


def delta_dur_abr(xm_class: Literal['XM1', 'XM2', 'XM3']) -> float:
    """Default addition of minimum cover for abraion unless specified in
    National Annex.

    EN1992-1-1:2023 (6.5.2.2(6)).

    Args:
        xm_class (str): The XM class.

    Returns:
        float: Added distance in mm.
    """
    data = {
        'XM1': 5.0,
        'XM2': 10.0,
        'XM3': 15.0,
    }
    return data[xm_class.upper()]
