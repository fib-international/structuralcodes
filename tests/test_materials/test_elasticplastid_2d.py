import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws import ElasticPlastic2D


@pytest.mark.parametrize(
    'E, fy, strain',
    [
        (210000, 410, np.array([0.001, 0.0, 0.0])),
        (210000, 410, np.array([-0.002, -0.1, -0.002])),
        (200000, 450, np.array([0.003, 0.005, 0.010])),
    ],
)
def test_elasticplastic_2d_get_stress(E, fy, strain):
    """Test elasticplastic 2D material."""
    mat = ElasticPlastic2D(E, fy)
    actual = mat.get_stress(strain)

    expected = np.array(
        [
            np.clip(E * strain[0], -fy, +fy),
            np.clip(E * strain[1], -fy, +fy),
            0.0,
        ]
    )

    assert np.allclose(actual, expected, atol=1e-8)


def test_elasticplastic_2d_input_correct():
    """Test invalid input values for ElasticPlastic2D."""
    with pytest.raises(ValueError) as excinfo:
        ElasticPlastic2D(-210000, 500)
    assert str(excinfo.value) == 'Elastic modulus E must be greater than zero'


@pytest.mark.parametrize(
    'E, fy',
    [
        (210000, 500),
        (200000, 500),
        (195000, 500),
    ],
)
def test_elasticplastic_get_secant(E, fy):
    """Test the elasticPlastic2D material."""
    assert np.allclose(
        ElasticPlastic2D(E, fy).get_secant(),
        np.array([[E, 0, 0], [0, E, 0], [0, 0, 0]]),
    )
