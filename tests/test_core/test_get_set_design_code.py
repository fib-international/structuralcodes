"""Tests for getting and setting design codes."""

import types

import pytest

import structuralcodes


@pytest.mark.parametrize(
    'design_code_to_set',
    ['mc2010', 'MC2010', 'mC2010'],
)
def test_set_design_code(design_code_to_set):
    """Test setting the design code."""
    # Arrange
    expected_design_code_title = 'fib Model Code 2010'

    # Act
    structuralcodes.set_design_code(design_code_to_set)

    # Assert
    assert isinstance(structuralcodes.codes._CODE, types.ModuleType)
    assert structuralcodes.codes._CODE.__title__ == expected_design_code_title


def test_get_design_codes():
    """Test get a list of implemented design codes."""
    # Arrange
    expected_list_of_codes = list(structuralcodes.codes._DESIGN_CODES.keys())

    # Act
    available_codes = structuralcodes.get_design_codes()

    # Assert
    assert available_codes == expected_list_of_codes


@pytest.mark.parametrize(
    'na_to_set',
    ['NO', 'no', 'No', 'nO'],
)
def test_set_national_annex(na_to_set):
    """Test setting the national annex."""
    # Arrange
    expected_na = 'no'

    # Act
    structuralcodes.set_national_annex(na_to_set)

    # Assert
    assert expected_na == structuralcodes.codes._NATIONAL_ANNEX


@pytest.mark.parametrize(
    'design_code_to_user',
    ['mc2010', 'MC2010', 'mC2010'],
)
def test_use_design_code(design_code_to_user):
    """Test return a design code for use."""
    # Arrange
    expected_design_code_title = 'fib Model Code 2010'

    # Act
    code_to_use = structuralcodes.codes._use_design_code(design_code_to_user)

    # Assert
    assert isinstance(code_to_use, types.ModuleType)
    assert code_to_use.__title__ == expected_design_code_title


def test_use_design_code_none():
    """Test return _CODE if called with None."""
    # Arrange
    design_code_to_set = 'mc2010'
    expected_design_code_title = 'fib Model Code 2010'
    structuralcodes.set_design_code(design_code_to_set)

    # Act
    code_to_use = structuralcodes.codes._use_design_code()

    # Assert
    assert isinstance(code_to_use, types.ModuleType)
    assert code_to_use.__title__ == expected_design_code_title
