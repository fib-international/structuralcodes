"""tests for the functions in structuralcodes.development.corrosion"""

import pytest

from structuralcodes.codes.mc2020._corrosion import calculate_velocity_of_corrosion

# Should return the velocity of corrosion (average representative value).
@pytest.mark.parametrize(
    "corrosion_type, exposure_class, fractile, expected",
    [
        ("carbonation_induced", "Sheltered", 0.5, 2),
        ("carbonation_induced", "Unsheltered",0.5, 5),
        ("chloride_induced", "Wet",0.5, 4),
        ("chloride_induced", "Cyclic_dry_wet",0.5, 30),
        ("chloride_induced", "Airborn_seawater",0.5, 30),
        ("chloride_induced", "Submerged",0.5, 4),
        ("chloride_induced", "Tidal_zone",0.5, 50),
    ]
)
def test_calculate_velocity_of_corrosion_without_fractile(corrosion_type, exposure_class, fractile, expected):
    Pcorr_rep = calculate_velocity_of_corrosion(corrosion_type=corrosion_type,exposure_class=exposure_class)
    assert Pcorr_rep == expected

# Should return the velocity of corrosion for the relative fractile.
@pytest.mark.parametrize(
    "corrosion_type, exposure_class, fractile, expected",
    [
        ("carbonation_induced", "Sheltered", 0.8413, 2+3),
        ("carbonation_induced", "Unsheltered",0.8413, 5+1),
        ("chloride_induced", "Wet",0.8413, 4+6),
        ("chloride_induced", "Cyclic_dry_wet",0.8413, 30+40),
        ("chloride_induced", "Airborn_seawater",0.8413, 30+40),
        ("chloride_induced", "Submerged",0.8413, 4+7),
        ("chloride_induced", "Tidal_zone",0.8413, 50+100),
        ("carbonation_induced", "Sheltered", 0.5, 2),
        ("carbonation_induced", "Unsheltered",0.5, 5),
        ("chloride_induced", "Wet",0.5, 4),
        ("chloride_induced", "Cyclic_dry_wet",0.5, 30),
        ("chloride_induced", "Airborn_seawater",0.5, 30),
        ("chloride_induced", "Submerged",0.5, 4),
        ("chloride_induced", "Tidal_zone",0.5, 50),
    ]
)
def test_calculate_velocity_of_corrosion_with_fractile(corrosion_type, exposure_class, fractile, expected):
    Pcorr_rep = calculate_velocity_of_corrosion(corrosion_type=corrosion_type,exposure_class=exposure_class,fractile=fractile)
    assert abs(Pcorr_rep - expected) < expected*0.001 #<0.1% error

# Should raise error because wrong combinations of corrosion_type and exposure_class are given as input.
@pytest.mark.parametrize(
    "corrosion_type, exposure_class",
    [
        ("carbonation_induced", "Wet"),
        ("carbonation_induced", "Cyclic_dry_wet"),
        ("carbonation_induced", "Airborn_seawater"),
        ("carbonation_induced", "Submerged"),
        ("carbonation_induced", "Tidal_zone"),
        ("chloride_induced", "Sheltered"),
        ("chloride_induced", "chloride_induced")
    ]
)
def test_wrong_corrosion_type_and_exposure_class_combinations(corrosion_type, exposure_class):
    with pytest.raises(Exception):
        calculate_velocity_of_corrosion(corrosion_type=corrosion_type,exposure_class=exposure_class)

# Should raise error because too low values of fractile gives negative velocity of corrosion.
def test_low_values_of_fractile():
    with pytest.raises(Exception):
        calculate_velocity_of_corrosion(corrosion_type="carbonation_induced",exposure_class="Sheltered",fractile=0.001)

from structuralcodes.codes.mc2020._corrosion import calculate_minimum_area_after_corrosion

# Round bar of diameter=16mm.
InitialArea=8*8*3.141592

# Should return the remaining area after corrosion.
@pytest.mark.parametrize(
    "mass_loss, pitting_factor, expected",
    [
        (0, 1, InitialArea),
        (0, 2, InitialArea),
        (0.5, 1, 0.5*InitialArea),
        (0.5, 2, (2*0.70710678118-1)**2*InitialArea),
        (0.75, 1, 0.25*InitialArea),
        (1, 1, 0),
    ]
)
def test_minimum_area_after_corrosion_given_mass_loss_and_pitting_factor(mass_loss, pitting_factor, expected):
    InitialArea=8*8*3.141592
    Minimum_area_after_corrosion=calculate_minimum_area_after_corrosion(uncorroded_area=InitialArea,pitting_factor=pitting_factor,mass_loss=mass_loss)
    assert abs(Minimum_area_after_corrosion-expected) <= expected *0.00000001

# Should return the remaining area after corrosion.
@pytest.mark.parametrize(
    "velocity_of_corrosion,time_of_corrosion, pitting_factor, expected",
    [
        (0,0, 1, InitialArea),
        (0,0, 2, InitialArea),
        (100,10, 1, (7**2/8**2)*InitialArea),
        (100,10, 2, (6**2/8**2)*InitialArea),
        (100,40, 1, (4**2/8**2)*InitialArea),
    ]
)
def test_minimum_area_after_corrosion_given_velocity_of_corrosion_and_time_of_corrosion(velocity_of_corrosion,time_of_corrosion, pitting_factor, expected):
    InitialArea=8*8*3.141592
    Minimum_area_after_corrosion=calculate_minimum_area_after_corrosion(
        uncorroded_area=InitialArea,pitting_factor=pitting_factor,
        velocity_of_corrosion=velocity_of_corrosion,
        time_of_corrosion=time_of_corrosion)
    assert abs(Minimum_area_after_corrosion-expected) <= expected *0.0001

# Should raise different kind of errors (wrong combinations of input values or negative area after corrosion).
@pytest.mark.parametrize(
    "velocity_of_corrosion, time_of_corrosion, mass_loss, pitting_factor",
    [
        (100,10, 0.5, 1),
        (100,None, 0.5, 1),
        (None,10, 0.5, 1),
        (100,100, None, 1),
        (100,50, None, 2),
        (None,None, 1.5, 1),
        (None,None, 0.9, 2),
    ]
)
def test_wrong_combinations_or_negative_area_after_corrosion(velocity_of_corrosion,time_of_corrosion, pitting_factor, mass_loss):
    InitialArea=8*8*3.141592
    with pytest.raises(Exception):
        Minimum_area_after_corrosion=calculate_minimum_area_after_corrosion(
            mass_loss=mass_loss,
            uncorroded_area=InitialArea,
            pitting_factor=pitting_factor,
            velocity_of_corrosion=velocity_of_corrosion,
            time_of_corrosion=time_of_corrosion)



