import typing as t
from structuralcodes.codes import _use_design_code

def calculate_velocity_of_corrosion(
    corrosion_type: t.Literal["carbonation_induced","chloride_induced"],
    exposure_class: t.Literal["Sheltered","Unsheltered","Wet","Cyclic_dry_wet","Airborn_seawater","Submerged","Tidal_zone"],
    fractile: t.Optional[float] = 0.5,
) -> float:
    """A function to calculate the representative velocity of corrosion given certain environmental conditions.
    Actually (02/12/2025) only MC2020 is implemented, see table 30.1-6a and 30.1-6b.

    Keyword Arguments:
        corrosion_type (str): Corrosion type, which can be "carbonation_induced" or "chloride_induced" (default: "chloride_induced").
        exposure_class (str): Exposure class, which can be "Sheltered" or "Unsheltered" for corrosion_type=="carbonation_induced"
            and "Wet", "Cyclic_dry_wet", "Airborn_seawater", "Submerged", "Tidal_zone" for corrosion_type=="chloride_induced"
            (default: "Unsheltered" for corrosion_type=="carbonation_induced" and "Cyclic_dry_wet" for corrosion_type=="chloride_induced").
        fractile (float): A statistic parameter used to obtain the velocity_of_corrosion at a given number of standard deviation from the
            average value. For example, the velocity of corrosion calculacted with fractile=0.5 can be interpreted as "50% of rebars
            have a velocity of corrosion lower than Pcorr_rep (the returning value of this function). For example, the velocity of corrosion
            calculacted with fractile=0.99 can be interpreted as "99% of rebars have a velocity of corrosion lower than Pcorr_rep, giving a 
            higher degree of safety of the returned value Pcorr_rep. (default: "0.5")

    Return:
        Pcorr_rep (float): representative velocity of corrosion of the rebars, defined as the ratio (asbolute value) between the variation of the rebar's diameter
            and the time of exposure. Unit of measurement: μm/yr.

    Raises:
        ValueError: if the corrosion type is not valid.
        ValueError: if the exposure class is not valid for the relative corrosion type.
        ValueError: if the fractile value is not valid.
        ValueError: if the fractile value used causes Pcorr_rep to obtain values lower than 0.
    """

    # Check if corrosion_type is valid.
    ValidCorrosionTypes=["carbonation_induced","chloride_induced"]
    if corrosion_type not in ValidCorrosionTypes:
        raise ValueError(
            'The variable corrosion_type is not valid, either use ' \
            '"carbonation_induced" or "chloride_induced"'
        )
    
    # Set exposure_class to default values if None. Then, check if exposure_class is valid.
    ValidExposureClassesForCarbonationInducedCorrosion=["Sheltered","Unsheltered"]
    ValidExposureClassesForChlorideInducedCorrosion=["Wet","Cyclic_dry_wet","Airborn_seawater","Submerged","Tidal_zone"]
    if corrosion_type=="carbonation_induced" and exposure_class not in ValidExposureClassesForCarbonationInducedCorrosion:
        raise ValueError(
            'The variable exposure_class is not valid for the current corrosion_type, either use ' \
            '"Sheltered" or "Unsheltered"'
        )
    if corrosion_type=="chloride_induced" and exposure_class not in ValidExposureClassesForChlorideInducedCorrosion:
        raise ValueError(
            'The variable exposure_class is not valid for the current corrosion_type, either use ' \
            '"Wet", "Cyclic_dry_wet", "Airborn_seawater", "Submerged" or "Tidal_zone"'
        )
    
    # Check if fractile value is valid.
    if fractile <=0 or fractile >=1:
        raise ValueError(
            'The value of fractile is not valid, use a value between 0 and 1 (both excluded)'
        )

    # Calculate MeanValue and StandardDeviation
    if corrosion_type == "carbonation_induced":
        table_mean = {
            "Sheltered": 2,
            "Unsheltered": 5,
        }
        table_sd = {
            "Sheltered": 3,
            "Unsheltered": 1,   # the table shows "t", but assuming 1 μm/yr (fix if needed)
        }
    elif corrosion_type == "chloride_induced":
        table_mean = {
            "Wet": 4,
            "Cyclic_dry_wet": 30,
            "Airborn_seawater": 30,
            "Submerged": 4,
            "Tidal_zone": 50,
        }
        table_sd = {
            "Wet": 6,
            "Cyclic_dry_wet": 40,
            "Airborn_seawater": 40,
            "Submerged": 7,
            "Tidal_zone": 100,
        }
    mean_value = table_mean[exposure_class]
    sd_value = table_sd[exposure_class]

    # Calculate Pcorr_rep for the correct fractile value (use Acklam's algorithm for the approximation of the inverse
    # standard normal CDF in order to avoid to add other dependencies).
    def inverse_normal_cdf(p: float) -> float:
        """Approximation of the inverse standard normal CDF (Acklam's algorithm)."""
        if p <= 0.0 or p >= 1.0:
            raise ValueError("p must be in (0,1)")
        
        import math
        from math import  sqrt

        # Coefficients for the approximation
        a = [ -3.969683028665376e+01,
            2.209460984245205e+02,
            -2.759285104469687e+02,
            1.383577518672690e+02,
            -3.066479806614716e+01,
            2.506628277459239e+00 ]

        b = [ -5.447609879822406e+01,
            1.615858368580409e+02,
            -1.556989798598866e+02,
            6.680131188771972e+01,
            -1.328068155288572e+01 ]

        c = [ -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
            4.374664141464968e+00,
            2.938163982698783e+00 ]

        d = [ 7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996e+00,
            3.754408661907416e+00 ]

        # Define break-points
        plow  = 0.02425
        phigh = 1 - plow

        if p < plow:
            q = sqrt(-2 * math.log(p))
            return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
        elif p > phigh:
            q = sqrt(-2 * math.log(1 - p))
            return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                    ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
        else:
            q = p - 0.5
            r = q * q
            return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
                (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
        
    z = inverse_normal_cdf(fractile)
    Pcorr_rep = mean_value + z * sd_value

    # Validate the result, which could be <=0 if fractile is too low.
    if Pcorr_rep<=0:
        raise ValueError(
            'The calculated value of Pcorr_rep is less or equal to 0 because the value of fractile is too low. '\
            'Try using a higher value of fractile. Note that using a low value of fractile could lead to an '\
            'underestimention of the effect of corrosion.'
        )

    return Pcorr_rep

def calculate_minimum_area_after_corrosion(
    uncorroded_area: float,
    pitting_factor: t.Optional[float] = 1,
    mass_loss: t.Optional[float] = None,
    velocity_of_corrosion: t.Optional[float] = None,
    time_of_corrosion: t.Optional[float] = None,
) -> float:
    
    """A function to calculate the minimum residual steel area after corrosion.
    The function allows two alternative approaches: (i) using the total mass loss,
    or (ii) using the velocity of corrosion together with the time of corrosion.
    A pitting factor is used to relate the maximum pit depth to the average pit depth.
    Actually (02/12/2025) only MC2020 is implemented, see chapter 30.1.11.3.5.
    Function is valid for round bars!.

    Keyword Arguments:
        uncorroded_area (float): Original (uncorroded) cross-sectional area of the rebar [mm²].
        pitting_factor (float): Ratio between maximum pit depth and average pit depth
            (default: 1). Must be ≥ 1.
        mass_loss (float): Fractional mass loss (0–1). If provided, velocity_of_corrosion
            and time_of_corrosion must be None. (default: None)
        velocity_of_corrosion (float): Corrosion rate [μm/yr]. Must be ≥ 0. Used together
            with time_of_corrosion. If provided, mass_loss must be None. (default: None)
        time_of_corrosion (float): Time of corrosion [years]. Must be ≥ 0. Used with
            velocity_of_corrosion. (default: None)

    Return:
        corroded_minimum_area (float): Minimum residual area of the corroded rebar [mm²],
            computed assuming axisymmetric corrosion with maximum pit depth defined by
            `pitting_factor * average_pit_depth`.

    Raises:
        ValueError: if both mass_loss and (velocity_of_corrosion + time_of_corrosion)
            are provided simultaneously.
        ValueError: if only one between velocity_of_corrosion and time_of_corrosion is provided.
        ValueError: if pitting_factor < 1.
        ValueError: if uncorroded_area < 0.
        ValueError: if mass_loss is not in the range [0, 1].
        ValueError: if velocity_of_corrosion < 0 or time_of_corrosion < 0.
        ValueError: if the computed minimum residual radius becomes negative.
    """

    import math
    from math import  sqrt
    uncorroded_diameter=sqrt(uncorroded_area/math.pi)*2
    
    # Input data validation
    if (mass_loss is not None) and (velocity_of_corrosion is not None or time_of_corrosion is not None):
        raise ValueError(
            'Too many input arguments have been given to the function. '
            'Either use mass_loss or velocity_of_corrosion + time_of_corrosion.'
        )

    if (velocity_of_corrosion is not None and time_of_corrosion is None) or \
       (velocity_of_corrosion is None and time_of_corrosion is not None):
        raise ValueError(
            'velocity_of_corrosion or time_of_corrosion is missing.'
        )

    if pitting_factor < 1:
        raise ValueError(
            'The pitting_factor must be greater than or equal to 1.'
        )

    if uncorroded_area < 0:
        raise ValueError(
            'The uncorroded_area must be non-negative.'
        )

    if mass_loss is not None and (mass_loss < 0 or mass_loss > 1):
        raise ValueError(
            'mass_loss must be between 0 and 1.'
        )

    if velocity_of_corrosion is not None and (velocity_of_corrosion < 0):
        raise ValueError(
            'velocity_of_corrosion must be non-negative.'
        )

    if time_of_corrosion is not None and (time_of_corrosion < 0):
        raise ValueError(
            'time_of_corrosion must be non-negative.'
        )
    
    # Calculate corroded area if mass_loss is given
    if mass_loss is not None:
        corroded_average_area=uncorroded_area*(1-mass_loss)
        corroded_average_pit=uncorroded_diameter/2-sqrt(corroded_average_area/math.pi)
        corroded_maximum_pit=pitting_factor*corroded_average_pit
        corroded_minimum_area=(uncorroded_diameter/2-corroded_maximum_pit)**2*math.pi
        if (uncorroded_diameter/2-corroded_maximum_pit) < 0:
            raise ValueError(
                'The combination of mass_loss and pitting_factor gave a corroded_minimum_area'
                'lower than 0. Consider changing these values.'
            )
        return corroded_minimum_area
    
    # Calculate corroded area if time_of_corrosion and velocity_of_corrosion are given
    if velocity_of_corrosion is not None and time_of_corrosion is not None:
        velocity_of_corrosion=velocity_of_corrosion/1000 #convert from micrometers/year to millimeters/year
        corroded_average_pit=velocity_of_corrosion*time_of_corrosion
        if corroded_average_pit>uncorroded_diameter:
            raise ValueError(
                'The combination of velocity_of_corrosion and time_of_corrosion gave a corroded_minimum_area'
                'lower than 0. Consider changing these values.'
            )
        corroded_maximum_pit=pitting_factor*corroded_average_pit
        corroded_minimum_area=(uncorroded_diameter/2-corroded_maximum_pit)**2*math.pi
        if uncorroded_diameter/2-corroded_maximum_pit < 0:
            raise ValueError(
                'The combination of velocity_of_corrosion and time_of_corrosion and pitting_factor'
                'gave a corroded_minimum_area lower than 0. Consider changing these values.'
            )
        return corroded_minimum_area
    return


