"""Example code to calculate the velocity_of_corrosion and the area_after_corrosion of rebars."""

from structuralcodes.codes.mc2020._corrosion import calculate_velocity_of_corrosion,calculate_minimum_area_after_corrosion

# Calculate the representative velocity of corrosion (defined Pcorr_rep according to MC2020)
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="carbonation_induced",exposure_class="Unsheltered")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="carbonation_induced",exposure_class="Sheltered")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Wet")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Airborn_seawater")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Submerged")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Tidal_zone")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Cyclic_dry_wet")
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Cyclic_dry_wet",fractile=0.5) #should give 30
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Cyclic_dry_wet",fractile=0.8413) #should give 30+1*40=70 (1 stdev)
Pcorr_rep=calculate_velocity_of_corrosion(corrosion_type="chloride_induced",exposure_class="Cyclic_dry_wet",fractile=0.9772) #should give 30+2*40=110 (2 stdev)
print("Calculated velocity of corrosion Pcorr_rep = "+str(round(Pcorr_rep,2))+" μm/yr")

# Calculate the minimum area after corrosion of a rebar of diameter 16mm.
InitialArea=8*8*3.141592
print("Area before corrosion = "+str(round(InitialArea,2))+" mm2")
# MethodA: by indicating the mass_loss
Minimum_area_after_corrosion=calculate_minimum_area_after_corrosion(uncorroded_area=InitialArea,pitting_factor=1.2,mass_loss=0.3)
# MethodB: by indicating the velocity_of_corrosion and the time_of_corrosion
Minimum_area_after_corrosion=calculate_minimum_area_after_corrosion(uncorroded_area=InitialArea,pitting_factor=1,velocity_of_corrosion=100,time_of_corrosion=10)

print("Calculated Minimum_area_after_corrosion = "+str(round(Minimum_area_after_corrosion,2))+" mm2")



