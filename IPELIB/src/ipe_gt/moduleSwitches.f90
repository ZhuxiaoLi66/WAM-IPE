MODULE moduleSwitches

IMPLICIT NONE

PRIVATE

PUBLIC :: SwitchesType

TYPE SwitchesType
  logical :: electro
  logical :: External_model_provides_NO_N4S_densities
  logical :: External_model_provides_low_lat_E_fields
  logical :: input_Auroral_production_is_single_overall_rate
END TYPE SwitchesType


END MODULE moduleSwitches
