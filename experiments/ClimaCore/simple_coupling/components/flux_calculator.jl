# - Surface Flux Calculation (coarse bulk formula)
calculate_flux(T_sfc, T1, parameters) = -parameters.λ * (T_sfc - T1);
