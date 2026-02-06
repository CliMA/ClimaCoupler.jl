# compute_spectrum for ClimaAnalysis.OutputVar using ClimaCoreSpectra
# TODO: Generalize using ClimaAnalysis

"""
    compute_spectrum(var::ClimaAnalysis.OutputVar; mass_weight = nothing)

Compute the spherical power spectrum associated to the given variable.
Returns a `ClimaAnalysis.OutputVar` with dimensions (Log10 Spherical Wavenumber, z) or
(time, Log10 Spherical Wavenumber, z).

Uses `ClimaCoreSpectra.power_spectrum_2d`, matching the usage in ClimaAtmos.jl
[ci_plots.jl](https://github.com/CliMA/ClimaAtmos.jl/blob/main/post_processing/ci_plots.jl).

The input `var` must have spatial dimensions (lon or long, lat, z) with `len(lon) == 2 * len(lat)`.
Optionally a time dimension may be present. It is advisable to window the OutputVar to a
narrow set of times before passing it to this function. With this function, you can also
take time-averages of spectra (e.g. via `ClimaAnalysis.average_time`).

# Arguments
- `var`: A ClimaAnalysis.OutputVar with 3D (lon, lat, z) or 4D (time, lon, lat, z) data.
- `mass_weight`: Optional vector of length `length(var.dims["z"])` for mass-weighted spectra.
  Defaults to ones.

# Returns
- A ClimaAnalysis.OutputVar with log10(spectrum) and attributes/short_name prefixed by "log_spectrum_".
"""
function compute_spectrum(var::ClimaAnalysis.OutputVar; mass_weight = nothing)
    # power_spectrum_2d works only when the two horizontal dimensions have precisely
    # one twice as many points as the other (len1 == 2*len2)
    if "time" in keys(var.dims)
        time, dim1, dim2, dim3 = var.index2dim[1:4]
        times = var.dims[time]
    else
        dim1, dim2, dim3 = var.index2dim[1:3]
        times = []
    end

    len1 = length(var.dims[dim1])
    len2 = length(var.dims[dim2])

    len1 == 2 * len2 || error("Cannot take this spectrum ($len1 != 2 $len2)")

    (dim1 == "lon" || dim1 == "long") ||
        error("First dimension has to be longitude (found $dim1)")
    dim2 == "lat" || error("Second dimension has to be latitude (found $dim2)")
    (dim3 == "z") || dim3 == "z_reference" || error("Third dimension has to be altitude (found $dim3)")

    FT = eltype(var.data)

    mass_weight =
        isnothing(mass_weight) ? ones(FT, length(var.dims[dim3])) : mass_weight

    # Number of spherical wave numbers, excluding the first and the last
    # (reverse-engineered from ClimaCoreSpectra)
    num_output = Int((floor((2 * len1 - 1) / 3) + 1) / 2 - 1)

    mesh_info = nothing

    if !isempty(times)
        output_spectrum =
            zeros(FT, (length(times), num_output, length(var.dims[dim3])))
        dims = Dict(time => times)
        dim_attributes = Dict(time => var.dim_attributes[time])
        for index in 1:length(times)
            spectrum_data, _wave_numbers, _spherical, mesh_info =
                power_spectrum_2d(FT, var.data[index, :, :, :], mass_weight)
            output_spectrum[index, :, :] .=
                dropdims(sum(spectrum_data; dims = 1); dims = 1)[
                    (begin + 1):(end - 1),
                    :,
                ]
        end
    else
        dims = Dict{String, Vector{FT}}()
        dim_attributes = Dict{String, Dict{String, String}}()
        output_spectrum = zeros(FT, (num_output, length(var.dims[dim3])))
        spectrum_data, _wave_numbers, _spherical, mesh_info =
            power_spectrum_2d(FT, var.data[:, :, :], mass_weight)
        output_spectrum[:, :] .=
            dropdims(sum(spectrum_data; dims = 1); dims = 1)[
                (begin + 1):(end - 1),
                :,
            ]
    end

    w_numbers = collect(1:1:(mesh_info.num_spherical - 1))

    dims["Log10 Spherical Wavenumber"] = log10.(w_numbers)
    dims[dim3] = var.dims[dim3]
    dim_attributes["Log10 Spherical Wavenumber"] = Dict("units" => "")
    dim_attributes[dim3] = var.dim_attributes[dim3]

    attributes = Dict(
        "short_name" => "log_spectrum_" * var.attributes["short_name"],
        "long_name" => "Spectrum of " * var.attributes["long_name"],
        "units" => "",
    )

    return ClimaAnalysis.OutputVar(
        attributes,
        dims,
        dim_attributes,
        log10.(output_spectrum),
    )
end
