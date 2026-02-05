const ConservativeRegriddingCCExt =
    Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)

### Extensions of Interfacer.jl functions for Oceananigans fields/grids
# Non-allocating ClimaCore -> Oceananigans remap
function Interfacer.remap!(target_field::OC.Field, source_field::CC.Fields.Field, remapping)
    ConservativeRegriddingCCExt.get_value_per_element!(
        remapping.value_per_element_cc,
        source_field,
        remapping.field_ones_cc,
    )

    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(target_field, 3)
    dst = vec(OC.interior(target_field, :, :, z))
    src = remapping.value_per_element_cc

    # Regrid the source field to the target field
    CR.regrid!(dst, transpose(remapping.remapper_oc_to_cc), src)
    return nothing
end
# Allocating ClimaCore -> Oceananigans remap
function Interfacer.remap(
    target_space::Union{
        OC.OrthogonalSphericalShellGrid,
        OC.ImmersedBoundaryGrid,
        OC.LatitudeLongitudeGrid,
    },
    source_field::CC.Fields.Field,
    remapping,
)
    target_field = OC.Field{OC.Center, OC.Center, Nothing}(target_space)
    Interfacer.remap!(target_field, source_field, remapping)
    return target_field
end

# Non-allocating Oceananigans Field -> ClimaCore remap
function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field, remapping)
    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(source_field, 3)
    src = vec(OC.interior(source_field, :, :, z))

    # Store the remapped FV values in a vector of length equal to the number of elements in the target space
    dst = remapping.value_per_element_cc

    # Regrid the source field to the target field
    CR.regrid!(dst, remapping.remapper_oc_to_cc, src)

    # Convert the vector of remapped values to a ClimaCore Field with one value per element
    ConservativeRegriddingCCExt.set_value_per_element!(target_field, dst)
    return nothing
end
# Allocating Oceananigans Field -> ClimaCore remap
function Interfacer.remap(
    target_space::CC.Spaces.AbstractSpace,
    source_field::OC.Field,
    remapping,
)
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, source_field, remapping)
    return target_field
end

# Non-allocating Oceananigans operation -> ClimaCore remap
function Interfacer.remap!(
    target_field::CC.Fields.Field,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    Interfacer.remap!(target_field, evaluated_field, remapping)
    return nothing
end
# Allocating Oceananigans operation -> ClimaCore remap
function Interfacer.remap(
    target_space::CC.Spaces.AbstractSpace,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, operation, remapping)
    return target_field
end

# Handle the case of remapping a scalar number to a ClimaCore space
Interfacer.remap!(target_field::CC.Fields.Field, source_field::Number, remapping) =
    Interfacer.remap!(target_field, source_field)
Interfacer.remap(target_space::CC.Spaces.AbstractSpace, source_num::Number, remapping) =
    Interfacer.remap(target_space, source_num, remapping)

# Handle the case of remapping the area fraction field, which is a ClimaCore Field
Interfacer.remap!(target_field::CC.Fields.Field, source_field::CC.Fields.Field, remapping) =
    Interfacer.remap!(target_field, source_field)

# Extend Interfacer.get_field to allow automatic remapping to the target space
function Interfacer.get_field!(target_field, sim::Union{OceananigansSimulation, ClimaSeaIceSimulation}, quantity)
    Interfacer.remap!(target_field, Interfacer.get_field(sim, quantity), sim.remapping)
    return nothing
end
# TODO see if we can remove this allocating version
function Interfacer.get_field(
    target_space::CC.Spaces.AbstractSpace,
    sim::OceananigansSimulation,
    quantity,
)
    return Interfacer.remap(
        target_space,
        Interfacer.get_field(sim, quantity),
        sim.remapping,
    )
end
