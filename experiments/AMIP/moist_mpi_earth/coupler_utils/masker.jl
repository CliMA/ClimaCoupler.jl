# TODO: build abstraction
# abstract type AbstractMaskType end
# struct LandSeaMask{I, FT, D} <: AbstractMaskType
#     infile:: String
#     ne:: I 
#     Nq:: I 
#     R:: FT 
#     data:: D
# end

function LandSeaMask(FT, infile, varname, h_space; outfile =  "data_cc.nc")
    R = h_space.topology.mesh.domain.radius
    ne = h_space.topology.mesh.ne
    Nq = Spaces.Quadratures.polynomial_degree(h_space.quadrature_style) + 1

    mask = ncreader_rll_to_cgll(FT, infile,  varname, ne = ne, R = R, Nq = Nq)    
    mask = clean_mask.(FT, mask) 
    # LandSeaMask(infile, ne, Nq, R, mask)
end

"""
clean_mask(FT, mask)
- convert to integer values after interpolation (but keep type as floats foe easier calculation (TODO))
"""
clean_mask(FT, mask) = mask > FT(0.7) ? FT(1) : FT(0)

"""
combine_surface(mask, sfc_1, sfc_2)
- combine two masked surfaces on the same horizontal topology (TODO: generalize to more surfaces and different resolutions)
"""
combine_surface(mask, sfc_1, sfc_2, value = 0.5) = (mask > FT(value) ? sfc_1 : FT(0)) + (mask <= FT(value) ? sfc_2 : FT(0)) 
combine_surface(mask, sfc_1, sfc_2, sfc_3, value1 = -0.5, value2 = 0.5) = (mask < FT(value1) ? sfc_3 : FT(0)) + ((mask >= FT(value1) && (mask <= FT(value2))) ? sfc_2 : FT(0)) + (mask > FT(value2) ? sfc_1 : FT(0)) 

"""
apply_mask(mask, condition, yes, no, value = 0.5) 
- apply mask mased on a threshold value in the mask
"""
apply_mask(mask, condition, yes, no, value = 0.5) = condition(mask, value) ? yes : no