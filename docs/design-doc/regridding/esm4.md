# **ESM4 (FMS) Coupler**

-  https://github.com/NOAA-GFDL/FMScoupler
- [architecture intro](https://extranet.gfdl.noaa.gov/~vb/talks/coupler.pdf)

# Regridding
- modes
    - `simple` (AMIP) - sequential coupler - all exchange fluxes / fields use first order (or, if specified, second order) conservative remap
    - `full` (CMIP) - concurrent - regridding possible

## Implementation

1. calculate areas using [mosaic](https://data1.gfdl.noaa.gov/~arl/pubrel/r/mom4p1/src/mom4p1/doc/mosaic_tool.html) or regular grid
    - vectors are transformed onto an [A grid](http://indico.ictp.it/event/a12235/material/0/2.pdf)

2. [generate xmap](https://github.com/NOAA-GFDL/FMS/blob/main/exchange/xgrid.F90#L1512):

        setup_xmap(xmap, grid_ids, grid_domains, grid_file, atm_grid, lnd_ug_domain)

3. put source_var onto xgrid:

        interface put_to_xgrid

        !> Scatters data to exchange grid
        subroutine put_side1_to_xgrid(d, grid_id, x, xmap, remap_method, complete)
        real(r8_kind), dimension(:,:), intent(in   )    :: d !< data to send
        character(len=3),     intent(in   )    :: grid_id !< 3 character grid ID
        real(r8_kind), dimension(:),   intent(inout)    :: x !< xgrid data
        type (xmap_type),     intent(inout)    :: xmap !< exchange grid
        integer, intent(in), optional          :: remap_method !< exchange grid interpolation method can be FIRST_ORDER(=1) or SECOND_ORDER(=2) 


        put_1_to_xgrid_order_1(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)

- OR if the exchange grid interpolation method remap_method = SECOND_ORDER

        put_1_to_xgrid_order_2(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)
            ! - uses values of the neighboring data points to calculate gradients which will then be passed and add them onto the xgrid

4. apply remapping:
        
        get_from_xgrid
        
        |
        V

            get_side1_from_xgrid(d, grid_id, x, xmap, complete)

            real(r8_kind), dimension(:,:), intent(  out) :: d !< recieved xgrid data
            character(len=3),     intent(in   ) :: grid_id !< 3 character grid ID
            real(r8_kind), dimension(:),   intent(in   ) :: x !< xgrid data
            type (xmap_type),     intent(inout) :: xmap !< exchange grid
            logical, intent(in), optional     :: complete
                |
                V
                get_1_from_xgrid(d_addrs, x_addrs, xmap, isize, jsize, xsize, lsize)

        |
        V

            get_side2_from_xgrid(d, grid_id, x, xmap)
                |
                V
                get_2_from_xgrid(d, grid, x, xmap)


## Refs
- [Balaji et al 06-18 - gridspec](https://arxiv.org/pdf/1911.08638.pdf)



