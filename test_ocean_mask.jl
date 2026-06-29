"""
Standalone test script for ocean_fraction masking BFS.
Reproduces the _init_ocean_isolation_mask! logic with visualization.
"""

using ClimaCore
using ClimaComms
using Plots

const CC = ClimaCore
gr()

function test_ocean_mask(h_elem=8, show_diagnostics=true)
    # Create a minimal cubed-sphere space (radius in meters)
    horz_space = CC.CommonSpaces.CubedSphereSpace(
        Float64;
        radius = 6.371e6,
        n_quad_points = 4,
        h_elem = h_elem,
    )

    # Realistic synthetic ocean_fraction: binary (~1 over water, ~0 over land), as in
    # the actual remapped data. Geography:
    #   - Open ocean everywhere by default (val = 1.0)
    #   - Antarctica: lat < -70°  → 0.0
    #   - Greenland-ish polar landmass: lat > 75° AND lon ∈ [-60°, -20°] → 0.0
    #   - Caspian-style isolated inland lake at lat≈42°, lon≈50° (radius 4°) → 0.6
    #     (this MUST be excluded by the mask; it's not connected to the global ocean)
    coords = CC.Fields.coordinate_field(horz_space)
    lats_field = coords.lat
    lons_field = coords.long
    ocean_frac_initial = CC.Fields.map(lats_field, lons_field) do lat, lon
        # Antarctica
        if lat < -70.0
            return 0.0
        end
        # Greenland-ish
        if lat > 75.0 && -60.0 <= lon <= -20.0
            return 0.0
        end
        # Isolated inland lake (Caspian-ish)
        dlat = lat - 42.0
        dlon = lon - 50.0
        if dlat^2 + (dlon * cosd(lat))^2 <= 16.0   # ~4° radius
            return 0.6
        end
        return 1.0
    end

    # Extract raw data for masking logic
    local_lats = vec(Array(CC.Fields.field2array(coords.lat)))
    local_lons = vec(Array(CC.Fields.field2array(coords.long)))
    local_vals = vec(Array(CC.Fields.field2array(ocean_frac_initial)))

    @info "Grid info" h_elem n_local=length(local_lats)
    @info "Ocean_frac range before mask" minimum(local_vals) maximum(local_vals)
    @info "Ocean_frac non-zero count" count(x -> x > 0, local_vals)

    # Compute the mask using the same logic as _init_ocean_isolation_mask!
    n = length(local_vals)
    keep = falses(n)
    perm = sortperm(local_lats)
    sorted_lats = local_lats[perm]

    # Adaptive radius — based on the LARGEST inter-cell lat gap so the BFS can bridge
    # the cubed-sphere face boundary. Mirrors the production code.
    lat_gaps = diff(sorted_lats)
    positive_gaps = filter(g -> g > 1e-8, lat_gaps)
    max_gap = isempty(positive_gaps) ? 1.0 : maximum(positive_gaps)
    r = 1.5 * Float64(max_gap)

    @info "Radius computation (deg)" max_gap r
    @info "Gap distribution (deg)" min_positive=minimum(positive_gaps) max_positive=max_gap

    # Seed BFS (all lat/lon values are in degrees, matching production)
    seed_thr = 0.8
    strong_thr = 0.8
    min_thr = 0.25
    seed_lat_band = 40.0  # degrees

    queue = Int[]
    for i in 1:n
        if Float64(local_vals[i]) > Float64(seed_thr) && abs(Float64(local_lats[i])) <= Float64(seed_lat_band)
            keep[i] = true
            push!(queue, i)
        end
    end

    n_seeds = count(keep)
    @info "BFS seeding" n_seeds_initial=n_seeds

    # Diagnostic: lat structure right at the seed band edge
    border_lats = sort(unique(round.(local_lats[abs.(local_lats) .>= 39.0 .&& abs.(local_lats) .<= 42.0], digits=3)))
    @info "Lat values near ±40° boundary (rounded to 0.001°)" border_lats
    n_outside_strong = count(i -> abs(local_lats[i]) > 40.0 && local_vals[i] > min_thr, 1:n)
    @info "Non-seed cells with val > min_thr (BFS targets)" n_outside_strong

    # BFS expansion
    n_expanded = 0
    while !isempty(queue)
        i = pop!(queue)
        lat_i = Float64(local_lats[i])
        lon_i = Float64(local_lons[i])
        lo = searchsortedfirst(sorted_lats, lat_i - r)
        hi = searchsortedlast(sorted_lats, lat_i + r)

        for k in lo:hi
            j = perm[k]
            keep[j] && continue
            Float64(local_vals[j]) <= Float64(min_thr) && continue
            dlon = abs(Float64(local_lons[j]) - lon_i)
            dlon > 180.0 && (dlon = 360.0 - dlon)
            dlat = Float64(sorted_lats[k]) - lat_i
            if dlat^2 + (dlon * cosd(lat_i))^2 ≤ r^2
                keep[j] = true
                n_expanded += 1
                if Float64(local_vals[j]) >= Float64(strong_thr)
                    push!(queue, j)
                end
            end
        end
    end

    n_final_keep = count(keep)
    @info "BFS expansion" n_expanded n_final_keep

    if show_diagnostics
        # Create output fields for plotting
        keep_field = CC.Fields.map(_ -> 0.0, ocean_frac_initial)
        keep_array = CC.Fields.field2array(keep_field)
        for i in 1:n
            if keep[i]
                keep_array[i] = 1.0
            end
        end

        # Create lat/lon grids for plotting (already in degrees)
        lats_plot = Float64.(local_lats)
        lons_plot = Float64.(local_lons)
        vals_plot = Float64.(local_vals)
        keep_plot = Float64.(keep)

        # Diagnostic plots
        p1 = scatter(lons_plot, lats_plot, marker_z=vals_plot,
                    title="Ocean_frac before mask", ylabel="latitude",
                    markersize=5, markerstrokewidth=0, legend=false, colorbar=true)
        p2 = scatter(lons_plot, lats_plot, marker_z=keep_plot,
                    title="Mask (1=keep, 0=zero)", ylabel="latitude",
                    markersize=5, markerstrokewidth=0, legend=false, colorbar=true)
        p3 = scatter(lons_plot, lats_plot, marker_z=vals_plot .* keep_plot,
                    title="Ocean_frac after mask", ylabel="latitude",
                    markersize=5, markerstrokewidth=0, legend=false, colorbar=true)

        plot(p1, p2, p3, layout=(3,1), size=(1200,1500))
        savefig("/tmp/ocean_mask_test.png")
        @info "Saved plot to /tmp/ocean_mask_test.png"

        # Latitude slice diagnostics (degrees)
        @info "Max latitude in keep=true cells (deg):" maximum(lats_plot[keep_plot .> 0.5])
        @info "Min latitude in keep=true cells (deg):" minimum(lats_plot[keep_plot .> 0.5])

        # Polar-ocean preservation check
        antarctic_kept = count(i -> -70.0 <= local_lats[i] < -60.0 && keep[i], 1:n)
        antarctic_total = count(i -> -70.0 <= local_lats[i] < -60.0 && local_vals[i] > min_thr, 1:n)
        arctic_kept = count(i -> 60.0 <= local_lats[i] <= 90.0 && local_vals[i] > 0.8 && keep[i], 1:n)
        arctic_total = count(i -> 60.0 <= local_lats[i] <= 90.0 && local_vals[i] > 0.8, 1:n)
        @info "Southern Ocean (-70 to -60°) coverage" antarctic_kept antarctic_total
        @info "Arctic ocean (>=60°N, val>0.8) coverage" arctic_kept arctic_total

        # Isolated lake exclusion check (Caspian at lat≈42, lon≈50, r≈4°)
        n_lake_total = count(1:n) do i
            (local_lats[i] - 42.0)^2 + ((local_lons[i] - 50.0) * cosd(local_lats[i]))^2 <= 16.0
        end
        n_lake_kept = count(1:n) do i
            in_lake = (local_lats[i] - 42.0)^2 + ((local_lons[i] - 50.0) * cosd(local_lats[i]))^2 <= 16.0
            in_lake && keep[i]
        end
        @info "Inland lake exclusion (should be 0/$n_lake_total)" n_lake_kept
    end

    return keep, local_vals, local_lats
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_ocean_mask(8, true)
end
