# Regime attribution of the clt pattern residual, and the B9 alignment
# test. Zero GPU: existing member output + artifacts.
#
# 1. Residual r = model clt (posterior of record: edmf iteration-1
#    member 1, ClimaAtmos 0.42.3) minus CALIPSO clt (Sep 2006-2010
#    mean), on the common 2.5 degree mask. r_anom = its zonal anomaly.
# 2. Regime predictors from ERA5 monthly means (Sep 2006-2010): LTS =
#    theta(700 hPa) - T(1000 hPa), and omega at 500 hPa. Ocean mask from
#    MAC lwp coverage. Composites: share of the cos-weighted squared
#    residual per predictor decile.
# 3. Regression r ~ 1 + LTS + omega500 (cos-weighted) for explained
#    variance = feasibility of regime-dependent parameters.
# 4. B9 alignment: signed per-parameter response maps from the release
#    and edmf iteration-1 ensembles, projected onto r and r_anom in
#    noise-whitened space. This is the number scoping should gate on.
import NCDatasets, Statistics, Dates, LinearAlgebra

const ART = "/glade/campaign/univ/ucit0011/ClimaArtifacts2/artifacts"
const OBS_FP = joinpath(ART, "calipso_cloudsat", "radarlidar_monthly_2.5x2.5.nc")
const ERA_FP = joinpath(
    ART, "era5_monthly_averages_pressure_levels_1979_2024",
    "era5_monthly_averages_pressure_levels_197901-202410.nc",
)
const MAC_FP = joinpath(ART, "mac_lwp", "mac_lwp.nc")

const ENSEMBLES = [
    (
        name = "release (closure+optics)",
        dir = "/glade/derecho/scratch/nefrathe/amip_calibration_release_out/iteration_001",
        simsub = "amip_calibration_release/output_active/clima_atmos",
    ),
    (
        name = "edmf (EDMF family)",
        dir = "/glade/derecho/scratch/nefrathe/amip_calibration_edmf_out/iteration_001",
        simsub = "amip_calibration_edmf/output_active/clima_atmos",
    ),
]

# ---- CALIPSO clt (fraction) per September, on (lon, lat) -------------------
function load_obs()
    NCDatasets.NCDataset(OBS_FP) do ds
        v = ds["cloud_cover_in_column"]
        dims = collect(NCDatasets.dimnames(v))
        lon = Float64.(Array(ds["lon"][:]))
        lat = Float64.(Array(ds["lat"][:]))
        t = Array(ds["time"][:])
        li = findfirst(==("lon"), dims)
        ji = findfirst(==("lat"), dims)
        fields = Dict{Int, Matrix{Float64}}()
        for (k, tt) in enumerate(t)
            d = Dates.DateTime(tt)
            (Dates.month(d) == 9 && Dates.year(d) in 2006:2010) || continue
            sel = Dict("lon" => Colon(), "lat" => Colon(), "type" => 1, "doop" => 1, "time" => k)
            x = Array(v[(sel[d2] for d2 in dims)...])
            f = Float64.(replace(x, missing => NaN)) ./ 100.0
            li > ji && (f = permutedims(f))
            fields[Dates.year(d)] = f
        end
        return lon, lat, fields
    end
end

function load_member(mdir)
    fp = joinpath(mdir, "clt_1M_average.nc")
    isfile(fp) || return nothing, nothing, nothing
    NCDatasets.NCDataset(fp) do ds
        v = ds["clt"]
        dims = NCDatasets.dimnames(v)
        li = findfirst(d -> startswith(lowercase(d), "lon"), dims)
        ji = findfirst(d -> startswith(lowercase(d), "lat"), dims)
        ti = findfirst(d -> occursin("time", lowercase(d)), dims)
        lon = Float64.(Array(ds[dims[li]][:]))
        lat = Float64.(Array(ds[dims[ji]][:]))
        nt = length(Array(ds[dims[ti]][:]))
        idx = Any[Colon(), Colon(), Colon()]
        idx[ti] = nt
        x = Array(v[idx...])
        li > ji && (x = permutedims(x))
        return lon, lat, Float64.(x) ./ 100.0
    end
end

"Nearest-neighbor map of model grid onto the obs 2.5 degree grid."
function align(mlon, mlat, x, olon, olat)
    lonmap = [argmin([min(mod(ml - ol, 360.0), mod(ol - ml, 360.0)) for ml in mlon]) for ol in olon]
    latmap = [argmin(abs.(mlat .- ol)) for ol in olat]
    return [x[lonmap[i], latmap[j]] for i in eachindex(olon), j in eachindex(olat)]
end

"Block-average a fine (lon, lat) field onto the obs 2.5 degree grid."
function to_obs_grid(flon, flat, x, olon, olat)
    out = fill(NaN, length(olon), length(olat))
    lonbin = [argmin([min(mod(fl - ol, 360.0), mod(ol - fl, 360.0)) for ol in olon]) for fl in flon]
    latbin = [argmin(abs.(olat .- fl)) for fl in flat]
    num = zeros(length(olon), length(olat))
    den = zeros(length(olon), length(olat))
    for i in eachindex(flon), j in eachindex(flat)
        v = x[i, j]
        isfinite(v) || continue
        num[lonbin[i], latbin[j]] += v
        den[lonbin[i], latbin[j]] += 1
    end
    for k in eachindex(out)
        den[k] > 0 && (out[k] = num[k] / den[k])
    end
    return out
end

"ERA5 Sep 2006-2010 mean of one variable at one pressure level (Pa)."
function era5_sep_mean(varname, plev_pa)
    NCDatasets.NCDataset(ERA_FP) do ds
        lon = Float64.(Array(ds["longitude"][:]))
        lat = Float64.(Array(ds["latitude"][:]))
        pl = Array(ds["pressure_level"][:])
        pk = argmin(abs.(pl .- plev_pa))
        abs(pl[pk] - plev_pa) < 1 || @warn "nearest level $(pl[pk]) for $plev_pa"
        t = Array(ds["time"][:])
        ks = [k for (k, tt) in enumerate(t) if
              Dates.month(Dates.DateTime(tt)) == 9 && Dates.year(Dates.DateTime(tt)) in 2006:2010]
        length(ks) == 5 || error("expected 5 Septembers in ERA5, got $(length(ks))")
        v = ds[varname]  # (time, pressure_level, latitude, longitude)
        acc = zeros(length(lon), length(lat))
        for k in ks
            x = Array(v[:, :, pk, k])  # NCDatasets order: (lon, lat, level, time) reversed
            acc .+= Float64.(replace(x, missing => NaN))
        end
        return lon, lat, acc ./ length(ks)
    end
end

"Ocean mask on the obs grid from MAC lwp coverage (ocean-only product)."
function ocean_mask(olon, olat)
    NCDatasets.NCDataset(MAC_FP) do ds
        v = ds["cloudlwp"]
        lon = Float64.(Array(ds["lon"][:]))
        lat = Float64.(Array(ds["lat"][:]))
        nt = size(v, 3)
        ks = round.(Int, range(1, nt; length = min(12, nt)))
        finite01 = zeros(length(lon), length(lat))
        for k in ks
            x = Float64.(replace(Array(v[:, :, k]), missing => NaN))
            finite01 .+= map(t -> isfinite(t) ? 1.0 : 0.0, x)
        end
        finite01 ./= length(ks)
        frac = to_obs_grid(lon, lat, finite01, olon, olat)
        return [isfinite(frac[k]) && frac[k] >= 0.5 for k in CartesianIndices(frac)]
    end
end

zonal_anom(f) = begin
    g = copy(f)
    for j in axes(f, 2)
        col = filter(isfinite, f[:, j])
        m = isempty(col) ? NaN : Statistics.mean(col)
        g[:, j] .= f[:, j] .- m
    end
    g
end

# ---- assemble ---------------------------------------------------------------
olon, olat, obsf = load_obs()
mask = trues(length(olon), length(olat))
for f in values(obsf); global mask = mask .& isfinite.(f); end
obsmean = [Statistics.mean(obsf[y][k] for y in 2006:2010) for k in CartesianIndices(mask)]
obssd = [Statistics.std([obsf[y][k] for y in 2006:2010]) for k in CartesianIndices(mask)]

pdir = joinpath(ENSEMBLES[2].dir, "member_001", ENSEMBLES[2].simsub)
ml, mt, mclt = load_member(pdir)
model = align(ml, mt, mclt, olon, olat)
r = [mask[k] && isfinite(model[k]) ? model[k] - obsmean[k] : NaN for k in CartesianIndices(mask)]
ra = zonal_anom(r)

elon, elat, t700 = era5_sep_mean("t", 70000.0)
_, _, t1000 = era5_sep_mean("t", 100000.0)
_, _, w500 = era5_sep_mean("w", 50000.0)
theta700 = t700 .* (100000.0 / 70000.0)^0.2854
lts = to_obs_grid(elon, elat, theta700 .- t1000, olon, olat)
om500 = to_obs_grid(elon, elat, w500, olon, olat)
ocean = ocean_mask(olon, olat)
nfin(x) = count(k -> isfinite(x[k]), CartesianIndices(x))
println("diag: ocean cells ", count(ocean), " / ", length(ocean),
    "  lts finite ", nfin(lts), "  om500 finite ", nfin(om500),
    "  r finite ", nfin(r))

w = [cosd(olat[j]) for i in eachindex(olon), j in eachindex(olat)]
gm = Statistics.mean(filter(isfinite, vec(obsmean)))
println("residual: global cos-mean r = ",
    round(sum(w[k] * r[k] for k in CartesianIndices(r) if isfinite(r[k])) /
          sum(w[k] for k in CartesianIndices(r) if isfinite(r[k])); digits = 4),
    "  (obs mean clt $(round(gm; digits=3)); negative = too few clouds)")

# ---- composites -------------------------------------------------------------
function composite(field, pred, sel, label)
    ks = [k for k in CartesianIndices(field) if
          sel[k] && isfinite(field[k]) && isfinite(pred[k])]
    isempty(ks) && return println("$label: no points")
    ps = [pred[k] for k in ks]
    edges = Statistics.quantile(ps, 0:0.1:1)
    tot = sum(w[k] * field[k]^2 for k in ks)
    println("\n$label (10 bins by predictor decile; share of cos-weighted sum r^2, mean r):")
    for b in 1:10
        lo, hi = edges[b], edges[b + 1]
        kb = [k for k in ks if lo <= pred[k] <= (b == 10 ? hi : hi - eps(hi))]
        isempty(kb) && continue
        share = sum(w[k] * field[k]^2 for k in kb) / tot
        meanr = sum(w[k] * field[k] for k in kb) / sum(w[k] for k in kb)
        println("  bin ", lpad(b, 2), "  [", rpad(round(lo; digits = 2), 7), ",",
            rpad(round(hi; digits = 2), 7), "]  share ",
            rpad(round(share; digits = 3), 6), "  mean r ", round(meanr; digits = 3))
    end
end

sel_ocean = [ocean[k] && abs(olat[k[2]]) <= 60 for k in CartesianIndices(r)]
sel_trop = [ocean[k] && abs(olat[k[2]]) <= 30 for k in CartesianIndices(r)]
composite(r, lts, sel_ocean, "r vs LTS (ocean 60S-60N)")
composite(ra, lts, sel_ocean, "r_anom vs LTS (ocean 60S-60N)")
composite(r, om500, sel_trop, "r vs omega500 (ocean 30S-30N)")
composite(ra, om500, sel_trop, "r_anom vs omega500 (ocean 30S-30N)")

# ---- regression -------------------------------------------------------------
function wregress(field, preds, sel, label)
    ks = [k for k in CartesianIndices(field) if
          sel[k] && isfinite(field[k]) && all(isfinite(p[k]) for p in preds)]
    isempty(ks) && return println(label, ": no points")
    y = Float64[field[k] for k in ks]
    ww = Float64[w[k] for k in ks]
    X = hcat(ones(length(ks)), (Float64[p[k] for k in ks] for p in preds)...)
    W = LinearAlgebra.Diagonal(ww)
    beta = (X' * W * X) \ (X' * W * y)
    yhat = X * beta
    ybar = sum(ww .* y) / sum(ww)
    r2 = 1 - sum(ww .* (y .- yhat) .^ 2) / sum(ww .* (y .- ybar) .^ 2)
    println(label, ": R2 = ", round(r2; digits = 3), "  beta = ",
        [round(b; sigdigits = 3) for b in beta])
end
println()
wregress(r, [lts], sel_ocean, "r ~ LTS (ocean)")
wregress(r, [lts, om500], sel_ocean, "r ~ LTS + omega500 (ocean)")
wregress(ra, [lts, om500], sel_ocean, "r_anom ~ LTS + omega500 (ocean)")

# ---- B9 alignment -----------------------------------------------------------
println("\nB9 alignment: cos angle between signed parameter response and residual,")
println("noise-whitened (sigma = interannual obs std), cos(lat)-weighted.")
sig = [max(obssd[k], 1e-3) for k in CartesianIndices(mask)]
function proj(delta, target)
    ks = [k for k in CartesianIndices(target) if
          isfinite(delta[k]) && isfinite(target[k]) && mask[k]]
    num = sum(w[k] * delta[k] * target[k] / sig[k]^2 for k in ks)
    d2 = sum(w[k] * delta[k]^2 / sig[k]^2 for k in ks)
    t2 = sum(w[k] * target[k]^2 / sig[k]^2 for k in ks)
    return num / sqrt(d2 * t2)
end
import TOML
for ens in ENSEMBLES
    nmem = length(filter(d -> startswith(d, "member_"), readdir(ens.dir)))
    memdir(m) = joinpath(ens.dir, "member_" * lpad(m, 3, '0'), ens.simsub)
    _, _, c0 = load_member(memdir(1))
    center = align(ml, mt, c0, olon, olat)
    t1 = Dict(k => v["value"] for (k, v) in TOML.parsefile(joinpath(ens.dir, "member_001", "parameters.toml")))
    println("\n  ", ens.name, " (vs r | vs r_anom):")
    # Sigma points are ordered +p1..+pP then -p1..-pP: parameter k is
    # the pair (member 1+k, member 1+P+k). Verified by the TOML labels
    # of both members agreeing.
    P = (nmem - 1) ÷ 2
    for k in 1:P
        _, _, xa = load_member(memdir(1 + k))
        _, _, xb = load_member(memdir(1 + P + k))
        (isnothing(xa) || isnothing(xb)) && continue
        fa = align(ml, mt, xa, olon, olat)
        fb = align(ml, mt, xb, olon, olat)
        delta = (fa .- fb) ./ 2
        diffof(m) = begin
            tm = Dict(kk => v["value"] for (kk, v) in TOML.parsefile(joinpath(ens.dir, "member_" * lpad(m, 3, '0'), "parameters.toml")))
            [kk for kk in keys(t1) if abs(tm[kk] - t1[kk]) > 1e-9 * max(abs(t1[kk]), 1e-300)]
        end
        da, db = diffof(1 + k), diffof(1 + P + k)
        pname = (length(da) == 1 && da == db) ? da[1] : "pair$(k) MISMATCH $(da)/$(db)"
        println("    ", rpad(pname, 48),
            rpad(round(proj(delta, r); digits = 3), 8), " | ",
            round(proj(delta, zonal_anom(r)); digits = 3))
    end
end
println()
