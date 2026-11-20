# Helpers `include`d by each warmup package (rather than a shared dependency, which would
# change the package graph described in README.md).

"""
    warmup_dir()

Return the path to `experiments/warmup`.
"""
warmup_dir() = @__DIR__

"""
    repo_dir()

Return the repository root.
"""
repo_dir() = normpath(joinpath(@__DIR__, "..", ".."))

"""
    config_path(name)

Return the path to the config `name` in `experiments/warmup/configs`.
"""
config_path(name) = joinpath(warmup_dir(), "configs", name)

"""
    build_sim(outdir)

Build the atmosphere simulation from `configs/atmos_warmup.yml`, writing output to `outdir`.
"""
function build_sim(outdir)
    cfg = YAML.load_file(config_path("atmos_warmup.yml"))
    cfg["output_dir"] = outdir
    # Parameter files are given relative to the repository root.
    haskey(cfg, "toml") &&
        (cfg["toml"] = [isabspath(f) ? f : joinpath(repo_dir(), f) for f in cfg["toml"]])
    return CA.get_simulation(CA.AtmosConfig(cfg; job_id = "atmos_warmup"))
end

"""
    small_gpu_atmos_config(path, outdir)

Return `path`, or, when `WARMUP_H_ELEM` is set, a copy in `outdir` with `h_elem`
overridden, so the coupled workloads fit on a small GPU.
"""
function small_gpu_atmos_config(path, outdir)
    h_elem = get(ENV, "WARMUP_H_ELEM", "")
    isempty(h_elem) && return path
    cfg = YAML.load_file(path)
    cfg["h_elem"] = parse(Int, h_elem)
    out = joinpath(outdir, "atmos_warmup_h$(h_elem).yml")
    YAML.write_file(out, cfg)
    return out
end

"""
    warmup_config(outdir)

Load `configs/cmip_warmup.yml` as a coupler config dict, with relative paths resolved and
output written to `outdir`.
"""
function warmup_config(outdir)
    cfg = YAML.load_file(config_path("cmip_warmup.yml"))
    cfg["coupler_output_dir"] = outdir
    cfg["job_id"] = "warmup"
    for key in ("atmos_config_file",)
        v = get(cfg, key, nothing)
        v isa AbstractString && !isabspath(v) && (cfg[key] = joinpath(repo_dir(), v))
    end
    cfg["atmos_config_file"] = small_gpu_atmos_config(cfg["atmos_config_file"], outdir)
    tmp = joinpath(outdir, "cmip_warmup_resolved.yml")
    YAML.write_file(tmp, cfg)
    return Input.get_coupler_config_dict(tmp)
end

"""
    write_report(path, status; h_elem_override = "")

Write the workload `status`, any `h_elem_override` and the device to `path`. Never throws.
"""
function write_report(path, status; h_elem_override = "")
    try
        open(path, "w") do io
            println(io, "status = ", status)
            isempty(h_elem_override) || println(
                io,
                "h_elem override = ",
                h_elem_override,
                " (type-neutral; see README)",
            )
            println(io, "device = ", get(ENV, "CLIMACOMMS_DEVICE", "(unset)"))
        end
    catch
    end
end
