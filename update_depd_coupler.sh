# Update checked-in Manifests for Julia v1.10
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.10.9 --project=docs/ -e 'using Pkg; Pkg.update()'
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.10.9 --project=experiments/ClimaCore/ -e 'using Pkg; Pkg.update()'
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.10.9 --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.update()'

# Update checked-in Manifests for Julia v1.11
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.11.8 --project=docs/ -e 'using Pkg; Pkg.update()'
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.11.8 --project=experiments/ClimaCore/ -e 'using Pkg; Pkg.update()'
JULIA_PKG_PRECOMPILE_AUTO=0 julia +1.11.8 --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.update()'
