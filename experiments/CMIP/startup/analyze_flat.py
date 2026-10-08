#!/usr/bin/env python3
"""Attribute a Julia `Profile.print(format=:flat, C=true)` dump to phases of work.

The Count column of a flat dump is the number of samples in which that frame
appeared anywhere on the stack, so picking one *anchor* frame per activity gives
that activity's share directly. Anchors are chosen to sit at the top of their
subtree, so nested work (e.g. the Julia inference that GPUCompiler itself drives)
is counted inside the enclosing anchor rather than separately.
"""
import re
import sys
from collections import defaultdict

# anchor label -> (regex on function name, regex on file path); either may be None
ANCHORS = [
    ("GPU compile: total",      r"cufunction|cached_compilation|GPUCompiler\.compile", None),
    ("GPU compile: ptxas/link", r"run_and_collect|\blink\b|ptxas|nvlink",              None),
    ("GPU compile: LLVM/IR",    r"emit_llvm|emit_asm|irgen|optimize_module",           None),
    ("Julia inference",         r"jl_type_infer|typeinf_ext_toplevel|typeinf_edge",    None),
    ("Julia codegen",           r"jl_create_native|_jl_compile_codeinst|jl_generate_fptr|jl_compile_method_internal|emit_function", None),
    ("GC",                      r"_jl_gc_collect|jl_gc_collect|gc_mark_loop",          None),
    ("NetCDF / HDF5 I/O",       r"^nc_|^H5",                      r"NCDatasets|HDF5|NetCDF_jll"),
    ("Regridding (polygons)",   None,                             r"ConservativeRegridding"),
    ("Oceananigans setup",      None,                             r"/Oceananigans/"),
    ("ClimaAtmos setup",        None,                             r"/ClimaAtmos/"),
    ("ClimaLand setup",         None,                             r"/ClimaLand/"),
    ("ClimaCore",               None,                             r"/ClimaCore/"),
]

LINE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(.*)$")


def parse(path):
    rows, total = [], None
    for line in open(path):
        if line.startswith("Total snapshots:"):
            m = re.search(r"Total snapshots:\s*(\d+)", line)
            if m:
                total = int(m.group(1))
            continue
        m = LINE.match(line)
        if not m or line.lstrip().startswith("====") or "Count" in line:
            continue
        count, _overhead, fil, _ln, func = m.groups()
        rows.append((int(count), fil, func))
    return rows, total


def main(path, delay=0.01):
    rows, total = parse(path)
    if not rows:
        print(f"{path}: no parseable rows")
        return
    if total is None:
        total = max(c for c, _, _ in rows)
    print(f"\n### {path}")
    print(f"total samples: {total}   (~{total * delay:.0f} s at {delay}s sampling)\n")
    print(f"{'activity':28s}{'samples':>10s}{'% of phase':>12s}{'~seconds':>11s}")
    print("-" * 61)
    for label, fpat, ppat in ANCHORS:
        best = 0
        for count, fil, func in rows:
            if fpat and re.search(fpat, func):
                best = max(best, count)
            elif ppat and re.search(ppat, fil):
                best = max(best, count)
        if best:
            print(f"{label:28s}{best:>10d}{100*best/total:>11.1f}%{best*delay:>10.0f}s")
    print()


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
