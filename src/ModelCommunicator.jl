"""
    ModelComms

This module contains functions to regrid information between spaces.
Many of the functions used in this module call TempestRemap functions
via ClimaCoreTempestRemap wrappers.
"""
module ModelComms

using ..Utilities
using ..TimeManager
using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using NCDatasets
using ClimaCoreTempestRemap
using Dates
using JLD2

export AbstractComponentModel

abstract type AbstractComponentModel end





end # Module
