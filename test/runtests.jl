using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using Pkg.Artifacts
using PkgUtility
using Test


ENV["JULIA_LOG_LEVEL"] = "WARN"


@testset verbose = true "CliMA Land v0.1" begin
    include("test_CanopyLayers.jl");
    include("test_Photosynthesis.jl");
    include("test_PlantHydraulics.jl");
    include("test_StomataModels.jl");
    include("test_SPAC.jl");
end;
