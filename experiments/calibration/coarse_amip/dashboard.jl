using Bonito
using WGLMakie
import GeoMakie as GM
using ClimaAnalysis

# Load data via ClimaAnalysis
simdir = SimDir("/Users/nefrathenrici/Downloads/clima_atmos")
rsut = get(simdir, "rsut")
pr = get(simdir, "pr")
lon = pr.dims["lon"]
lat = pr.dims["lat"]

function rsut_land(slider, menu, fig, ax)
    s = slider.value
    m = menu.value
    var = @lift(Dict("rsut" => rsut, "pr" => pr)[$m])
    varslice = @lift(slice($var, time = $var.dims["time"][$s])) # Get rsut
    slice_land = @lift(apply_oceanmask($varslice)) # Apply oceanmask
    slice_land_data = @lift($slice_land.data) # Get data
    p = heatmap!(ax, lon, lat, slice_land_data) # Plot
    return fig
end

app = App() do
    fig = Figure()
    ax = GM.GeoAxis(fig[1, 1]; dest = "+proj=wintri")
    lines!(ax, GM.coastlines())
    menu = Dropdown(["rsut", "pr"])
    slider = StylableSlider(1:12)
    landmap = rsut_land(slider, menu, fig, ax)
    return DOM.div(slider, menu, landmap)
end
