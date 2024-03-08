using Plots

num_gpus = [1, 2, 4]

# h_elem: 42, 60, 84
ws_sypd = [0.8, 1.2, 2.1]
ws_res = ["77 km" "55 km" "39 km"]

# h_elem: 30
ss_30_sypd = [1.0, 1.7, 3.2]
ss_30_res = 110

# h_elem: 60
ss_60_sypd = [0.7, 1.5, 2.8]
ss_60_res = 55


scatter(num_gpus,
    ws_sypd,
    title = "AMIP GPU scaling",
    xlabel = "Number of GPUs",
    ylabel = "Simulated years per day (SYPD)",
    labels = "Weak Scaling [various res.]",
    right_margin = 20Plots.mm,
)
annotate!(num_gpus[1] + 0.38, ws_sypd[1], text(ws_res[1], 1, :right, 10))
annotate!(num_gpus[2] + 0.38, ws_sypd[2], text(ws_res[2], 1, :right, 10))
annotate!(num_gpus[3] + 0.38, ws_sypd[3], text(ws_res[3], 1, :right, 10))

scatter!(num_gpus,
    ss_30_sypd,
    label = "Strong Scaling [$ss_30_res km]",
    color = 3,
)
# annotate!(num_gpus[1] + 0.52, ss_sypd[1], text(ss_res[1], 1, :right, 10))
# annotate!(num_gpus[2] + 0.52, ss_sypd[2], text(ss_res[2], 1, :right, 10))
# annotate!(num_gpus[3] + 0.52, ss_sypd[3], text(ss_res[3], 1, :right, 10))

scatter!(num_gpus,
    ss_60_sypd,
    label = "Strong Scaling [$ss_60_res km]",
    color = 4,
)
