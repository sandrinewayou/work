include("MainFunction.jl")
include("segond_trial,jl")

import .Gillespie
import .NewTrial

using Plots






num_types = 5
num_events = 9
hist = zeros(Int,(length(t),5))

Gillespie.run_gillespie!(
        t,x0,par,
        SandrineSIR.execute!,SandrineSIR.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "V" "PI" "I" "R"])
