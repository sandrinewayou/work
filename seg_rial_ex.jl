include("MainFunctions.jl")
include("segond_trial.jl")

import .Gillespie
import .NewTrial

using Plots

f_dead1 = 0.03
f_dead2 = 0.02
p       = 0.5


#contact parameter
t_end = 300
a     = 0.5
r_0   = 3.2
ind   = [6666,3333]
infect_period = [5,10]
# Paramters
par = (
     α   = 1/50 + 1/15,
     σ       = 1/5,
        # γ_1 = σ * f_dead1,
        # γ_2 = σ * f_dead2,
        # δ_1 = σ * (1-f_dead1),
        # δ_2 = σ * (1-f_dead2),
    ϕ   = 1/200,
    β   = 1/500,
    contact_par = [t_end,a,r_0,infect_period,ind]
    )

#FRACTION




#initial condition
x0 = [sum(ind),0,1,0,0,0]
t = 0:1500


num_events = 9
num_types = 6

hist = zeros(Int,(length(t),6))

Gillespie.run_gillespie!(
        t,x0,par,
        NewTrial.execute!,NewTrial.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "pi" "I" "Iv" "D" "R"])
