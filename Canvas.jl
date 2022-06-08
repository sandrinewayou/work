include("MainFunctions.jl")
include("SandrineSIR.jl")

import .Gillespie
import .SandrineSIR

using Plots

#---

#contact = [0.2,0.04,0.1,0.15]
t_end=500
a=0.5
r_0=0.1
infect_period=[5,10]
ind=[6666,3333]

par = (
    α = 1/15,
    β = 1/50,
    #β = contact_rate(0,t_end,a,r_0,infect_period,ind),
    γ = 0.05,
    δ = 1/270,
    contact_par = [t_end,a,r_0,infect_period,ind]
    #contact = contact
    )

x0 = [sum(ind),0,0,1,0]
t = 0:1500

# 1 S | 2 V | 3 PI | 4 I | 5 R
num_types = 5
num_events = 6
hist = zeros(Int,(length(t),5))

Gillespie.run_gillespie!(
        t,x0,par,
        SandrineSIR.execute!,SandrineSIR.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "V" "PI" "I" "R"])
