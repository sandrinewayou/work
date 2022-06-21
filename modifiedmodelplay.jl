include("MainFunctions.jl")
include("modifiedmodel.jl")

import .Gillespie
import .modified

using Plots


#contact parameter
t_end = 300
a     = 0.5
r_0   = 3.2
ind   = [6666,3333]
infect_period = [5,10]

#parameters

par = (
      α_1 = 1/50 + 1/15,
      α_2 = (1/50 + 1/15)*2,
      p = 0.5,
      η = 0.5,
      σ_1 = 1/5,
      σ_2 = 1/5 * 2
      f_dead1 = 0.03
      f_dead2 = 0.02
      ϕ = 1/200,
      β = 1/500,
      μ = 1/1000
      contact_par = [t_end,a,r_0,infect_period,ind]
)


#initial condition
x0 = [sum(ind),0,0,1,0,0,0,0,0,0,0]
t = 0:200


num_events = 19
num_types = 11

hist = zeros(Int,(length(t),11))

Gillespie.run_gillespie!(
        t,x0,par,
       modified.execute!,modified.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "pi" "I" "Iv" "D" "R"])
