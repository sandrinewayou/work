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
      α_1 =
      α_2 =
      p = 0.5
      η = 0.5
      σ_1 =
      σ_2 =
      f_dead1 = 0.03
      f_dead2 = 0.02
      ϕ =
      β =
      μ =
      contact_par = [t_end,a,r_0,infect_period,ind]
)


#initial condition
x0 = [sum(ind),0,1,0,0,0]
t = 0:200


num_events = 19
num_types = 11

hist = zeros(Int,(length(t),6))

Gillespie.run_gillespie!(
        t,x0,par,
       modified.execute!,modified.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "pi" "I" "Iv" "D" "R"])
