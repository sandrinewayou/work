include("MainFunctions.jl")
include("modifiedmodel.jl")

import .Gillespie
import .modified

using Plots


#contact parameter
t_end = 300
a     = 0.5
r_0   = 3.2
ind   = [666,333]
infect_period = [5,10]
p_d1 = 0.2
p_d2 = 0.5
p_d3 = 0.7

#parameters

par = (
      α_1 = 1/50 + 1/15,
      α_2 = (1/50 + 1/15)*2,
      p = 0.5,
      η = 0.5,
      σ_1 = 1/5,
      σ_2 = 1/5 * 2,
      f_dead1 = 0.03,
      f_dead2 = 0.02,
      β = 1/200,
      μ =1/500,
      contact_par = [t_end,a,r_0,infect_period,ind],
      gen_par = [p_d1,p_d2,p_d3]
)


#initial condition
x0 = [0,0,0,665,1,0,333,0,0,0,0,0,0]
t = 0:400

# 1 Sva | 2 Iva | 3 Rva | 4 Sa | 5 Ia | 6Ra
#7 Sb | 8 Ib| 9 Rb | 10 Svb | 11 Ivb | 12 Rvb|13 D
num_events = 20
num_types = 13

hist = zeros(Int,(length(t),13))

Gillespie.run_gillespie!(
        t,x0,par,
       modified.execute!,modified.rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "pi" "I" "Iv" "D" "R"])
#plot(hist[i][1] for i in t,label=["S"])
