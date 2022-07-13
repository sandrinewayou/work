include("MainFunctions.jl")
include("modifiedmodel.jl")

import .Gillespie
import .modified

using Plots


#contact parameter
t_end = 300
a     = 0.3
r_0   = 3
ind   = [666,333]
infect_period = [5,10]
p_d1 = 0.2
p_d2 = 0.5
p_d3 = 0.7
n = 4
λ_n = 2

t_start = 20
t_1 = 50
α_1 = 1/65
α_2 = 1/35
α_3 = 1/25
#parameters

par = (
	t_start = 20,
	t_1 = 50,
	α_1 = 1/65,
	α_2 = 1/35,
	α_3 = 1/25,
      p = 0.5,
      η = 0.5,
      σ_1 = 1/10,
      σ_2 = 2/10 ,
      f_dead1 = 0.03,
      f_dead2 = 0.02,
      β = 1/3,
      μ =1/200,
      contact_par = [t_end,a,r_0,infect_period,ind],
      gen_par = [p_d1,p_d2,p_d3],
      ext = [n,λ_n],
	va = [t_start,t_1,α_1,α_2],
	vb =[t_start,t_1,α_2,α_3]
)


#initial condition
x0 = [0,0,0,665,1,0,333,0,0,0,0,0,0]
t = 0:50

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
plot([hist[:,1],hist[:,10]],label=["S" "pi" "I" "Iv" "D" "R"])
#plot(hist[i][1] for i in t,label=["S"])
