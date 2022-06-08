#rates functions

function plusminusone(x,plus_index,minus_index)
    x[plus_index] += 1
    x[minus_index] -= 1
    nothing
end
vaccination!(x) = plusminusone(x,2,1)
immunity!(x) = plusminusone(x,3,2)
infection!(x,j) = plusminusone(x,4,j)
recovery!(x) = plusminusone(x,5,4)
susceptibility!(x) = plusminusone(x,3,5)

function execute!(i,x,par)
    if i == 1
        vaccination!(x)
    elseif i == 2
        immunity!(x)
    elseif i == 3
        infection!(x,3)
    elseif i == 4
        recovery!(x)
    elseif i == 5
        susceptibility!(x)
    elseif i == 6
        infection!(x,1)
    else
        error("Unknown event number i = $i")
    end
    nothing
end

function rates!(rates,x,par,t)
    λ = contact_rate(t,par.contact_par...)
    #rate of vaccination
    rates[1] = par.α * x[1]
    #rate of immunity
    rates[2] = par.β * x[2]
    #rate of infection
    rates[3] = λ * x[3] * x[4]
    #rate of recovery
    rates[4] = par.γ * x[4]
    #rate of suscepbility
    rates[5] = par.δ * x[5]
    #rate of infection of susceptibil individuals
    rates[6] = λ * x[1] * x[4]
    nothing
end

function basic_reproduction(t,t_end,a,r_0)
    return r_0*(1+a*cos(2π*((t-t_end)/365)))
end

function contact_rate(t,t_end,a,r_0,infect_period,ind)
    br = basic_reproduction(t,t_end,a,r_0)
    return sum( br/(i*j) for i in infect_period for j in ind)
end

#----

include("MainFunctions.jl")
import .Gillespie

using Plots

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
        execute!,rates!,
        Vector{Float64}(undef,num_events),hist
        )

#plot simulation
plot(hist,label=["S" "V" "PI" "I" "R"])
