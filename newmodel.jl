#try something
function vaccination!(x)
    x[1] +=-1
    x[2] +=1
    nothing
end

function immunity!(x)
    x[2] +=-1
    x[3] +=1
    nothing
end

function infection!(x)
    x[3] +=-1
    x[4] +=1
    nothing
end

function recovery!(x)
    x[4] +=-1
    x[5] +=1
    nothing
end

function susceptibility!(x)
    x[5] +=-1
    x[3] +=1
    nothing
end

function execute!(i,x,par)
    if i == 1
        vaccination!(x)
    elseif i == 2
        immunity!(x)
    elseif i == 3
        infection!(x)
    elseif i == 4
        recovery!(x)
    elseif i == 5
        susceptibility!(x)
    else
        error("Unknown event number i = $i")
    end
    nothing
end

function rates!(rates,x,par,t)
    #rate of vaccination
    rates[1] = par.a * x[1]
    #rate of immunity
    rates[2] = par.b * x[2]
    #rate of infection
    rates[3] = par.β * x[3] * x[4]
    #rate of recovery
    rates[4] = par.γ * x[4]
    #rate of suscepbility
    rates[4] = par.α * x[5]
    nothing
end

function basic_reproduction(t,t_end,a,r_0)
    return r_0(1+a*cos(2*pi*((t-t_end)/365)))
end

function contact_rate(t,t_end,a,r_0,infect_period,ind)
    return (basic_reproduction(t,t_end,a,r_0)/infect_period[1]*ind[1],
            basic_reproduction(t,t_end,a,r_0)/infect_period[1]*ind[2],
            basic_reproduction(t,t_end,a,r_0)/infect_period[2]*ind[1],
            basic_reproduction(t,t_end,a,r_0)/infect_period[2]*ind[2])
end





include("MainFunctions.jl")
import .Gillespie

using Plots

#contact = [0.2,0.04,0.1,0.15]
t_end=500
a=0,5
r_0=4
infect_period=[5,10]
ind=[6666,3333]

par = (
    a = 1/15,
    b=1/50,
    β = sum(contact_rate(t,t_end,a,r_0,infect_period,ind)),
    γ = 0.005,
    α = 1/270,
    #contact = contact
    )

x0 = [sum(ind),0,0,1,0]
t = 0:500

hist = zeros(Int,(length(t),5))

Gillespie.run_gillespie!(
        t,x0,par,
        execute!,rates!,
        Vector{Float64}(undef,4),hist
        )

#plot simulation
plot(hist,label=["S" "V" "PI" "I" "R"])
