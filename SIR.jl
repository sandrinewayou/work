"""
    Implementation of the SIR model to compare the speed of the implementation
    of Gillespies algorithms with others.

    The population state is represented as vector with three entries.
    First the number of susceptible Individuals, second the infected and
    third the recovered individuals.
"""

function infection!(x)
    x[1] += -1
    x[2] += 1
    nothing
end

function recovery!(x)
    x[2] += -1
    x[3] += 1
    nothing
end

function execute!(i,x,par)
    if i == 1
        infection!(x)
    elseif i == 2
        recovery!(x)
    else
        error("Unknown event number i = $i")
    end
    nothing
end

function rates!(rates,x,par)
    #rate of infection
    rates[1] = par.β * x[1] * x[2]
    #rate of recovery
    rates[2] = par.γ * x[2]
    nothing
end

#------

include("../MainFunctions.jl")
import .Gillespie

using Plots

par = (
    β = 5/1000000,
    γ = 0.005
    )

x0 = [9999,1,0]
t = 0:1000

hist = zeros(Int,(length(t),3))

result = Gillespie.run_gillespie!(
        t,x0,par,
        execute!,rates!,
        Vector{Float64}(undef,2),hist
        )

#plot simulation
plot(hist,label=["S" "I" "R"])
