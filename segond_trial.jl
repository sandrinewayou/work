#new trial
module NewTrial

#population moving
function  plus_minus_one(x,plus_index, minus_index)
    x[plus_index] += 1
    x[minus_index] -= 1
    nothing
end

    #events
vaccination!(x) =  plus_minus_one(x,2, 1)
infection!(x,1) =  plus_minus_one(x,3, 1)
infection!(x,2) =  plus_minus_one(x,4, 2)
death!(x,3)     =  plus_minus_one(x,5, 3)
death!(x,4)     =  plus_minus_one(x,5, 4)
recovery!(x,3)  =  plus_minus_one(x,6, 3)
recovery!(x.4)  =  plus_minus_one(x,6, 4)
immunity!(x)    =  plus_minus_one(x,2, 6)
inimmunity!(x)  =  plus_minus_one(x,1, 2)

#execution function
function execute!(i,x,par,)
    if i == 1
        vaccination!(x)
    elseif i == 2
        infection!(x,1)
    elseif i == 3
        infection!(x,2)
    elseif i == 4
        death!(x,3)
    elseif i == 5
        death!(x,4)
    elseif i == 6
        recovery!(x,3)
    elseif i == 7
        recovery!(x,4)
    elseif i == 8
        immunity!(x)
    elseif i == 9
        innimunity!(x)
    else
        error("unknown event number i = $i")
    end
    nothing
end

#basic reproduction number
function basic_reproduction(t,t_end,a,r_0)
    return r_0*(1+a*cos(2π*((t-t_end)/365)))
end

#contact rate
function contact_rate(t,t_end,a,r_0,infect_period,ind)
    r0_t = basic_reproduction(t,t_end,a,r_0)
    return sum( br/(i*j) for i in infect_period for j in ind)
end

#rate function
function rates!(rates,x,par,t)
    λ = contact_rate(t,par.contact_par...)
    #vaccination rate
    rate[1] = par.α x[1]
    #infection rate
    rate[2] = λ * (x[3] + x[4]) * x[1]
    rate[3] = p * λ * (x[3] + x[4]) * x[2]
    #death rate
    rate[4] = par.γ_1 * x[3]
    rate[5] = par.γ_2 * x[4]
    #recovery rate
    rate[6] = par.δ_1 * x[3]
    rate[7] = par.δ_2 * x[4]
    #immunity rate
    rate[8] = par.ϕ * x[6]
    #inimmunity rate
    rate[9] = par.β * x[2]
    nothing
end


end #end module NewTrial
