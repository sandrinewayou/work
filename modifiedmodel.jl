#  Modify model

module modified

 #Population moving

 function  plus_minus_one(x,minus_index,plus_index)
     x[minus_index] -= 1
     x[plus_index] += 1
     nothing
 end


 #events

grow!(x) = plus_minus_one(x,2,6)

#  (m,j)_a = (2,1), (m,j)_b = (6,7)
vaccination!(x,m,j) =  plus_minus_one(x,m,j)

# (m,j)_a = [(1,3),(2,4)], (m,j)_b = [(6,8),(7,9)]
infection!(x,m,j) =  plus_minus_one(x,m,j)

 # (m,j)_a = [(3,5),(4,5)], (m,j)_b = [(8,10),(9,10)]
recovery!(x,m,j) = plus_minus_one(x,m,j)

# (m,j)_a = [(3,11),(4,11)], (m,j)_b = [(8,11),(9,11)]
dead!(x,m,j) = plus_minus_one(x,m,j)

#  (m,j)_a = (5,1), (m,j)_b = (10,7)]
suscepbility!(x,m,j) = plus_minus_one(x,m,j)

#  (m,j)_a = (1,2), (m,j)_b = (7,6)
lost_imunity!(x,m,j) = plus_minus_one(x,m,j)


# events execution

function execute!(i,x,par)
    if i == 1
        vaccination!(x,2,1)
    elseif i == 2
        vaccination!(x,6,7)
     elseif i == 3
         infection!(x,1,3)
    elseif i == 4
        infection!(x,2,4)
    elseif i == 5
        infection!(x,6,8)
    elseif i == 6
        infection!(x,7,9)
    elseif i == 7
        recovery!(x,3,5)
    elseif i == 8
        recovery!(x,4,5)
    elseif i == 9
        recovery!(x,8,10)
    elseif i == 10
        recovery!(x,9,10)
    elseif i == 11
        dead!(x,3,11)
    elseif i == 12
        dead!(x,4,11)
    elseif i == 13
        dead!(x,8,11)
    elseif i == 14
        dead!(x,9,11)
    elseif i == 15
        suscepbility!(x,5,1)
    elseif i == 16
        suscepbility!(x,10,7)
    elseif i == 17
        lost_imunity!(x,1,2)
    elseif i == 18
        lost_imunity!(x,7,6)
    elseif i == 19
        grow!(x)
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
    return sum(r0_t/(i*j) for i in infect_period for j in ind)
end

# events rates

function rates!(rates,x,par,t)
    λ = contact_rate(t,par.contact_par...)
    #vaccination rate
    rates[1] = par.α_1 * x[2]
    rates[2] = par.α_2 * x[6]

    #infection rate
    rates[3] = par.p * λ * ( x[4] + x[8] + par.η * (x[3] + x[9]) ) * x[1]
    rates[4] = λ * ( x[4] + x[8] + par.η * (x[3] + x[9]) ) * x[2]
    rates[5] = λ * ( x[4] + x[8] + par.η * (x[3] + x[9]) ) * x[6]
    rates[6] = par.p * λ * ( x[4] + x[8] + par.η * (x[3] + x[9]) ) * x[7]

    #recovery rate
    rates[7] = par.σ_1 * (1-par.f_dead2) * x[3]
    rates[8] = par.σ_1 * (1-par.f_dead1) * x[4]

    rates[9]  = par.σ_2 * (1-par.f_dead1) * x[8]
    rates[10] = par.σ_2 * (1-par.f_dead2) * x[9]

    #death rate young
    rates[11] = par.σ_1 * par.f_dead2 * x[3]
    rates[12] = par.σ_1 * par.f_dead1 * x[4]
    #death rate old
    rates[13] = par.σ_2 * par.f_dead1 * x[8]
    rates[14] = par.σ_2 * par.f_dead2 * x[9]

    #suscepbility rate
    rates[15] = par.ϕ * x[5]
    rates[16] = par.ϕ * x[10]

    #lost immunity rate
    rates[17] = par.β * x[1]
    rates[18] = par.β * x[7]

    #grow rate
    rates[19] = par.μ * x[2]

    nothing
end









end
