#  Modify model

module modified

 #Population moving

 function  plus_minus_one(x,minus_index,plus_index)
     x[minus_index] -= 1
     x[plus_index] += 1
     nothing
 end


 #events

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
        vaccination!(x,4,1)
    elseif i == 2
        vaccination!(x,7,10)
     elseif i == 3
         infection!(x,1,2)
    elseif i == 4
        infection!(x,4,5)
    elseif i == 5
        infection!(x,7,8)
    elseif i == 6
        infection!(x,10,11)
    elseif i == 7
        recovery!(x,2,3)
    elseif i == 8
        recovery!(x,5,6)
    elseif i == 9
        recovery!(x,8,9)
    elseif i == 10
        recovery!(x,11,12)
    elseif i == 11
        dead!(x,2,13)
    elseif i == 12
        dead!(x,5,13)
    elseif i == 13
        dead!(x,8,13)
    elseif i == 14
        dead!(x,11,13)
    elseif i == 15
        suscepbility!(x,3,1)
    elseif i == 16
        suscepbility!(x,6,4)
    elseif i == 17
        suscepbility!(x,9,7)
    elseif i == 18
        suscepbility!(x,12,10)
    elseif i == 19
        lost_imunity!(x,1,4)
    elseif i == 20
        lost_imunity!(x,7,10)
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
function contact_rate!(x,t,t_end,a,r_0,infect_period)
	ind = [x[1]+x[4],x[7]+x[10]]
    r0_t = basic_reproduction(t,t_end,a,r_0)
    return sum(r0_t/(i*j) for i in infect_period for j in ind)
end

#general contact reduction for population b

function contact_reduction(p_d2,p_d3,t,t_c1, t_c2)
    if t< t_c1
        return  p_d3
    elseif t_c1 ≤ t ≤ t_c2
        return  p_d2
    else
       return  0
    end
end
	#external force of infection

	function ext_force(n,λ_n)
		return (ℯ^(-λ_n)*λ_n^n)/factorial(n)
	end


# vaccination rate
function vaccination_rate_a(t,par)
	if par.t_start ≤ t < par.t_1
		return  par.α_1
	elseif t ≥ par.t_1
		return par.α_2
	else
		return 0
	end
end


function vaccination_rate_b(t,par)
	if par.t_start ≤ t < par.t_1
		return par.α_2
	elseif  t ≥ par.t_1
		return par.α_3
	else
		return 0
	end
end

# events rates

function rates!(rates,x,par,t)
    λ = contact_rate!(x,t,par.contact_par...)
	p_n = ext_force(par.ext...)
	α_a = vaccination_rate_a(t,par)
	α_b = vaccination_rate_b(t,par)
	k = (1-contact_reduction(t,par.gen_par...)) * λ * ( x[5] + x[8] + par.η * (x[2] + x[11])+p_n )
    #vaccination rate
    rates[1] = par.α_1 * x[4]
    rates[2] = par.α_2 * x[7]

    #infection rate
    rates[3] = par.p * (1-contact_reduction(t,par.gen_par...)) * λ * ( x[5] + x[8] + par.η * (x[2] + x[11])+p_n ) * (x[1])
    rates[4] = λ * (1-contact_reduction(t,par.gen_par...)) * ( x[5] + x[8] + par.η * (x[2] + x[11])+p_n  ) * (x[4])
    rates[5] = λ * (1-contact_reduction(t,par.gen_par...)) * ( x[5] + x[8] + par.η * (x[2] + x[11]) +p_n ) * (x[7])
    rates[6] = par.p * (1-contact_reduction(t,par.gen_par...)) * λ * ( x[5] + x[8] + par.η * (x[2] + x[11])+p_n  ) * (x[10])

    #recovery rate
    rates[7] = par.σ_1 * (1-par.f_dead2) * x[2]
    rates[8] = par.σ_1 * (1-par.f_dead1) * x[5]

    rates[9]  = par.σ_2 * (1-par.f_dead1) * x[8]
    rates[10] = par.σ_2 * (1-par.f_dead2) * x[11]

    #death rate young
    rates[11] = par.σ_1 * par.f_dead2 * x[2]
    rates[12] = par.σ_1 * par.f_dead1 * x[5]
    #death rate old
    rates[13] = par.σ_2 * par.f_dead1 * x[8]
    rates[14] = par.σ_2 * par.f_dead2 * x[11]

    #suscepbility rate
    rates[15] = par.β * x[3]
    rates[16] = par.β * x[6]
	rates[17] = par.β * x[9]
    rates[18] = par.β * x[12]

    #lost immunity rate
    rates[19] = par.μ * x[1]
    rates[20] = par.μ * x[10]

    nothing
end


# separation graph








end
