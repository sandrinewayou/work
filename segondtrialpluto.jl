### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 42e093ec-dbda-4474-bc31-8ad3bdfda9e7
begin
	using Pkg; Pkg.add("Distributions")
	using Plots
	using PlutoUI
end

# ╔═╡ c5576946-69fd-4ac9-986f-0afb74b70137
md"""
**Enter full path to locally saved files**
"""

# ╔═╡ e4edc158-42b3-4bf1-8420-14b2bd19cd05
#Local Path to MainFunctions.jl"
main_functions_path = "/home/larocca/github/work/MainFunctions.jl";

# ╔═╡ 774add99-f8b9-4ae5-8a08-cc6cc71f64a7
#Local Path to Modelname.jl"
model_path = "/home/larocca/github/work/SandrineSIR.jl";

# ╔═╡ 59e73279-85bf-4c77-a207-cfd1106c1451
md"""---
###### Fixed parameters
"""

# ╔═╡ 85412cf0-362c-4fa5-8c37-eb261db37876
md""" Time $(@bind t_end Slider(500:100:10_000,show_value=true,default=1500))
Population Size $(@bind s0 Slider(1000:1000:100_000,show_value=true,default=10_000))
"""

# ╔═╡ 2684693e-e731-40e0-a22d-35c62110009c
initial_population = [s0-1,0,1,0,0,0];

# ╔═╡ 3b6e0504-40a1-409e-a4ef-b39be2bef971
md"""---
**Model Parameter**
"""

# ╔═╡ c02d5705-f136-433a-98d4-eab309926694
md""" Basic infection rate: $\quad r_0$ $(@bind r_0  Slider(0.001:0.001:5.0,show_value=true,default=0.1)) ;
"""

# ╔═╡ 5fa15756-0c38-4f9f-b36d-2f8d3a5f9e6a
md"""Strength of Seasonal Effect $\quad a$ $(@bind a  Slider(0.001:0.001:1.0,show_value=true,default=0.1))
"""

# ╔═╡ 9d5c0f86-1e41-47f4-9dea-503d3bfbd3d5
md"""Rate of vaccination $\quad \alpha$ $(@bind al  Slider(0.001:0.001:1.0,show_value=true,default=0.087))
"""

# ╔═╡ 397ac091-4c45-4fa1-9b9b-24a966661a13
md"""Rate of suscepbility $\quad \phi$ $(@bind b  Slider(0.001:0.001:1.0,show_value=true,default=0.005))
"""

# ╔═╡ 75002645-eeb0-43dd-bf53-25cc9a6b089d
md"""Rate of recovery or dead $\quad \delta$ $(@bind g  Slider(0.001:0.001:1.0,show_value=true,default=0.2))

Rate of pi infectin $\quad p$ $(@bind Pi  Slider(0.001:0.001:1.0,show_value=true,default=0.5))

Rate of  dead1 $\quad f_dead1$ $(@bind f1  Slider(0.001:0.001:1.0,show_value=true,default=0.03))

Rate of dead 2 $\quad f_dead2$ $(@bind f2  Slider(0.001:0.001:1.0,show_value=true,default=0.02))
"""

# ╔═╡ 193be177-08a1-42ed-98df-8d58031cd0ae
md"""Rate of loss of immunity $\quad \beta$ $(@bind d  Slider(0.001:0.001:1.0,show_value=true,default=0.002))
"""

# ╔═╡ 5ba2c008-5020-4fcd-82be-ec8089bcc6fe
md"""Proportion of young to individuals $\quad$ $(@bind pr  Slider(0.0:0.001:1.0,show_value=true,default=0.6666))
"""

# ╔═╡ a63a83f1-b37b-4918-804b-28d4bf04ae5a
md""" Infection Period young Individuals $(@bind ipy Slider(0:1:14,show_value=true,default=5)) days

Infection Period old Individuals $(@bind ipo Slider(0:1:14,show_value=true,default=10)) days
"""

# ╔═╡ 585df063-eebf-4367-a176-930d819b8c17
md"""---
**Execute and Plot**
"""

# ╔═╡ beee6289-5971-452b-9561-a40868f0874a
md"""---

**General Model Settings**
"""

# ╔═╡ 76c1915b-c6af-4e6d-94e8-59dc4608f8a9
begin
	#number of different types (knots)
	num_types = 6
	# 1 S | 2 V | 3 PI | 4 I | 5 R
	types = ["S" "pi" "I" "Iv" "D" "R"]

	#number of different events (arrows)
	num_events = 9
end;

# ╔═╡ a80a9b8c-c738-4f88-aefb-29fcdde1ecf5
begin #Setup Parameter
	infect_period = [ipy,ipo]
	ind = [ceil(Integer,s0*pr),floor(Integer,s0*(1-pr))]

	par = (
		α       = al,
		σ       = g,
   		ϕ       = b,
   		β       = d,
   		f_dead1 = f1,
   		f_dead2 = f2,
   		p       = Pi,
	    contact_par = [t_end,a,r_0,infect_period,ind]
	    #contact = contact
	)
end;









# ╔═╡ 1ac6bfb7-a740-4508-9615-065ba809bb0e
md"""---
**Function Library**
"""

# ╔═╡ 8471abab-a994-413c-b272-294ee818ed37
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 51a474ce-7a49-405a-8aeb-632f819557ce
Gillespie = ingredients(main_functions_path).Gillespie

# ╔═╡ 69e4276c-b804-44f4-bc5e-7298db351242
SIR = ingredients(model_path).NewTrial

# ╔═╡ eb32b9a4-dfd2-4a24-8a75-e10bf9eceda4
begin
	#Execute Simulation
	x0 = copy(initial_population)
	hist = zeros(Int,(t_end+1,num_types))

	Gillespie.run_gillespie!(
	        0:t_end,x0,par,
	        SIR.execute!,SIR.rates!,
	        Vector{Float64}(undef,num_events),hist
	)
	#Plot Simulation
	plot(hist,label=types)
end

# ╔═╡ Cell order:
# ╟─c5576946-69fd-4ac9-986f-0afb74b70137
# ╠═e4edc158-42b3-4bf1-8420-14b2bd19cd05
# ╠═774add99-f8b9-4ae5-8a08-cc6cc71f64a7
# ╟─59e73279-85bf-4c77-a207-cfd1106c1451
# ╟─85412cf0-362c-4fa5-8c37-eb261db37876
# ╠═2684693e-e731-40e0-a22d-35c62110009c
# ╟─3b6e0504-40a1-409e-a4ef-b39be2bef971
# ╟─c02d5705-f136-433a-98d4-eab309926694
# ╟─5fa15756-0c38-4f9f-b36d-2f8d3a5f9e6a
# ╟─9d5c0f86-1e41-47f4-9dea-503d3bfbd3d5
# ╟─397ac091-4c45-4fa1-9b9b-24a966661a13
# ╟─75002645-eeb0-43dd-bf53-25cc9a6b089d
# ╟─193be177-08a1-42ed-98df-8d58031cd0ae
# ╟─5ba2c008-5020-4fcd-82be-ec8089bcc6fe
# ╟─a63a83f1-b37b-4918-804b-28d4bf04ae5a
# ╟─585df063-eebf-4367-a176-930d819b8c17
# ╠═eb32b9a4-dfd2-4a24-8a75-e10bf9eceda4
# ╟─beee6289-5971-452b-9561-a40868f0874a
# ╠═76c1915b-c6af-4e6d-94e8-59dc4608f8a9
# ╠═a80a9b8c-c738-4f88-aefb-29fcdde1ecf5
# ╟─1ac6bfb7-a740-4508-9615-065ba809bb0e
# ╟─8471abab-a994-413c-b272-294ee818ed37
# ╟─51a474ce-7a49-405a-8aeb-632f819557ce
# ╟─69e4276c-b804-44f4-bc5e-7298db351242
# ╠═42e093ec-dbda-4474-bc31-8ad3bdfda9e7
