### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 90c397a8-2256-11ee-32ea-b51af7e1cac3
begin
	using Pkg
	Pkg.activate(@__DIR__)
end;

# ╔═╡ 8af9f185-2934-4560-9e5c-035267e25009
begin
	using LessUnitful
	using ExtendableGrids,GridVisualize
	using VoronoiFVM
	using LiquidElectrolytes
	using PyPlot,Colors 
	using StaticArrays
	using InteractiveUtils
	using PlutoUI
	using ForwardDiff
	using PlutoVista
	using CatmapInterface
end;

# ╔═╡ 456c01bc-a5bd-4c42-bd8a-dff19f17182b
begin
	ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
	ENV["PYTHON"] = Sys.which("python3")
	using PyCall
	#Pkg.build("PyCall")
end;

# ╔═╡ 13bdd9d6-c85d-41f7-ac68-cc74534172ac
md"""
## Microkinetic Model of the Surface Reactions

An interface to CatMap is used.
"""

# ╔═╡ b7087c4c-5e6e-4b43-8d6f-9116d6915e7a
md"""
Enable surface reactions: $(@bind allow_surfacereactions PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 7fee9783-2bb3-4825-ac95-c58c3d175fc7
mutable struct CatmapData
    rate_constant_fns
    compute_rates
end

# ╔═╡ 90f50ffd-84b1-4ebb-ba87-32f49fac9e22
function instantiate_catmap_template(template_file_path, ϕ_we)
    
    input_instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

    # replace descriptor range
    input_instance_string = replace(input_instance_string, r"descriptor_ranges.*" => "descriptor_ranges = [[$ϕ_we, $ϕ_we], [298, 298]]")

    input_instance_file_name = "CO2R.mkm"
    open(input_instance_file_name, "w") do input_instance_file
        write(input_instance_file, input_instance_string)
    end

    return input_instance_file_name
end

# ╔═╡ c8fc68d0-93fe-44ed-9b6e-ff81aa4489aa
begin
	catmap_setup_file   = "catmap_CO2R_template.mkm"
    catmap_energy_file  = "catmap_CO2R_energies.txt"
end;

# ╔═╡ 22e93bf9-7c5c-4a55-bbf3-06623bc37b66
@unitfactors mol dm eV μA μF cm μm bar;

# ╔═╡ 2db7bc4c-e768-4093-abee-f6545c4d9d1c
begin
	@phconstants c_0 N_A e k_B
	F = N_A * e
end;

# ╔═╡ 62238672-6989-46d5-bea9-cda73c398f34
md"""
## System Parameters
"""

# ╔═╡ f9f9ff55-4cd1-404e-9fff-a227f1e72872
pH = 6.8;

# ╔═╡ 3cbf0396-336a-4cd3-8301-0776b8a9e8db
T = 273.15 + 25 * ufac"K";

# ╔═╡ 1adc2d10-8a14-4679-8c55-c748f7ee39d2
begin
	C_gap = 20 * μF/cm^2
	ϕ_pzc = 0.16 * ufac"V"
end;

# ╔═╡ cd9d945e-418a-4fc6-9c48-b839d56c2925
md"""
#### Diffusion constants
"""

# ╔═╡ e48990a3-f734-487f-b32c-72a56378a409
D = [1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s"; # from Ringe paper

# ╔═╡ 276b3a57-333c-4f7d-80f7-34b95a256076
md"""
Henry's Law
"""

# ╔═╡ 7b88963b-500c-4a6a-bb1b-e5bf05275e74
begin # Henry constants
	Hcp_CO  = 9.7e-6 * ufac"mol/(m^3 * Pa)"
	Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"
end;

# ╔═╡ ac612516-58b9-4cdd-9fb4-389ef1c17402
md"""
## Kinetic Model for the Buffer Reactions
"""

# ╔═╡ 19b9d9e2-a0e3-4da6-96de-3ea180fd00e8
md"""
Consider the bicarbonate buffer system in base:

$CO_2 + OH^- \rightleftharpoons HCO_3^-$

$HCO_3^- + OH^- \rightleftharpoons CO_3^{2-} + H_2O$

and acid:

$CO_2 + H_2O \rightleftharpoons HCO_3^- + H^+$

$HCO_3^- \rightleftharpoons CO_3^{2-} + H^+$

Moreover, the autoprotolysis of water is considered:

$H_2O \rightleftharpoons H^+ + OH^-$

"""

# ╔═╡ 1e1a135c-b3c2-439b-9647-de369809a785
begin # bulk species
	ikplus      = 1
    ihco3       = 2
    ico3        = 3
    ico2        = 4
    ico         = 5
    iohminus    = 6
    ihplus      = 7
	nc 			= 7
end;

# ╔═╡ 96dd7667-db0d-458b-9dec-cb95f9f7a0dd
md"""
### Reaction Constants
"""

# ╔═╡ b72b2600-7a8c-4fad-b3fd-c955ad1953a2
begin
	# buffer equations
	## in base
	## CO2 + OH- <=> HCO3-
	kbe1 = 4.44e7 / (mol/dm^3)
	kbf1 = 5.93e3 / (mol/dm^3) / ufac"s"
	kbr1 = kbf1 / kbe1
	## HCO3- + OH- <=> CO3-- + H2O
	kbe2 = 4.66e3 / (mol/dm^3)
	kbf2 = 1.0e8 / (mol/dm^3) / ufac"s"
	kbr2 = kbf2 / kbe2
	
	## in acid
	## CO2 + H20 <=> HCO3- + H+
	kae1 = 4.44e-7 * (mol/dm^3)
	kaf1 = 3.7e-2 / ufac"s"
	kar1 = kaf1 / kae1
	## HCO3- <=> CO3-- + H+ 
	kae2 = 4.66e-5 / (mol/dm^3)
	kaf2 = 59.44e3 / (mol/dm^3) / ufac"s"
	kar2 = kaf2 / kae2
	## autoprotolyse
	kwe  = 1.0e-14 * (mol/dm^3)^2
	kwf  = 2.4e-5 * (mol/dm^3) / ufac"s"
	kwr  = kwf / kwe
end;

# ╔═╡ d783cd12-e7b5-499f-924c-89ca5e59fcf1
function reaction(f, u::VoronoiFVM.NodeUnknowns{Tv, Tc, Tp, Ti}, node, data) where {Tv, Tc, Tp, Ti}  
	# buffer reactions
	rates       = zeros(Tv, 5)
	## in base
	## CO2 + OH- <=> HCO3-
	rates[1]    = kbf1 * u[ico2] * u[iohminus]  - kbr1 * u[ihco3]  
	## HCO3- + OH- <=> CO3-- + H2O
	rates[2]    = kbf2 * u[ihco3] * u[iohminus] - kbr2 * u[ico3]

	## in acid
	## CO2 + H20 <=> HCO3- + H+
	rates[3]    = kaf1 * u[ico2] - kar1 * u[ihco3] * u[ihplus]  
	## HCO3- <=> CO3-- + H+ 
	rates[4]    = kaf2 * u[ihco3] - kar2 * u[ico3] * u[ihplus]  

	## autoprotolyse
	rates[5]    = kwf - kwr * u[ihplus] * u[iohminus]  

	f[ihco3]    -= rates[1] - rates[2] + rates[3] - rates[4]
	f[ico3]     -= rates[2] + rates[4]
	f[ihplus]   -= rates[3] + rates[4] + rates[5]
	f[iohminus] -= -rates[1] -rates[2] + rates[5]

	nothing
end;

# ╔═╡ 4b12420f-b7b1-48e2-949b-2028f43702ee
begin # surface species
	ico_ad      = 8
    ico2_ad     = 9
    icooh_ad    = 10
	na 			= 3
end;

# ╔═╡ 31d33cb2-6814-444f-8358-9d84588a3e4d
begin
  	instance_input_file_name = instantiate_catmap_template(catmap_setup_file, 0.0)
    
	~, rate_constant_fns, ~, compute_rates  = get_catmap_output(
		instance_input_file_name, 
		catmap_energy_file, 
		true
	)
	catmap_data = CatmapData(rate_constant_fns, compute_rates)    
end;

# ╔═╡ 6ab26fc1-9e89-439c-8fb6-0463b424cdad
function we_breactions(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
	(; iϕ, ϕ_we) = data

	
	σ 		= C_gap / (μF/cm^2) * (ϕ_we - u[iϕ] - ϕ_pzc)
	rates 	= catmap_data.compute_rates(
		catmap_data.rate_constant_fns[1], 
		[u[ico2_ad], u[icooh_ad], u[ico_ad]], 
		[u[ico2] / Hcp_CO2 / bar, u[ico] / Hcp_CO / bar, 1.0, 1.0, u[iohminus], 1.0], σ
	)

	S   = 9.61e-5 / N_A * (1.0e10)^2 * ufac"mol/m^2"

	# bulk species
	f[ico] 		+= -rates[4] * S
	f[ico2] 	+= rates[1] * S
	f[iohminus] += -rates[2] * S - rates[3] * S
	
	# surface species
	f[ico2_ad] 	+= -rates[1] + rates[2]
	f[ico_ad] 	+= -rates[3] + rates[4]
	f[icooh_ad] += -rates[2] + rates[3]

	nothing
end;

# ╔═╡ 5dda5494-51e6-47b0-9e32-b59fc7126904
function halfcellbc(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
		(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)

		#boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
		# Robin b.c. for the Poisson equation
		boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))

        if bnode.region==Γ_we
			if allow_surfacereactions
				we_breactions(f, u, bnode, data)
			end
        end
        nothing
end;

# ╔═╡ 4dec7b8d-9ab2-40c0-a48d-cb358f9c49a4
function pre(sol, t)
	
	instance_input_file_name = instantiate_catmap_template(catmap_setup_file, -t)
	
	~, rate_constant_fns, ~, compute_rates  = get_catmap_output(
		instance_input_file_name, 
		catmap_energy_file, 
		true
	)
	
	catmap_data.rate_constant_fns           = rate_constant_fns
	catmap_data.compute_rates               = compute_rates

	nothing
end;

# ╔═╡ 2edd8816-3c9f-415c-b78a-225ea6b9cdc0
md"""
## Simulation of the $CO_2$ reduction
"""

# ╔═╡ 5e8cc51e-18c1-429d-a173-9fbbd14190ef


# ╔═╡ 1ab940f5-d5ea-4c13-a282-1e1a758920c8
function simulate_CO2R(;nref 		= 0,
              			voltages 	= (-1.5:0.1:-0.0)*ufac"V",
              			scheme 		= :μex,
              			κ 			= 10.0,
              			kwargs...)

    defaults = (; 	max_round 	= 3,
              		tol_round 	= 1.0e-9,
              		verbose 	= "e",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)

    kwargs 	= merge(defaults, kwargs) 

    hmin    = 1.0e-6 * μm * 2.0^(-nref)
    hmax    = 1.0 * μm * 2.0^(-nref)
    L       = 80.0 * μm
    X       = geomspace(0,L,hmin,hmax)
    grid    = simplexgrid(X)


    celldata 	= ElectrolyteData(; nc 		= 7,
                             		na 		= 3,
                             		z 		= [1,-1,-2,0,0,-1,1],
                             		D 		= D,
                             		T 		= T,
                             		eneutral=false,
                             		κ 		=fill(κ,7),
                             		Γ_we 	= 1,
                             		Γ_bulk 	= 2,
                             		scheme)
    
    celldata.c_bulk[ikplus]         = 0.0 * mol/dm^3
    celldata.c_bulk[ihco3]          = 0.091 * mol/dm^3
    celldata.c_bulk[ico3]           = 2.68e-5 * mol/dm^3
    celldata.c_bulk[ico2]           = 0.033 * mol/dm^3
    celldata.c_bulk[ico]            = 0.0 * mol/dm^3
    celldata.c_bulk[iohminus]       = 10^(pH - 14) * mol/dm^3
    celldata.c_bulk[ihplus]         = 10^(-pH) * mol/dm^3

    celldata.c_bulk[ikplus]         = -celldata.c_bulk' * celldata.z

    @assert isapprox(celldata.c_bulk' * celldata.z, 0, atol = 1.0e-10)
    
    cell 	= PNPSystem(grid; bcondition=halfcellbc, reaction=reaction, celldata)     
    ivresult= ivsweep(cell; voltages, more_pre=pre, more_post=((x...) -> nothing), store_solutions=true, kwargs...)

	cell, ivresult
end


# ╔═╡ 7ce683b4-b876-4c7f-9d51-809fb1580c50
(cell, ivresult) = simulate_CO2R(; κ = 4.0);

# ╔═╡ 1723ade4-1073-49de-b19c-cf045ad32e4b
md"""
## Visualization
"""

# ╔═╡ 2d89a813-9ed5-4ab4-875f-3ff6f7e234e6
md"""
Choose applied voltage: $(@bind ϕ_we_index PlutoUI.Slider(1:10:length(ivresult.voltages)))
"""

# ╔═╡ d25ee1bb-d0f0-4955-8d9c-d2166f4c1c8e
md"""
Potential at working electrode = $(ivresult.voltages[ϕ_we_index]) 
"""

# ╔═╡ 82e95927-9e09-4fda-8a34-e609fe1e4337
begin
	vis = GridVisualizer(; Plotter = PlutoVista, layout=(1,2), resolution = (1200, 400))
	
	# current-voltage plot
    currs = [j[iohminus] * F for j in ivresult.j_we]
	
	scalarplot!(
			vis[1,1], 
			ivresult.voltages[ivresult.voltages .< -0.6], 
			currs[ivresult.voltages .< -0.6] * ufac"cm^2/mA",
			color 		= "red",
			markershape = :utriangle,
			markersize 	= 7, 
			markevery 	= 10,
			label 		= "PNP",
			clear 		= true,
			legend 		= :lt,
			xlabel 		= "Δϕ [V]",
			ylabel 		= "I [mA/cm^2]", 
			yscale 		= :log
		)
	
	pHs 		= -log10.(max.(0.0, ivresult.solutions[ϕ_we_index][ihplus, :]) / (mol/dm^3))
    scalarplot!(
		vis[1,2], 
		cell.grid, 
		pHs, 
		xlabel 	="Distance from working electrode [m]", 
		ylabel 	="pH-Value", 
		xscale 	= :log,
	)
	
	reveal(vis)
end

# ╔═╡ Cell order:
# ╠═90c397a8-2256-11ee-32ea-b51af7e1cac3
# ╠═8af9f185-2934-4560-9e5c-035267e25009
# ╠═456c01bc-a5bd-4c42-bd8a-dff19f17182b
# ╟─13bdd9d6-c85d-41f7-ac68-cc74534172ac
# ╟─b7087c4c-5e6e-4b43-8d6f-9116d6915e7a
# ╠═7fee9783-2bb3-4825-ac95-c58c3d175fc7
# ╠═90f50ffd-84b1-4ebb-ba87-32f49fac9e22
# ╠═c8fc68d0-93fe-44ed-9b6e-ff81aa4489aa
# ╠═22e93bf9-7c5c-4a55-bbf3-06623bc37b66
# ╠═2db7bc4c-e768-4093-abee-f6545c4d9d1c
# ╟─62238672-6989-46d5-bea9-cda73c398f34
# ╠═f9f9ff55-4cd1-404e-9fff-a227f1e72872
# ╠═3cbf0396-336a-4cd3-8301-0776b8a9e8db
# ╠═1adc2d10-8a14-4679-8c55-c748f7ee39d2
# ╟─cd9d945e-418a-4fc6-9c48-b839d56c2925
# ╠═e48990a3-f734-487f-b32c-72a56378a409
# ╟─276b3a57-333c-4f7d-80f7-34b95a256076
# ╠═7b88963b-500c-4a6a-bb1b-e5bf05275e74
# ╟─ac612516-58b9-4cdd-9fb4-389ef1c17402
# ╟─19b9d9e2-a0e3-4da6-96de-3ea180fd00e8
# ╠═1e1a135c-b3c2-439b-9647-de369809a785
# ╟─96dd7667-db0d-458b-9dec-cb95f9f7a0dd
# ╠═b72b2600-7a8c-4fad-b3fd-c955ad1953a2
# ╠═d783cd12-e7b5-499f-924c-89ca5e59fcf1
# ╠═4b12420f-b7b1-48e2-949b-2028f43702ee
# ╠═6ab26fc1-9e89-439c-8fb6-0463b424cdad
# ╠═5dda5494-51e6-47b0-9e32-b59fc7126904
# ╠═31d33cb2-6814-444f-8358-9d84588a3e4d
# ╠═4dec7b8d-9ab2-40c0-a48d-cb358f9c49a4
# ╟─2edd8816-3c9f-415c-b78a-225ea6b9cdc0
# ╠═5e8cc51e-18c1-429d-a173-9fbbd14190ef
# ╠═1ab940f5-d5ea-4c13-a282-1e1a758920c8
# ╠═7ce683b4-b876-4c7f-9d51-809fb1580c50
# ╟─1723ade4-1073-49de-b19c-cf045ad32e4b
# ╟─2d89a813-9ed5-4ab4-875f-3ff6f7e234e6
# ╟─d25ee1bb-d0f0-4955-8d9c-d2166f4c1c8e
# ╠═82e95927-9e09-4fda-8a34-e609fe1e4337
