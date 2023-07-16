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

# ╔═╡ 19bec5c0-23ca-11ee-00fb-c3a65c901c37
begin
	using Pkg
	Pkg.activate(@__DIR__)
end;

# ╔═╡ 5efbbf3d-9f95-47f4-bc18-4888cce1f1f7
begin
	using LessUnitful
	using ExtendableGrids,GridVisualize
	using VoronoiFVM
	using LiquidElectrolytes 
	using StaticArrays
	using InteractiveUtils
	using PlutoUI
	using ForwardDiff
	using PlutoVista
end;

# ╔═╡ d19615b1-e862-4a0d-b087-a4f07f6b4982
begin
	ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
	ENV["PYTHON"] = Sys.which("python3")
	using PyCall
	#Pkg.build("PyCall")
end;

# ╔═╡ c3575a60-3cc8-4f2a-a378-53140ca2d582
md"""
### Some Help from Python
"""

# ╔═╡ c75e1263-2ee2-48e1-8ecd-c499e5d4388d
py"""
from ase.thermochemistry import HarmonicThermo

def get_thermal_correction_adsorbate(T, frequencies):
	thermo = HarmonicThermo(frequencies)
	return thermo.get_helmholtz_energy(T, verbose=False)
"""

# ╔═╡ 98153553-0551-4d4f-abad-d866b3ff3b1d
md"""
### Units and Constants
"""

# ╔═╡ 5a37f25b-9b9d-4f38-8c71-b9639bea74fd
@unitfactors cm μF mol dm bar Pa μA μm eV m;

# ╔═╡ 22c8693e-49be-45fe-9bc1-5a1eccba58a6
begin
	@phconstants N_A e R ε_0 c_0 h
	F = N_A * e
end;

# ╔═╡ 1b28939e-1ffd-4e8d-9f04-0ac7e125f844
md"""
### System Parameters
"""

# ╔═╡ 40bdb2bd-44f1-4ea5-9187-0a43afabd37b
pH = 6.8;

# ╔═╡ b1bede24-45e2-4d80-8d01-f68e9c9fedb1
T = 273.15 + 25 * ufac"K";

# ╔═╡ fdf5263b-8393-409e-97bf-6e4ad6870299
D = [1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s"; # from Ringe paper

# ╔═╡ 9890dc42-2f3e-4da3-a628-9cb1d93441b4
md"""
Henry's Law
"""

# ╔═╡ 06856204-9aec-40e3-9f22-7bb33abdbd5a
Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"

# ╔═╡ 668d1f5d-1a17-44d2-a676-6eb7df2acd49
md"""
## Kinetic Model for the Buffer Reactions
"""

# ╔═╡ a7bfda5f-9414-4259-a635-ad9831efcd5c
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

# ╔═╡ ae7c1791-b548-4e27-8e6a-a2b458e90b81
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

# ╔═╡ fec21ff8-4434-4379-a3b3-9144a7654d86
md"""
### Reaction Constants
"""

# ╔═╡ d9968efe-7ba1-4520-9a66-395dd9c4398e
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

# ╔═╡ 8c3e8dc8-a205-4a5a-8f1e-fb15f0065272
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

# ╔═╡ 7f19d93d-f275-42ee-8ba3-cc3f3bf5f99f
md"""
## Kinetic Model for Surface Reactions: Rate-determining Step
"""

# ╔═╡ 68c5cfed-7fb3-41c3-bb68-d76c7866383e
E_ads_CO2 = 0.657600203 * eV;

# ╔═╡ 0a03fb7b-2fa9-48bf-987b-4b8d6ca5d61a
md"""
#### Harmonic Adsorbate Correction
"""

# ╔═╡ 6d0de240-935e-45e7-b94a-a77d898097cb
frequencies = [136.85, 183.6, 212.95, 250.7, 306.0, 510.55, 562.25, 1176.05, 1889.85] * h * c_0 / cm / eV;

# ╔═╡ a5e0c589-1794-4493-a316-5b7a8ac22705
harmonic_adsorbate_correction = py"get_thermal_correction_adsorbate"(T, frequencies) * eV

# ╔═╡ 0d06826e-710e-4ea3-9157-ca6b21f10640
md"""
#### Electrochemical Corrections

The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$

The surface charging relation $σ = σ(U)$ is given by the Robin boundary condition

$σ(U) = C_{gap} (ϕ_{we} - ϕ_{pzc} - ϕ^\ddagger)$

where the gap capacitance between the working electrode and the reaction plane ($\ddagger$) is given by $C_{gap} = 20~μF/cm^2$. The potential of zero current is measured to be $ϕ_{pzc} = 0.16~V$.
"""

# ╔═╡ fd1315f2-e996-4695-b6b1-9dbbb4d4fd28
begin
	C_gap = 20 * μF/cm^2
	ϕ_pzc = 0.16 * ufac"V"
end;

# ╔═╡ a7f36422-7d2f-4646-b250-e86931d42b6f
electro_corrections_params = [-0.000286600929 / (μA/cm^2)^2, 0.0297720125 / (μA/cm^2)];

# ╔═╡ f34235fc-5a1d-4a30-8188-c32d8e5b908f
md"""
### Boundary Condition: Surface Reaction
"""

# ╔═╡ 83d7afb2-332e-4548-a287-13ed3fef5b8c
md"""
Enable surface reactions: $(@bind allow_surfacereactions PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 47e8a25f-6390-4794-bf08-9cf7719c2b9a
function we_breactions(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
	(; iϕ, ϕ_we, RT) = data
	
	# Flux conditions for CO2 and CO
	prefactor   = 1.0e13
	γ_CO2       = 1.0
	a_CO2       = γ_CO2 * u[ico2] / Hcp_CO2 / bar
	θ_free      = 0.9999
	S           = 9.61e-5 / N_A * (1.0e10)^2 * mol/m^2
	σ           = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)

	# compute adsorption energy 
	electrochemical_correction = electro_corrections_params' * [σ^2, σ] * eV
	ΔG  = (E_ads_CO2 + harmonic_adsorbate_correction + electrochemical_correction) * N_A            
	
	# reaction rate
	r = prefactor * a_CO2 * θ_free * exp(-ΔG / RT - 2.3 * pH) * S
	
	f[ico2] 	= r
	f[ico]  	= -r
	f[iohminus] = -2 * r

	nothing
end

# ╔═╡ a1781717-d9ab-4731-a082-d7691f409738
function halfcellbc(f,u,bnode,data)
        (; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

        bulkbcondition(f, u, bnode, data; region=Γ_bulk)
		#boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
		# Robin b.c. for the Poisson equation
		boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))

        if bnode.region==Γ_we
			if allow_surfacereactions
				we_breactions(f, u, bnode, data)
			end
        end
        nothing
    end

# ╔═╡ ba725d6f-0271-4322-9405-42c81a045ef2
function simulate_CO2R(;nref 	=0,
              			voltages=(-1.5:0.1:-0.0)*ufac"V",
              			scheme 	= :μex,
              			κ 		= 10.0,
              			kwargs...)
   
    defaults =(;max_round 	= 3,
				tol_round 	= 1.0e-9,
				verbose 	= "e",
				reltol 		= 1.0e-8,
				tol_mono 	= 1.0e-10)

    kwargs 	= merge(defaults, kwargs) 

	# grid
    hmin 	= 1.0e-6 * μm * 2.0^(-nref)
    hmax 	= 1.0 * μm * 2.0^(-nref)
    L 		= 80.0 * μm
    X 		= geomspace(0, L, hmin, hmax)
    grid 	= simplexgrid(X)

	# electrolyte
    celldata = ElectrolyteData(;nc 		= 7,
                             	z 		= [1, -1, -2, 0, 0, -1, 1],
                             	D 		= D, # from Ringe paper
                             	T 		= T,
                             	eneutral= false,
                             	κ 		= fill(κ,7),
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
    
    cell        = PNPSystem(grid; bcondition=halfcellbc, reaction=reaction, celldata)
    ivresult    = ivsweep(cell; voltages, store_solutions=true, more_pre=((x...)->nothing), more_post=((x...)->nothing), kwargs...)

	cell, ivresult
end

# ╔═╡ e64f3ebb-4654-4569-9cf9-76bb38d1b133
(cell, ivresult) = simulate_CO2R(; κ = 4.0);

# ╔═╡ 52f8c956-a04b-466f-80cb-c62b92cdee31
md"""
## Visualization
"""

# ╔═╡ 4fb1bf76-936d-4603-b353-6489ef447211
md"""
Choose applied voltage: $(@bind ϕ_we_index PlutoUI.Slider(1:10:length(ivresult.voltages)))
"""

# ╔═╡ c652b5df-5693-44c5-b0e6-c6e23ebb70aa
md"""
Potential at working electrode = $(ivresult.voltages[ϕ_we_index]) 
"""

# ╔═╡ c1288f6a-7e10-4f30-aaf8-5799393b1885
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
	
	pHs 		= -log10.(ivresult.solutions[ϕ_we_index][ihplus, :] / (mol/dm^3))
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
# ╠═19bec5c0-23ca-11ee-00fb-c3a65c901c37
# ╠═5efbbf3d-9f95-47f4-bc18-4888cce1f1f7
# ╟─c3575a60-3cc8-4f2a-a378-53140ca2d582
# ╠═d19615b1-e862-4a0d-b087-a4f07f6b4982
# ╠═c75e1263-2ee2-48e1-8ecd-c499e5d4388d
# ╟─98153553-0551-4d4f-abad-d866b3ff3b1d
# ╠═5a37f25b-9b9d-4f38-8c71-b9639bea74fd
# ╠═22c8693e-49be-45fe-9bc1-5a1eccba58a6
# ╟─1b28939e-1ffd-4e8d-9f04-0ac7e125f844
# ╠═40bdb2bd-44f1-4ea5-9187-0a43afabd37b
# ╠═b1bede24-45e2-4d80-8d01-f68e9c9fedb1
# ╠═fdf5263b-8393-409e-97bf-6e4ad6870299
# ╟─9890dc42-2f3e-4da3-a628-9cb1d93441b4
# ╠═06856204-9aec-40e3-9f22-7bb33abdbd5a
# ╟─668d1f5d-1a17-44d2-a676-6eb7df2acd49
# ╟─a7bfda5f-9414-4259-a635-ad9831efcd5c
# ╠═ae7c1791-b548-4e27-8e6a-a2b458e90b81
# ╟─fec21ff8-4434-4379-a3b3-9144a7654d86
# ╠═d9968efe-7ba1-4520-9a66-395dd9c4398e
# ╠═8c3e8dc8-a205-4a5a-8f1e-fb15f0065272
# ╟─7f19d93d-f275-42ee-8ba3-cc3f3bf5f99f
# ╠═68c5cfed-7fb3-41c3-bb68-d76c7866383e
# ╟─0a03fb7b-2fa9-48bf-987b-4b8d6ca5d61a
# ╠═6d0de240-935e-45e7-b94a-a77d898097cb
# ╠═a5e0c589-1794-4493-a316-5b7a8ac22705
# ╟─0d06826e-710e-4ea3-9157-ca6b21f10640
# ╠═fd1315f2-e996-4695-b6b1-9dbbb4d4fd28
# ╠═a7f36422-7d2f-4646-b250-e86931d42b6f
# ╟─f34235fc-5a1d-4a30-8188-c32d8e5b908f
# ╟─83d7afb2-332e-4548-a287-13ed3fef5b8c
# ╠═47e8a25f-6390-4794-bf08-9cf7719c2b9a
# ╠═a1781717-d9ab-4731-a082-d7691f409738
# ╠═ba725d6f-0271-4322-9405-42c81a045ef2
# ╠═e64f3ebb-4654-4569-9cf9-76bb38d1b133
# ╟─52f8c956-a04b-466f-80cb-c62b92cdee31
# ╟─4fb1bf76-936d-4603-b353-6489ef447211
# ╟─c652b5df-5693-44c5-b0e6-c6e23ebb70aa
# ╠═c1288f6a-7e10-4f30-aaf8-5799393b1885
