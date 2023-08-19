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

# ╔═╡ a532a10a-9cb8-44a9-89f4-f3517177b602
begin
	using Pkg
	Pkg.activate(@__DIR__)
end;

# ╔═╡ 99c3f409-dd51-4694-8f78-87a2aca41d40
begin
	using LessUnitful
	using ExtendableGrids,GridVisualize
	using PlutoVista
	using PlutoUI
	using LiquidElectrolytes
	using VoronoiFVM
	using PyCall
	using CatmapInterface
	using ForwardDiff
	using DataFrames
end;

# ╔═╡ 92be2111-381e-4260-a514-1632c5a38355
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 53c1639f-58e0-4b28-b4b6-34d7fd1f969f
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA;

# ╔═╡ 02aefeeb-4ebb-4daf-b06f-80d537409844
@phconstants N_A c_0 k_B e;

# ╔═╡ 306dc5b4-29f6-11ee-3112-874796d88919
md"""
# Simulation of a $CO_2$ Reduction Mechanism
"""

# ╔═╡ aaa10928-59d2-4651-b9d1-c71bb672ae95
md"""
## Geometry
"""

# ╔═╡ 4d283cbb-6bd7-4f84-9e1b-9d5f0819aa4d
begin 
	hmin 	= 1.0e-6 	* μm
	hmax 	= 1.0 		* μm 
	L 		= 80.0 		* μm 
	X 		= geomspace(0, L, hmin, hmax)
    
	grid 	= simplexgrid(X)

	Γ_we 	= 1
	Γ_bulk  = 2

	gridplot(grid; Plotter = PlutoVista, resolution = (700, 200))
end

# ╔═╡ ebbd9eb3-7dfc-417b-a39c-7023c15361c2
md"""
## The Electrolyte
"""

# ╔═╡ 57bdf775-9974-4f1f-8466-37f28b0c80b5
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

# ╔═╡ 43b24817-648d-4e0c-a885-f5445a726b5e
pH = 6.8;

# ╔═╡ 58e784bc-f838-4ad2-87f5-c868313c0f91
T = 273.15 + 25 * K;

# ╔═╡ e3a2076f-f37b-4732-9238-49209253e3b8
struct BulkSpecies
	name::String
	index::Int32
	charge::Int32
	diffusivity::Float64
	κ::Float64
	c_bulk::Float64
end

# ╔═╡ 5391f6fc-ba0b-4fd1-bf4c-d26b38bd12dc
begin
	# Henry constants
    Hcp_CO  = 9.7e-6 * mol/(m^3 * Pa)
    Hcp_CO2 = 3.3e-4 * mol/(m^3 * Pa)
end;

# ╔═╡ e91db1f6-1388-4262-b0bc-235cae1867eb
md"""
## Surface Species
"""

# ╔═╡ bb331dd8-3990-4cb5-87d4-ed4a93e54098
S = 9.61e-5 / N_A * (1.0e10)^2 * mol/m^2;

# ╔═╡ f59ddfde-3eef-485c-b6a6-7bd050eb5b64
struct SurfaceSpecies
	name::String
	index::Int32
	initial_coverage::Float64
end;

# ╔═╡ 4a8da214-c8b0-4b73-8f14-3010d1bd6756
begin
	surface = Dict(
		"CO" 	=> SurfaceSpecies("CO", 8, 0.0),
		"CO₂" 	=> SurfaceSpecies("CO₂", 9, 0.0),
		"COOH" 	=> SurfaceSpecies("COOH", 10, 0.0)
	)

	na = length(surface)
end;

# ╔═╡ deffe3ec-6580-4222-9c9d-075e0a7cb958
md"""
Enable surface reactions: $(@bind allow_surfacereactions PlutoUI.CheckBox(default=true))
"""

# ╔═╡ a9d434d9-7f7e-4c2f-95fd-3e24476e34ce
md"""
## The Half-Cell
"""

# ╔═╡ c28be48d-1c3b-4581-bfa7-41c1ed61b51e
md"""
## Kinetic Model for the Buffer Systems
"""

# ╔═╡ ecc38283-877e-4a1b-9a0c-451c024c9683
md"""
### Reaction Constants
"""

# ╔═╡ 278cf5a7-da40-4106-b28e-b701359510d6
begin
	# buffer equations
	## in base
	## CO2 + OH- <=> HCO3-
	kbe1 = 4.44e7 / (mol/dm^3)
	kbf1 = 5.93e3 / (mol/dm^3) / s
	kbr1 = kbf1 / kbe1
	## HCO3- + OH- <=> CO3-- + H2O
	kbe2 = 4.66e3 / (mol/dm^3)
	kbf2 = 1.0e8 / (mol/dm^3) / s
	kbr2 = kbf2 / kbe2
	
	## in acid
	## CO2 + H20 <=> HCO3- + H+
	kae1 = 4.44e-7 * (mol/dm^3)
	kaf1 = 3.7e-2 / s
	kar1 = kaf1 / kae1
	## HCO3- <=> CO3-- + H+ 
	kae2 = 4.66e-5 / (mol/dm^3)
	kaf2 = 59.44e3 / (mol/dm^3) / s
	kar2 = kaf2 / kae2
	## autoprotolyse
	kwe  = 1.0e-14 * (mol/dm^3)^2
	kwf  = 2.4e-5 * (mol/dm^3) / s
	kwr  = kwf / kwe
end;

# ╔═╡ df16bc86-8041-402c-8b94-35ef5aac961b
md"""
### Mean-field Approach
"""

# ╔═╡ 04c5a531-f7bb-4a93-a619-2afd92586e2a
md"""
## Microkinetic Model for the Surface Reactions
"""

# ╔═╡ faedb03c-c4ee-45cc-9263-0fd8eaf09850
md"""
The reaction mechanism for the $CO_2$ reduction is divided into four elementary reactions at the electrode surface:

1. Adsorption of $CO_2$ molecules at the oxygen atoms
${CO_2}_{(aq)} + 2* \rightleftharpoons {CO_2 *}_{(ad)}$

2. First proton-coupled electron transfer
${CO_2*}_{(ad)} + H_2O_{(l)} + e^- \rightleftharpoons COOH*_{(ad)} + OH^-_{(aq)}$

3. Second proton-coupled electron transfer
$COOH*_{(ad)} + H_2O_{(l)} + e^- \rightleftharpoons CO*_{(ad)} + H_2O_{(l)} + OH^{-}_{(aq)} + *$
with the transition state: $*CO-OH^{TS}$

4. Desorption of $CO$
$*CO_{(ad)} \rightleftharpoons CO_{(aq)} + *$
"""

# ╔═╡ aed5431b-a89c-46c6-be7c-2d1c6fcd00c7
md"""
Choose calculation mode of the rate constants: $(@bind rate_const_mode PlutoUI.Select(["Explicit", "CatMAP"]))
"""

# ╔═╡ fb3df4df-2a66-4796-a7de-be46e8c8fd29
md"""
## The Simulation
"""

# ╔═╡ 11efa591-a50d-451b-b9bf-31c8d5c87a7c
md"""
## Visualization
"""

# ╔═╡ 20114f4a-527f-4782-8b67-8bb8fe824319
md"""
## Some Help
"""

# ╔═╡ c64685fb-4c1c-435d-af95-7dafdfc77e15
function is_electroneutral(bulk_species)
	isapprox(
		sum([bs.c_bulk * bs.charge for bs in values(bulk_species)]), 0.0; 
		atol = 1.0e-5
	)
end;

# ╔═╡ c0c3429a-a426-4623-aec7-977b7f9cb685
begin         #name #index #charge #diffusion #solvation number #bulk concentration
	bulk = Dict( 
	"K⁺"  => BulkSpecies("K⁺", 	 1,  1, 1.957e-9 * m^2/s, 8.0, 0.0910535 * mol/dm^3),
	"HCO₃⁻"=>BulkSpecies("HCO₃⁻",2, -1, 1.185e-9 * m^2/s, 4.0, 0.091     * mol/dm^3),
	"CO₃²⁻"=>BulkSpecies("CO₃²⁻",3, -2, 0.923e-9 * m^2/s, 4.0, 2.68e-5   * mol/dm^3),
	"CO₂" => BulkSpecies("CO₂",  4,  0, 1.91e-9  * m^2/s, 4.0, 0.033     * mol/dm^3),
	"CO"  => BulkSpecies("CO", 	 5,  0, 2.23e-9  * m^2/s, 4.0, 0.0       * mol/dm^3),
	"OH⁻" => BulkSpecies("OH⁻",  6, -1, 5.273e-9 * m^2/s, 4.0, 10^(pH-14)*mol/dm^3),
	"H⁺"  => BulkSpecies("H⁺", 	 7,  1, 9.310e-9 * m^2/s, 4.0, 10^(-pH)  * mol/dm^3)
	)
	
	@assert is_electroneutral(bulk)
	
	nc = length(bulk)
end;

# ╔═╡ 4896bcdb-8661-4b8c-aaf3-5d2ae6eee1d1
bulk_ordered = sort(collect(values(bulk)); by = x -> x.index);

# ╔═╡ a6e0b8fe-fa45-4c2e-a1c2-7c2710a65988
celldata = ElectrolyteData(;nc    = nc,
						  	na    = allow_surfacereactions ? na : 0,
						  	z     = getproperty.(bulk_ordered, :charge),
						  	D     = getproperty.(bulk_ordered, :diffusivity),
						  	T     = T,
						  	eneutral=false,
						  	κ     = getproperty.(bulk_ordered, :κ),
                            c_bulk= getproperty.(bulk_ordered, :c_bulk),
						  	Γ_we  = Γ_we,
						  	Γ_bulk= Γ_bulk,
						  	scheme= :μex);

# ╔═╡ 2a7a3783-846d-4715-960e-979b72ded8d9
function reaction(f, u::VoronoiFVM.NodeUnknowns{Tv, Tc, Tp, Ti}, node, data) where {Tv, Tc, Tp, Ti}  
	# buffer reactions
	rates       = zeros(Tv, 5)
	## in base
	## CO2 + OH- <=> HCO3-
	rates[1]    = kbf1 * u[bulk["CO₂"].index] * u[bulk["OH⁻"].index] 
	rates[1]   -= kbr1 * u[bulk["HCO₃⁻"].index]  
	## HCO3- + OH- <=> CO3-- + H2O
	rates[2]    = kbf2 * u[bulk["HCO₃⁻"].index] * u[bulk["OH⁻"].index] 
	rates[2]   -= kbr2 * u[bulk["CO₃²⁻"].index]

	## in acid
	## CO2 + H20 <=> HCO3- + H+
	rates[3]    = kaf1 * u[bulk["CO₂"].index] 
	rates[3]   -= kar1 * u[bulk["HCO₃⁻"].index] * u[bulk["H⁺"].index]
				
	## HCO3- <=> CO3-- + H+ 
	rates[4]    = kaf2 * u[bulk["HCO₃⁻"].index] 
	rates[4]   -= kar2 * u[bulk["CO₃²⁻"].index] * u[bulk["H⁺"].index] 

	## autoprotolyse
	rates[5]    = kwf 
	rates[5]   -= kwr * u[bulk["H⁺"].index] * u[bulk["OH⁻"].index]  

	f[bulk["HCO₃⁻"].index] 	-= rates[1] - rates[2] + rates[3] - rates[4]
	f[bulk["CO₃²⁻"].index] 	-= rates[2] + rates[4]
	f[bulk["H⁺"].index]   	-= rates[3] + rates[4] + rates[5]
	f[bulk["OH⁻"].index] 	-= -rates[1] -rates[2] + rates[5]

	nothing
end;

# ╔═╡ 8905005c-0d88-4e2b-ae3e-cbee952d9a86
md"""
### From Python
"""

# ╔═╡ 4fda1782-4f70-4007-9ac3-12c931c4bd54
begin
	py"""
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.build import molecule

def get_thermal_correction_adsorbate(T, frequencies):
    thermo = HarmonicThermo(frequencies)
    return thermo.get_helmholtz_energy(T, verbose=False)



def get_thermal_correction_ideal_gas(T, frequencies, symmetrynumber, geometry, spin, name):
    thermo = IdealGasThermo(frequencies, geometry, atoms=molecule(name), symmetrynumber=symmetrynumber, spin=spin) 
    H = thermo.get_enthalpy(T, verbose=False)
    S = thermo.get_entropy(T, 1.0e5, verbose=False)

    free_energy = H-T*S
    return free_energy
	"""

	get_thermal_correction_ideal_gas = py"get_thermal_correction_ideal_gas"
	get_thermal_correction_adsorbate = py"get_thermal_correction_adsorbate"
end;

# ╔═╡ 0742ba27-af65-4c9d-b28d-795837185b1b
md"""
### Compute Rates
"""

# ╔═╡ 564ac0ed-ad22-4515-8a32-a149cb965499
md"""
### Explicit Calculation
"""

# ╔═╡ 3c8b1d70-c614-4733-98bc-704c2c827169
begin
	bulk_species        = ["H2_g", "CO2_g", "CO_g", "OH_g", "H2O_g", "ele_g"]
	surface_species     = ["CO2*_t", "COOH*_t", "CO*_t"]
	transition_state    = ["COOH-H2O-ele_t"]
	surface_site 		= ["*_t"]
end;

# ╔═╡ e59911c0-e234-459b-a6b0-804bee22cb1d
md"""
#### Formation Energies at 0 K from DFT calculations
"""

# ╔═╡ cb4fb272-1b8c-41b6-b264-881ec0e6e312
E_f = Dict( 
	"CO2_g" => 0.0 * eV,
	"CO2*_t" => 0.657600203 * eV,
	"COOH*_t" => 0.128214079 * eV,
	"COOH-H2O-ele_t" => 0.95 * eV,
	"CO*_t" => -0.02145440850567823 * eV,
	"CO_g" => 0.2701852315 * eV,
	"ele_g" => 0.0 * eV,
	"OH_g" => 0.0 * eV,
	"H2O_g" => 0.0 * eV,
	"H2_g" => 0.0 * eV,
	"*_t" => 0.0 * eV,
);

# ╔═╡ 25e1369f-a636-45a2-97fd-5a1f126dc489
md"""
#### Vibrational Frequencies from DFT Calculations
"""

# ╔═╡ 5814c5ff-ed3d-450a-be5a-c2f1d9dc742c
frequencies = Dict(
	        "CO2_g" => [24.1, 70.7, 635.8, 640.5, 1312.2, 2361.2] * ufac"h / cm / eV" * c_0,
	        "CO2*_t" => [136.85, 183.6, 212.95, 250.7, 306.0, 510.55, 562.25, 1176.05, 1889.85] * ufac"h / cm / eV" * c_0,
	        "COOH*_t" => [102.25, 172.05, 242.65, 265.55, 303.45, 520.8, 646.8, 794.1500000000001, 1032.6999999999998, 1344.35, 1658.15, 3551.35] * ufac"h  / cm / eV" * c_0,
	        "COOH-H2O-ele_t" => [],
	        "CO*_t" => [129.65, 155.55, 189.8, 227.5, 2073.25] * ufac"h / cm / eV" * c_0,
	        "CO_g" => [89.8, 127.2, 2145.5] * ufac"h / cm / eV" * c_0,
	        "ele_g" => [],
	        "OH_g" => [],
	        "H2O_g" => [103.0, 180.6, 245.1, 1625.9, 3722.8, 3830.3] * ufac"h / cm / eV" * c_0,
	        "H2_g" => [3.8, 13.5, 4444.5] * ufac"h / cm / eV" * c_0
);

# ╔═╡ a4acc450-bcf5-4938-85c2-f82492235890
md"""
#### Thermodynamical and Electrochemical Corrections

The energy calculations hold at 0 K and no excess surface charges. Therefore the energy terms have to be corrected according to the system parameters.
"""

# ╔═╡ 06d28618-e44c-4fd3-9940-3628d6816334
md"""
##### Correction of Bulk Species
"""

# ╔═╡ df663eb4-9e27-4ce8-9a93-3b8cd58870c9
ideal_gas_params = Dict(
	        "H2O_g" => (2, "nonlinear", 0, "H2O"),
	        "CO_g" => (1, "linear", 0, "CO"),
	        "CO2_g" => (2, "linear", 0, "CO2"),
	        "H2_g" => (2, "linear", 0, "H2")
);

# ╔═╡ ec5a8df2-c7a4-4665-982c-021251afdc3a
### ideal gases
function add_ideal_gas_correction!(thermo_corrections)
	for sp in bulk_species
		if sp ∉ ["OH_g", "ele_g"]
			thermo_corrections[sp] += get_thermal_correction_ideal_gas(T, frequencies[sp], ideal_gas_params[sp]...) * eV
		end
	end
	nothing
end

# ╔═╡ f98e4f84-1b0e-4bd9-a0e3-9b530f8f2afa
md"""
##### Correction of Adsorbates & BEB scaling of the Transition State
"""

# ╔═╡ 872c20f3-b9c1-40b7-8cc7-10eed29d6e16
### harmonic adsorbates
function add_harmonic_adsorbate_correction!(thermo_corrections)
	for sp ∈ surface_species
		thermo_corrections[sp] += get_thermal_correction_adsorbate(T, frequencies[sp]) * eV
		#println("$sp: $(E_f[sp] / eV)")
	end
	thermo_corrections["COOH-H2O-ele_t"] += (thermo_corrections["COOH*_t"] + thermo_corrections["CO*_t"]) / 2
	nothing
end

# ╔═╡ e900c3dd-874f-4d27-a01a-587861f69ab1
begin
	thermo_corrections = Dict(zip(keys(E_f), zeros(length(E_f))));
	add_ideal_gas_correction!(thermo_corrections)
	add_harmonic_adsorbate_correction!(thermo_corrections)
end;

# ╔═╡ 8fac9276-0fa1-4bb8-9cc8-55f9ede5994e
md"""
##### Electrochemical Corrections

The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$

The surface charging relation $σ = σ(U)$ is given by the Robin boundary condition

$σ(U) = C_{gap} (ϕ_{we} - ϕ_{pzc} - ϕ^\ddagger)$

where the gap capacitance between the working electrode and the reaction plane ($\ddagger$) is given by $C_{gap} = 20~μF/cm^2$. The potential of zero current is measured to be $ϕ_{pzc} = 0.16~V$.
"""

# ╔═╡ 4e0d3004-acaf-4098-b123-b378c03fbe99
begin
	C_gap = 20 * μF/cm^2
	ϕ_pzc = 0.16 * ufac"V"
end;

# ╔═╡ 0570d57c-5edb-41d2-a2a4-d0d2b4abe778
electro_correction_params = Dict(
	"CO2*_t"  => [-0.000286600929 / (μA/cm^2)^2, 0.0297720125 / (μA/cm^2)],
	"COOH*_t" => [-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)],
	"CO*_t"   => [-0.000189106972 / (μA/cm^2)^2,-0.00942574086 / (μA/cm^2)],
	"COOH-H2O-ele_t" => [-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)],
);

# ╔═╡ a0b39948-7364-4709-b0c2-f580797ac2aa
# ╠═╡ disabled = true
#=╠═╡
begin ## _get_echem_corrections
	G_H2O       = E_f["H2O_g"] + thermo_corrections["H2O_g"]
	G_H2        = E_f["H2_g"] + thermo_corrections["H2_g"]
	#G_H         = 0.5 * G_H2 - .0592 * local_pH / 298.14 * T * eV # wrong pH, use local pH
	G_OH        = G_H2O - G_H
	#thermo_corrections["OH_g"] += G_OH
end;
  ╠═╡ =#

# ╔═╡ f4413efe-7f2b-449a-8a91-91c3634e8737
md"""
#### Explicit Rate Calculations
"""

# ╔═╡ ddaba9e7-3af8-4bf8-82c5-298d37159040
function compute_rateconstants_explicit!(k, σ, ϕ_we, ϕ, local_pH)
	(kf, kr) = k
	β = 0.5 # hack
	
	Tval = eltype(kf)
	
	# compute corrections due to surface charge densities
	electro_corrections = Dict(zip(keys(E_f), zeros(Tval, length(E_f))))

	## simple_electrochemical
	#electro_corrections["ele_g"] -= (ϕ_we - (ϕ_we - u[iϕ])) * eV
	electro_corrections["ele_g"] -= (ϕ_we - ϕ) * eV
	#electro_corrections["COOH-H2O-ele_t"] += (-u[iϕ] + 0.5 * (u[iϕ] - ϕ_pzc)) * eV
	# Frumking correct: ϕ_we - u[iϕ]
	electro_corrections["COOH-H2O-ele_t"] += (-(ϕ_we - ϕ) + β * (ϕ_we - ϕ - ϕ_pzc)) * eV


	## hbond_electrochemical
	electro_corrections["COOH*_t"] += -0.25 * eV
	electro_corrections["CO2*_t"] += 0.0 * eV
	electro_corrections["CO*_t"] += -0.1 * eV

	## hbond_surface_charge_density
	for sp in [surface_species; transition_state]
		electro_corrections[sp] += electro_correction_params[sp]' * [σ^2, σ] * eV
	end

	## _get_echem_corrections
	G_H2O       = E_f["H2O_g"] + thermo_corrections["H2O_g"]
	G_H2        = E_f["H2_g"] + thermo_corrections["H2_g"]
	G_H         = 0.5 * G_H2 - .0592 * local_pH / 298.14 * T * eV
	G_OH        = G_H2O - G_H


	G_f = Dict(zip([bulk_species; surface_species; transition_state; surface_site], zeros(Tval, nc + na + length(transition_state) + length(surface_site))))
	
	for sp in [bulk_species; surface_species; transition_state; surface_site]
		G_f[sp] += E_f[sp] + thermo_corrections[sp] + electro_corrections[sp]
	end
	G_f["OH_g"] += G_OH

	# rate constants
	G_IS = zeros(Tval, 4)
	G_FS = zeros(Tval,4)
	G_TS = zeros(Tval, 4)
	
	# 'CO2_g + 2*_t <-> CO2*_t',	                  #1
	G_IS[1] = G_f["CO2_g"] + 2 * G_f["*_t"]
	G_FS[1] = G_f["CO2*_t"]
	G_TS[1] = max(G_IS[1], G_FS[1])

	kf[1] = 1.0e13 * exp(-(G_TS[1] - G_IS[1]) / (k_B * T))
	kr[1] = 1.0e13 * exp(-(G_TS[1] - G_FS[1]) / (k_B * T))

	# 'CO2*_t + H2O_g + ele_g <-> COOH*_t + OH_g',  #2            
	G_IS[2] = G_f["CO2*_t"] + G_f["H2O_g"] + G_f["ele_g"]
	G_FS[2] = G_f["COOH*_t"] + G_f["OH_g"]
	G_TS[2] = max(G_IS[2], G_FS[2])

	kf[2] = 1.0e13 * exp(-(G_TS[2] - G_IS[2]) / (k_B * T))
	kr[2] = 1.0e13 * exp(-(G_TS[2] - G_FS[2]) / (k_B * T))

	# 'COOH*_t + H2O_g + ele_g <-> COOH-H2O-ele_t <-> CO*_t + H2O_g + OH_g + *_t; beta=0.5', #3
	G_IS[3] = G_f["COOH*_t"] + G_f["H2O_g"] + G_f["ele_g"]
	G_FS[3] = G_f["CO*_t"] + G_f["H2O_g"] + G_f["OH_g"] + G_f["*_t"]
	G_TS[3] = G_f["COOH-H2O-ele_t"]

	kf[3] = 1.0e13 * exp(-(G_TS[3] - G_IS[3]) / (k_B * T))
	kr[3] = 1.0e13 * exp(-(G_TS[3] - G_FS[3]) / (k_B * T))

	# 'CO*_t <-> CO_g + *_t',	                      #4
	G_IS[4] = G_f["CO*_t"]
	G_FS[4] = G_f["CO_g"] + G_f["*_t"]
	G_TS[4] = max(G_IS[4], G_FS[4])

	kf[4] = 1.0e8 * exp(-(G_TS[4] - G_IS[4]) / (k_B * T))
	kr[4] = 1.0e8 * exp(-(G_TS[4] - G_FS[4]) / (k_B * T))
end;

# ╔═╡ ba1af1f7-8452-4288-900a-d6e305706d73
function compute_rates_explicit!(rates, u, σ, ϕ_we)
	
	Tval = eltype(u)
	
	kf = zeros(Tval, 4)
	kr = zeros(Tval, 4)
	
	local_pH = -log10(u[bulk["H⁺"].index] / (mol/dm^3))
	compute_rateconstants_explicit!((kf, kr), σ, ϕ_we, u[celldata.iϕ], local_pH)

	θ_free  = 1 - u[surface["CO₂"].index] - u[surface["CO"].index] - u[surface["COOH"].index]

	rates[1]  = kf[1] * (u[bulk["CO₂"].index] / Hcp_CO2 / bar) * θ_free^2 
	rates[1] -= kr[1] * u[surface["CO₂"].index]
	
	rates[2]  = kf[2] * u[surface["CO₂"].index] * 1.0 * 1.0 
	rates[2] -= kr[2] * u[surface["COOH"].index] * u[bulk["OH⁻"].index]
	
	rates[3]  = kf[3] * u[surface["COOH"].index] * 1.0 * 1.0 
	rates[3] -= kr[3] * u[surface["CO"].index] * 1.0 * u[bulk["OH⁻"].index] * θ_free
	
	rates[4]  = kf[4] * u[surface["CO"].index] 
	rates[4] -= kr[4] * (u[bulk["CO"].index] / Hcp_CO / bar) * θ_free

	nothing
end;

# ╔═╡ fdb9b695-6fcf-4535-989c-c2abdb9f69ab
md"""
### Catmap Interface
"""

# ╔═╡ efc51b73-609c-4d0e-a1ca-8224adb0f4b2
begin
	catmap_setup_file   = "catmap_CO2R_template.mkm"
    catmap_energy_file  = "catmap_CO2R_energies.txt"
end;

# ╔═╡ f55d15fd-ba9e-4dbb-9f6f-eba5c85f6139
mutable struct CatmapData
    rate_constant_fns
    compute_rates
end

# ╔═╡ a4191983-dedf-444d-9c2d-ac2bff89edd3
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
end;

# ╔═╡ f1a4de16-21d3-415a-b2ef-37436f04ad89
if rate_const_mode == "CatMAP"
  	instance_input_file_name = instantiate_catmap_template(catmap_setup_file, 0.0)
    
	(_, rate_constant_fns, _, compute_rates)  = get_catmap_output(
		instance_input_file_name, 
		catmap_energy_file, 
		true
	)
	catmap_data = CatmapData(rate_constant_fns, compute_rates)    
end;

# ╔═╡ 015d40f6-de19-4288-b5c3-5f85e0a8f745
md"""
#### CatMAP Rate Calculations
"""

# ╔═╡ 741d77a0-9376-4192-ba82-c9d7ae98001d
function compute_rates_catmap!(
	rates, 
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti},
	σ, 
	ϕ_we
) where {Tval, Tv, Tc, Tp, Ti}

	rates .= catmap_data.compute_rates(
		catmap_data.rate_constant_fns[1], 
		[u[surface["CO₂"].index], u[surface["COOH"].index], u[surface["CO"].index]], 
		[
			u[bulk["CO₂"].index] / Hcp_CO2 / bar, 
			u[bulk["CO"].index] / Hcp_CO / bar, 
			1.0, 
			1.0, 
			u[bulk["OH⁻"].index], 
			1.0
		],
		σ / (μF/cm^2) + 1.0e-8
	)[1:end-1]
	nothing
end;

# ╔═╡ 7d9a1531-c1d0-409b-9bf0-f361350bb7cb
function compute_rates!(
	rates, 
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	σ, 
	ϕ_we; mode = "Explicit"
) where {Tval, Tv, Tc, Tp, Ti}
	if mode == "Explicit"
		compute_rates_explicit!(rates, u, σ, ϕ_we)
	elseif mode == "CatMAP"
		compute_rates_catmap!(rates, u, σ, ϕ_we)
	else 
		throw(ArgumentError("The mode $mode for calculation of the rates is invalid."))
	end
end

# ╔═╡ 2fdaa88d-d04f-4dac-9362-b8996786d3a0
function we_breactions(
	f, 
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	bnode, 
	data
) where {Tval, Tv, Tc, Tp, Ti}
	
	(; iϕ, ϕ_we) = data


	σ = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)

	# rates of the elementary reactions
	rates 	 = zeros(Tval, 4)
	compute_rates!(rates, u, σ, ϕ_we; mode = rate_const_mode)

	# rates_cmp= zeros(Tval, 4)
	# compute_rates!(rates_cmp, u, σ, ϕ_we; mode = "Explicit")
	# @show ForwardDiff.value.(rates) .- ForwardDiff.value.(rates_cmp)

	#println("$(ForwardDiff.value.(rates))")
	
	# bulk species
	f[bulk["CO"].index]     += -rates[4] * S
	f[bulk["CO₂"].index]    += rates[1] * S
	f[bulk["OH⁻"].index] 	+= -rates[2] * S - rates[3] * S
	
	# surface species
	f[surface["CO₂"].index]  += -rates[1] + rates[2]
	f[surface["CO"].index]   += -rates[3] + rates[4]
	f[surface["COOH"].index] += -rates[2] + rates[3]
end;

# ╔═╡ f5231e9e-f00c-4b39-8fe7-b25f626e62d0
function halfcellbc(
	f,
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	bnode,
	data
) where {Tval, Tv, Tc, Tp, Ti}
	
	(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

	bulkbcondition(f, u, bnode, data; region = Γ_bulk)

	# Robin b.c. for the Poisson equation
	boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))

	if bnode.region == Γ_we
		if allow_surfacereactions
			we_breactions(f, u, bnode, data)
		end
	end
	nothing
end;

# ╔═╡ 0b3c0047-e8d2-49d6-9566-570fde1fad71
function pre(sol, t)
	
	instance_input_file_name = instantiate_catmap_template(catmap_setup_file, -t)
	
	(_, rate_constant_fns, _, compute_rates)  = get_catmap_output(
		instance_input_file_name, 
		catmap_energy_file, 
		true
	)
	
	catmap_data.rate_constant_fns           = rate_constant_fns
	catmap_data.compute_rates               = compute_rates

	nothing
end;

# ╔═╡ e77852ef-83de-44dc-8048-6bbea7aa7d9c
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:-0.0) * V, kwargs...)
    
    defaults = (; 	max_round 	= 3,
              		tol_round 	= 1.0e-9,
              		verbose 	= "e",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)
	
    kwargs 	= merge(defaults, kwargs) 

    
    cell        = PNPSystem(grid; bcondition=halfcellbc, reaction=reaction, celldata)
    
	more_pre 	= (rate_const_mode == "CatMAP") ? pre : ((args...) -> nothing)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, more_pre = more_pre, kwargs...)

	cell, ivresult
end;

# ╔═╡ bb1dcf91-578b-4440-a4b9-1617a4cd3032
(cell, ivresult) = simulate_CO2R(grid, celldata);

# ╔═╡ 0d945390-e2eb-427a-9879-016b77f85476
(_, default_index) = findmin(abs, ivresult.voltages .- -0.9 * V);

# ╔═╡ 7381aab1-e6a8-4e4f-9f82-b2114a31c011
md"""
Choose applied voltage: $(@bind ϕ_we_index PlutoUI.Slider(1:10:length(ivresult.voltages), default = default_index))
"""

# ╔═╡ ee06a76d-c191-4aa3-8e1f-75ce57513406
md"""
Potential at working electrode = $(ivresult.voltages[ϕ_we_index]) 
"""

# ╔═╡ 94e387ac-236e-46c6-9981-cbce1c973a1f
let
	F = N_A * e
	
	vis = GridVisualizer(; Plotter = PlutoVista, layout=(1,2), resolution = (1200, 300))
	
	# current-voltage plot
    currs = [j[bulk["OH⁻"].index] * F for j in ivresult.j_we]
	
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
	
	pHs = -log10.(ivresult.solutions[ϕ_we_index][bulk["H⁺"].index, :] / (mol/dm^3))
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

# ╔═╡ f4648583-9de7-47a2-8fae-ed84e8a687bd
md"""
### Save to csv-File
"""

# ╔═╡ a9992ee3-466c-44ce-a7dc-ee459729d57a
function save_ivresult(ivresult, mode)

	F = N_A * e
	currs = [j[bulk["OH⁻"].index] * F for j in ivresult.j_we] * ufac"cm^2/mA"
	
	open("data/CO2RCell_$(mode)_voltage-current-table.csv"; write=true) do file
		for (volt, curr) in zip(ivresult.voltages, currs)
			println(file, "$volt;$curr")
		end
	end

	open("data/CO2RCell_$(mode)_0point9volts.csv"; write=true) do file
		for (dist, vs) in zip(grid[Coordinates], eachcol(ivresult.solutions[ϕ_we_index]))
			println(file, "$dist,$(join(vs, ","))")
		end
	end
end

# ╔═╡ 755564b3-393c-47a5-917f-e52ffbdd01ec
save_ivresult(ivresult, rate_const_mode);

# ╔═╡ Cell order:
# ╟─92be2111-381e-4260-a514-1632c5a38355
# ╟─a532a10a-9cb8-44a9-89f4-f3517177b602
# ╠═99c3f409-dd51-4694-8f78-87a2aca41d40
# ╠═53c1639f-58e0-4b28-b4b6-34d7fd1f969f
# ╠═02aefeeb-4ebb-4daf-b06f-80d537409844
# ╟─306dc5b4-29f6-11ee-3112-874796d88919
# ╟─aaa10928-59d2-4651-b9d1-c71bb672ae95
# ╠═4d283cbb-6bd7-4f84-9e1b-9d5f0819aa4d
# ╟─ebbd9eb3-7dfc-417b-a39c-7023c15361c2
# ╟─57bdf775-9974-4f1f-8466-37f28b0c80b5
# ╠═43b24817-648d-4e0c-a885-f5445a726b5e
# ╠═58e784bc-f838-4ad2-87f5-c868313c0f91
# ╠═e3a2076f-f37b-4732-9238-49209253e3b8
# ╠═c0c3429a-a426-4623-aec7-977b7f9cb685
# ╠═4896bcdb-8661-4b8c-aaf3-5d2ae6eee1d1
# ╠═5391f6fc-ba0b-4fd1-bf4c-d26b38bd12dc
# ╟─e91db1f6-1388-4262-b0bc-235cae1867eb
# ╠═bb331dd8-3990-4cb5-87d4-ed4a93e54098
# ╠═f59ddfde-3eef-485c-b6a6-7bd050eb5b64
# ╠═4a8da214-c8b0-4b73-8f14-3010d1bd6756
# ╟─deffe3ec-6580-4222-9c9d-075e0a7cb958
# ╟─a9d434d9-7f7e-4c2f-95fd-3e24476e34ce
# ╠═a6e0b8fe-fa45-4c2e-a1c2-7c2710a65988
# ╟─c28be48d-1c3b-4581-bfa7-41c1ed61b51e
# ╟─ecc38283-877e-4a1b-9a0c-451c024c9683
# ╠═278cf5a7-da40-4106-b28e-b701359510d6
# ╟─df16bc86-8041-402c-8b94-35ef5aac961b
# ╠═2a7a3783-846d-4715-960e-979b72ded8d9
# ╟─04c5a531-f7bb-4a93-a619-2afd92586e2a
# ╟─faedb03c-c4ee-45cc-9263-0fd8eaf09850
# ╟─aed5431b-a89c-46c6-be7c-2d1c6fcd00c7
# ╠═2fdaa88d-d04f-4dac-9362-b8996786d3a0
# ╠═f5231e9e-f00c-4b39-8fe7-b25f626e62d0
# ╟─fb3df4df-2a66-4796-a7de-be46e8c8fd29
# ╠═e77852ef-83de-44dc-8048-6bbea7aa7d9c
# ╠═bb1dcf91-578b-4440-a4b9-1617a4cd3032
# ╟─11efa591-a50d-451b-b9bf-31c8d5c87a7c
# ╠═0d945390-e2eb-427a-9879-016b77f85476
# ╟─7381aab1-e6a8-4e4f-9f82-b2114a31c011
# ╟─ee06a76d-c191-4aa3-8e1f-75ce57513406
# ╠═94e387ac-236e-46c6-9981-cbce1c973a1f
# ╟─20114f4a-527f-4782-8b67-8bb8fe824319
# ╠═c64685fb-4c1c-435d-af95-7dafdfc77e15
# ╟─8905005c-0d88-4e2b-ae3e-cbee952d9a86
# ╠═4fda1782-4f70-4007-9ac3-12c931c4bd54
# ╟─0742ba27-af65-4c9d-b28d-795837185b1b
# ╠═7d9a1531-c1d0-409b-9bf0-f361350bb7cb
# ╟─564ac0ed-ad22-4515-8a32-a149cb965499
# ╠═3c8b1d70-c614-4733-98bc-704c2c827169
# ╟─e59911c0-e234-459b-a6b0-804bee22cb1d
# ╠═cb4fb272-1b8c-41b6-b264-881ec0e6e312
# ╟─25e1369f-a636-45a2-97fd-5a1f126dc489
# ╠═5814c5ff-ed3d-450a-be5a-c2f1d9dc742c
# ╟─a4acc450-bcf5-4938-85c2-f82492235890
# ╟─06d28618-e44c-4fd3-9940-3628d6816334
# ╠═df663eb4-9e27-4ce8-9a93-3b8cd58870c9
# ╠═ec5a8df2-c7a4-4665-982c-021251afdc3a
# ╟─f98e4f84-1b0e-4bd9-a0e3-9b530f8f2afa
# ╠═872c20f3-b9c1-40b7-8cc7-10eed29d6e16
# ╠═e900c3dd-874f-4d27-a01a-587861f69ab1
# ╟─8fac9276-0fa1-4bb8-9cc8-55f9ede5994e
# ╠═4e0d3004-acaf-4098-b123-b378c03fbe99
# ╠═0570d57c-5edb-41d2-a2a4-d0d2b4abe778
# ╠═a0b39948-7364-4709-b0c2-f580797ac2aa
# ╟─f4413efe-7f2b-449a-8a91-91c3634e8737
# ╠═ddaba9e7-3af8-4bf8-82c5-298d37159040
# ╠═ba1af1f7-8452-4288-900a-d6e305706d73
# ╟─fdb9b695-6fcf-4535-989c-c2abdb9f69ab
# ╠═efc51b73-609c-4d0e-a1ca-8224adb0f4b2
# ╠═f55d15fd-ba9e-4dbb-9f6f-eba5c85f6139
# ╠═a4191983-dedf-444d-9c2d-ac2bff89edd3
# ╠═f1a4de16-21d3-415a-b2ef-37436f04ad89
# ╟─015d40f6-de19-4288-b5c3-5f85e0a8f745
# ╠═741d77a0-9376-4192-ba82-c9d7ae98001d
# ╠═0b3c0047-e8d2-49d6-9566-570fde1fad71
# ╟─f4648583-9de7-47a2-8fae-ed84e8a687bd
# ╠═a9992ee3-466c-44ce-a7dc-ee459729d57a
# ╠═755564b3-393c-47a5-917f-e52ffbdd01ec
