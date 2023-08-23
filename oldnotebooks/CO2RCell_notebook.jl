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

# ╔═╡ e82ccb4e-2223-11ee-1325-ffa068db7408
begin
	using Pkg
	Pkg.activate(@__DIR__)
end;

# ╔═╡ e03a851c-e822-4318-bcd2-3cef1d5e6c2d
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
end;

# ╔═╡ 0da68bd3-6821-4cc2-8591-8f89f47337f6
begin
	ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
	ENV["PYTHON"] = Sys.which("python3")
	using PyCall
	#Pkg.build("PyCall")
end;

# ╔═╡ 4b4dc09e-ad77-4a30-801e-53ae8bdfbc2a
md"""
## Some Python code
"""

# ╔═╡ 6b72af93-cb91-40fe-9224-3ebb90087d75
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

# ╔═╡ 7e845d0f-f498-486e-ab68-7f40c077473f
@unitfactors mol dm eV μA μF cm μm bar;

# ╔═╡ 40a14821-8114-4fd6-8dac-34553525e466
begin
	@phconstants c_0 N_A e k_B
	F = N_A * e
end;

# ╔═╡ a08ee2a2-c160-42a9-ac8b-9c9c458a82dc
md"""
## System Parameters
"""

# ╔═╡ 9b22386b-840c-4fed-bf74-d3f9805efcbb
pH = 6.8;

# ╔═╡ de15dd2c-6286-4a63-9866-7d0e013bbb3b
T = 273.15 + 25 * ufac"K";

# ╔═╡ 689756df-805b-418d-9de4-030aebf00321
md"""
#### Diffusion constants
"""

# ╔═╡ ef681f56-bed0-4c7a-804c-3f0382aceb69
D = [1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s"; # from Ringe paper

# ╔═╡ 90d9c558-48d6-4568-8f96-5b64417217e3
md"""
## Kinetic Model for the Buffer Reactions
"""

# ╔═╡ d6133698-9ff9-44d4-ba23-05ff0e9ebd60
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

# ╔═╡ 4403fb73-7d3d-4cf3-b202-2f7539b86806
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

# ╔═╡ cacb012a-4847-46d9-94e5-42480038ba76
md"""
### Reaction Constants
"""

# ╔═╡ 9b98e9db-9856-46d4-8e7e-7492f2e43df5
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

# ╔═╡ 73cd7231-5847-41c0-95e3-290127545039
md"""
### Kinetic Model
"""

# ╔═╡ 153e10f8-ec31-401e-a5cd-21a7f4edcedf
md"""
For now the mean-field approach using the law of mass action is applied.
"""

# ╔═╡ 85c6dc01-7d7c-4f49-b4fa-0eb02bb958ad
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

# ╔═╡ 3085f165-6228-41dc-aa6b-9c1b76207c35
md"""
## Microkinetic Model for the Surface Reactions
"""

# ╔═╡ 757bc30e-d26e-474c-ad17-28cb864c83ed
md"""
Enable surface reactions: $(@bind allow_surfacereactions PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 38a29731-8e07-408c-90ec-1097765c81ab
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

# ╔═╡ 251a0e47-13e3-4052-8a00-56683b97df00
begin # surface species
	ico_ad      = 8
    ico2_ad     = 9
    icooh_ad    = 10
	na 			= 3
end;

# ╔═╡ c1c3b140-200f-4e28-8b9e-e5328e6954f4
begin
	bulk_species        = ["H2_g", "CO2_g", "CO_g", "OH_g", "H2O_g", "ele_g"]
	surface_species     = ["CO2*_t", "COOH*_t", "CO*_t"]
	transition_state    = ["COOH-H2O-ele_t"]
end;

# ╔═╡ 496ff70a-2606-4f99-8b42-32356ea66c68
md"""
#### Formation Energies at 0 K from DFT calculations
"""

# ╔═╡ 646ec6d2-29e4-44af-834c-57fc520e0fb6
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
	        "H2_g" => 0.0 * eV
	    )

# ╔═╡ af12d451-0209-42d9-a713-bcba05deb821
md"""
#### Vibrational Frequencies from DFT Calculations
"""

# ╔═╡ 2171c678-5ce9-47ba-b249-1eeb69d3b574
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
	    )

# ╔═╡ dcb2afb7-dadb-4f5b-8161-ed8e4a9b0b59
md"""
### Thermodynamical and Electrochemical Corrections

The energy calculations hold at 0 K and no excess surface charges. Therefore the energy terms have to be corrected according to the system parameters.
"""

# ╔═╡ c35ce560-b61e-4005-b35d-3ed8e675c972
thermo_corrections = Dict(zip(keys(E_f), zeros(length(E_f))));

# ╔═╡ 78ad5cf3-0776-46d5-b006-e66958f249b8
md"""
#### Correction of Bulk species
"""

# ╔═╡ 2fd46ff3-37ea-42cd-9344-21e3d82d927d
ideal_gas_params = Dict(
	        "H2O_g" => (2, "nonlinear", 0, "H2O"),
	        "CO_g" => (1, "linear", 0, "CO"),
	        "CO2_g" => (2, "linear", 0, "CO2"),
	        "H2_g" => (2, "linear", 0, "H2")
);

# ╔═╡ b6f97a02-f59b-4f69-a11a-d111eab8f203
### ideal gases
for sp in bulk_species
	if sp ∉ ["OH_g", "ele_g"]
		thermo_corrections[sp] += py"get_thermal_correction_ideal_gas"(T, frequencies[sp], ideal_gas_params[sp]...) * eV
	end
end

# ╔═╡ e4bafffd-5de0-4495-a0b3-f76e3d6eef33
md"""
#### Correction of Adsorbates
"""

# ╔═╡ e5a6c4b1-f2f8-474b-80e2-22ddae581e26
### harmonic adsorbates
for sp ∈ surface_species
	thermo_corrections[sp] += py"get_thermal_correction_adsorbate"(T, frequencies[sp]) * eV
	#println("$sp: $(E_f[sp] / eV)")
end;

# ╔═╡ be46a3bf-2dd1-405a-9829-f3a762edcf90
md"""
#### BEB scaling of the Transition State
"""

# ╔═╡ c4f2394e-39ef-4778-a876-f02f9a687d4f
thermo_corrections["COOH-H2O-ele_t"] += (thermo_corrections["COOH*_t"] + thermo_corrections["CO*_t"]) / 2;

# ╔═╡ 82a698df-15e9-4e75-93cf-b51a79eee8ca
md"""
#### Electrochemical Corrections

The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$

The surface charging relation $σ = σ(U)$ is given by the Robin boundary condition

$σ(U) = C_{gap} (ϕ_{we} - ϕ_{pzc} - ϕ^\ddagger)$

where the gap capacitance between the working electrode and the reaction plane ($\ddagger$) is given by $C_{gap} = 20~μF/cm^2$. The potential of zero current is measured to be $ϕ_{pzc} = 0.16~V$.
"""

# ╔═╡ 33be6880-3213-46d4-ae14-7aafa7f096d7
begin
	C_gap = 20 * μF/cm^2
	ϕ_pzc = 0.16 * ufac"V"
end;

# ╔═╡ 84ecc576-6d69-44e9-b70b-34773c765a01
electro_correction_params = Dict(
	"CO2*_t"  => [-0.000286600929 / (μA/cm^2)^2, 0.0297720125 / (μA/cm^2)],
	"COOH*_t" => [-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)],
	"CO*_t"   => [-0.000189106972 / (μA/cm^2)^2,-0.00942574086 / (μA/cm^2)]
);

# ╔═╡ 8395fadf-9f18-40ac-8b8e-ae577ff6e96a


# ╔═╡ bc74780a-69d9-4433-aa4a-f57c40e4e366
begin ## _get_echem_corrections
	G_H2O       = E_f["H2O_g"] + thermo_corrections["H2O_g"]
	G_H2        = E_f["H2_g"] + thermo_corrections["H2_g"]
	G_H         = 0.5 * G_H2 - .0592 * pH / 298.14 * T * eV
	G_OH        = G_H2O - G_H
	thermo_corrections["OH_g"] += G_OH
end;

# ╔═╡ 19e7678e-378a-4334-821e-6f24a9e478b5
md"""
Henry's Law
"""

# ╔═╡ 393d97a5-d2f1-4725-a831-701426228022
begin
	# Henry constants
    Hcp_CO  = 9.7e-6 * ufac"mol/(m^3 * Pa)"
    Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"
end;

# ╔═╡ c7e01f10-09b2-4fac-901e-c03ba440abb0
function compute_rateconstants!(k, σ, ϕ_we)
	(kf, kr) = k
	
	Tval = eltype(kf)
	
	# compute corrections due to surface charge densities
	electro_corrections = Dict(zip(keys(E_f), zeros(Tval, length(E_f))))

	## simple_electrochemical
	#electro_corrections["ele_g"] -= (ϕ_we - (ϕ_we - u[iϕ])) * eV
	electro_corrections["ele_g"] -= ϕ_we * eV
	#electro_corrections["COOH-H2O-ele_t"] += (-u[iϕ] + 0.5 * (u[iϕ] - ϕ_pzc)) * eV
	electro_corrections["COOH-H2O-ele_t"] += (-ϕ_we + 0.5 * ϕ_we) * eV

	## hbond_electrochemical
	electro_corrections["COOH*_t"] += -0.25 * eV
	electro_corrections["CO2*_t"] += 0.0 * eV
	electro_corrections["CO*_t"] += -0.1 * eV

	## hbond_surface_charge_density
	for sp in surface_species
		electro_corrections[sp] += electro_correction_params[sp]' * [σ^2, σ] * eV
	end


	G_f = Dict(zip([bulk_species; surface_species; transition_state], zeros(Tval, nc + na + length(transition_state))))
	
	for sp in [bulk_species; surface_species; transition_state]
		G_f[sp] += E_f[sp] + thermo_corrections[sp] + electro_corrections[sp]
	end

	# rate constants
	G_IS = zeros(Tval, 4)
	G_FS = zeros(Tval,4)
	G_TS = zeros(Tval, 4)


	# 'CO2_g + 2*_t <-> CO2*_t',	                  #1
	G_IS[1] = G_f["CO2_g"] + 2 * 0.0
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
	G_FS[3] = G_f["CO*_t"] + G_f["H2O_g"] + G_f["OH_g"] + 0.0
	G_TS[3] = G_f["COOH-H2O-ele_t"]

	kf[3] = 1.0e13 * exp(-(G_TS[3] - G_IS[3]) / (k_B * T))
	kr[3] = 1.0e13 * exp(-(G_TS[3] - G_FS[3]) / (k_B * T))

	# 'CO*_t <-> CO_g + *_t',	                      #4
	G_IS[4] = G_f["CO*_t"]
	G_FS[4] = G_f["CO_g"] + 0.0
	G_TS[4] = max(G_IS[4], G_FS[4])

	kf[4] = 1.0e8 * exp(-(G_TS[4] - G_IS[4]) / (k_B * T))
	kr[4] = 1.0e8 * exp(-(G_TS[4] - G_FS[4]) / (k_B * T))
end;

# ╔═╡ eefa6514-77bd-4759-a986-5a0c132e1f6a
function we_breactions(f, u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode, data) where {Tval, Tv, Tc, Tp, Ti}
	(; iϕ, ϕ_we) = data

	
	#Tval = eltype(u)
	kf = zeros(Tval, 4)
	kr = zeros(Tval, 4)
	σ = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)

	compute_rateconstants!((kf, kr), σ, ϕ_we)
			
	S       = 9.61e-5 / N_A * (1.0e10)^2 * ufac"mol/m^2"
	θ_free  = (1- u[ico2_ad] - u[ico_ad] - u[icooh_ad])

	# rates of the elementary reactions
	rates 	 = zeros(Tval, 4)
	rates[1] = kf[1] * (u[ico2] / Hcp_CO2 / bar) * θ_free^2 - kr[1] * u[ico2_ad]
	rates[2] = kf[2] * u[ico2_ad] * 1.0 * 1.0 - kr[2] * u[icooh_ad] * u[iohminus]
	rates[3] = kf[3] * u[icooh_ad] * 1.0 * 1.0 - kr[3] * u[ico_ad] * 1.0 * u[iohminus] * θ_free
	rates[4] = kf[4] * u[ico_ad] - kr[4] * (u[ico] / Hcp_CO / bar) * θ_free
	
	#println("rate constants: $(ForwardDiff.value.(kf)) and $(ForwardDiff.value.(kr))")
	#println("rates: $(ForwardDiff.value.(rates)): $(eltype(rates))")
	
	# bulk species
	f[ico]      += -rates[4] * S
	f[ico2]     += rates[1] * S
	f[iohminus] += -rates[2] * S - rates[3] * S
	
	# surface species
	f[ico2_ad]  += -rates[1] + rates[2]
	f[ico_ad]   += -rates[3] + rates[4]
	f[icooh_ad] += -rates[2] + rates[3]
end;

# ╔═╡ 380b4c07-1c4e-4cea-8754-75b785952697
function halfcellbc(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
	(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

	bulkbcondition(f, u, bnode, data; region = Γ_bulk)
	#boundary_dirichlet!(f, u, bnode; species=iϕ, region=Γ_we, value=ϕ_we*0.18)

	# Robin b.c. for the Poisson equation
	boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))

	if bnode.region == Γ_we

		if allow_surfacereactions
			we_breactions(f, u, bnode, data)
		end
	end
	nothing
end;

# ╔═╡ 1944404c-57ee-413c-8c44-56ca11032ddb
md"""
## Simulation of the $CO_2$ reduction
"""

# ╔═╡ 6b919937-582f-40e9-ad86-cdae31486547
function simulate_CO2R(; nref 		= 0,
              			 voltages 	= (-1.5:0.1:-0.0) *ufac"V",
              			 scheme 	= :μex,
              			 xmax 		= 1,
               			 κ 			= 4.0,
              			 kwargs...)
    
    defaults = (; 	max_round 	= 3,
              		tol_round 	= 1.0e-9,
              		verbose 	= "e",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)
    kwargs 	= merge(defaults, kwargs) 

    hmin 	= 1.0e-6 * μm * 2.0^(-nref)
    hmax 	= 1.0 * μm * 2.0^(-nref)
    L 		= 80.0 * μm
    X 		= geomspace(0, L, hmin, hmax)
    grid 	= simplexgrid(X)

    
    celldata=ElectrolyteData(;nc    = nc,
                              na    = allow_surfacereactions ? na : 0,
                              z     = [1, -1, -2, 0, 0, -1, 1],
                              D     = D,
                              T     = T,
                              eneutral=false,
                              κ     = fill(κ,7),
                              Γ_we  = 1,
                              Γ_bulk= 2,
                              scheme)

    (;iϕ::Int,ip::Int)=celldata
    
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
end;

# ╔═╡ 474f5a0f-a574-4b28-b5ad-71db62fc99d6
(cell, ivresult) = simulate_CO2R(; κ = 4.0);

# ╔═╡ ba055678-22cf-4ff5-a032-70c3c55686db
md"""
## Visualization
"""

# ╔═╡ f4ef9513-b18e-48cb-bd86-1dddbad8b561
md"""
Choose applied voltage: $(@bind ϕ_we_index PlutoUI.Slider(1:10:length(ivresult.voltages)))
"""

# ╔═╡ be6e227f-0cd4-4aa8-83ad-802bf2e3f7e6
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

# ╔═╡ 393bb858-2bd1-4731-ab7a-c6e480bca75e
md"""
Potential at working electrode = $(ivresult.voltages[ϕ_we_index]) 
"""

# ╔═╡ 4ea31264-8724-400a-b55a-fe3b31add833
open("CO2RCell.csv"; write=true) do file
	for (volt, curr) in zip(ivresult.voltages, currs)
		println(file, "$volt;$curr")
	end
end

# ╔═╡ e98393b9-1b48-47b9-8632-b251b3f234b7
θ_free = 1 - ivresult.solutions[ϕ_we_index][ico2_ad, 1] - ivresult.solutions[ϕ_we_index][ico_ad, 1] - ivresult.solutions[ϕ_we_index][icooh_ad, 1]

# ╔═╡ be056fec-dbaf-4b1a-8fd3-ba5b18b12584
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ Cell order:
# ╠═e82ccb4e-2223-11ee-1325-ffa068db7408
# ╠═e03a851c-e822-4318-bcd2-3cef1d5e6c2d
# ╟─4b4dc09e-ad77-4a30-801e-53ae8bdfbc2a
# ╠═0da68bd3-6821-4cc2-8591-8f89f47337f6
# ╠═6b72af93-cb91-40fe-9224-3ebb90087d75
# ╠═7e845d0f-f498-486e-ab68-7f40c077473f
# ╠═40a14821-8114-4fd6-8dac-34553525e466
# ╟─a08ee2a2-c160-42a9-ac8b-9c9c458a82dc
# ╠═9b22386b-840c-4fed-bf74-d3f9805efcbb
# ╠═de15dd2c-6286-4a63-9866-7d0e013bbb3b
# ╟─689756df-805b-418d-9de4-030aebf00321
# ╠═ef681f56-bed0-4c7a-804c-3f0382aceb69
# ╟─90d9c558-48d6-4568-8f96-5b64417217e3
# ╟─d6133698-9ff9-44d4-ba23-05ff0e9ebd60
# ╠═4403fb73-7d3d-4cf3-b202-2f7539b86806
# ╟─cacb012a-4847-46d9-94e5-42480038ba76
# ╠═9b98e9db-9856-46d4-8e7e-7492f2e43df5
# ╟─73cd7231-5847-41c0-95e3-290127545039
# ╟─153e10f8-ec31-401e-a5cd-21a7f4edcedf
# ╠═85c6dc01-7d7c-4f49-b4fa-0eb02bb958ad
# ╟─3085f165-6228-41dc-aa6b-9c1b76207c35
# ╟─757bc30e-d26e-474c-ad17-28cb864c83ed
# ╟─38a29731-8e07-408c-90ec-1097765c81ab
# ╠═251a0e47-13e3-4052-8a00-56683b97df00
# ╠═c1c3b140-200f-4e28-8b9e-e5328e6954f4
# ╟─496ff70a-2606-4f99-8b42-32356ea66c68
# ╠═646ec6d2-29e4-44af-834c-57fc520e0fb6
# ╟─af12d451-0209-42d9-a713-bcba05deb821
# ╠═2171c678-5ce9-47ba-b249-1eeb69d3b574
# ╟─dcb2afb7-dadb-4f5b-8161-ed8e4a9b0b59
# ╠═c35ce560-b61e-4005-b35d-3ed8e675c972
# ╟─78ad5cf3-0776-46d5-b006-e66958f249b8
# ╠═2fd46ff3-37ea-42cd-9344-21e3d82d927d
# ╠═b6f97a02-f59b-4f69-a11a-d111eab8f203
# ╟─e4bafffd-5de0-4495-a0b3-f76e3d6eef33
# ╠═e5a6c4b1-f2f8-474b-80e2-22ddae581e26
# ╟─be46a3bf-2dd1-405a-9829-f3a762edcf90
# ╠═c4f2394e-39ef-4778-a876-f02f9a687d4f
# ╟─82a698df-15e9-4e75-93cf-b51a79eee8ca
# ╠═33be6880-3213-46d4-ae14-7aafa7f096d7
# ╠═84ecc576-6d69-44e9-b70b-34773c765a01
# ╠═8395fadf-9f18-40ac-8b8e-ae577ff6e96a
# ╠═bc74780a-69d9-4433-aa4a-f57c40e4e366
# ╟─19e7678e-378a-4334-821e-6f24a9e478b5
# ╠═393d97a5-d2f1-4725-a831-701426228022
# ╠═c7e01f10-09b2-4fac-901e-c03ba440abb0
# ╠═eefa6514-77bd-4759-a986-5a0c132e1f6a
# ╠═380b4c07-1c4e-4cea-8754-75b785952697
# ╟─1944404c-57ee-413c-8c44-56ca11032ddb
# ╠═6b919937-582f-40e9-ad86-cdae31486547
# ╠═474f5a0f-a574-4b28-b5ad-71db62fc99d6
# ╟─ba055678-22cf-4ff5-a032-70c3c55686db
# ╠═be6e227f-0cd4-4aa8-83ad-802bf2e3f7e6
# ╟─f4ef9513-b18e-48cb-bd86-1dddbad8b561
# ╟─393bb858-2bd1-4731-ab7a-c6e480bca75e
# ╠═4ea31264-8724-400a-b55a-fe3b31add833
# ╠═e98393b9-1b48-47b9-8632-b251b3f234b7
# ╟─be056fec-dbaf-4b1a-8fd3-ba5b18b12584
