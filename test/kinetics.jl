### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 20f09fe8-3d2f-11ee-3749-550c6dbe9f0e
begin
	import Pkg as _Pkg
	_Pkg.activate(joinpath(@__DIR__, ".."))
	using PyCall
	using DataFrames
	using LessUnitful
end

# ╔═╡ 44a506d1-b0e8-4bb0-8265-fac97271efff
md"""
## Units
"""

# ╔═╡ 3afa3da7-d14b-43a9-bd43-c0ce0550086f
begin
	@unitfactors cm eV μF μA
	@phconstants c_0 h
end;

# ╔═╡ 1dc9ffb6-02e2-4278-99c8-d1495818256a
md"""
## Data
"""

# ╔═╡ b24253cc-e044-48ad-8d1c-f12a6765c398
begin
	const ϕ_we 	= -0.8
	const T 	= 298.0
	const C_gap = 20 * μF/cm^2
	const ϕ_pzc = 0.16 * ufac"V"
	const local_pH = 6.8
end;

# ╔═╡ 15b52261-ba3e-42ba-96a9-c772626e5496
const gases = DataFrame(
	:name => ["H₂O", "CO", "CO₂", "H₂"],
	:components => ["H2O", "CO", "CO2", "H2"],
	:symmetrynumber => [2, 1, 2, 2],
	:geometry => ["nonlinear", "linear", "linear", "linear"],
	:spin => [0, 0, 0, 0],
	:energy => [0.0, 0.2701852315, 0.0, 0.0] .* eV,
	:frequencies => [
		[103.0, 180.6, 245.1, 1625.9, 3722.8, 3830.3] * h / cm / eV * c_0,
		[89.8, 127.2, 2145.5] * h / cm / eV * c_0,
		[24.1, 70.7, 635.8, 640.5, 1312.2, 2361.2] * h / cm / eV * c_0,
		[3.8, 13.5, 4444.5] * h / cm / eV * c_0
	],

)

# ╔═╡ b66f7fe9-009d-4287-9d89-a36dbde1e434
const adsorbates = DataFrame(
	:name 		 => ["CO₂_t" 	, "COOH_t" 	 , "CO_t"],
	:energy 	 => [0.657600203, 0.128214079, -0.02145440850567823] .* eV,
	:frequencies => [
		[136.85, 183.6, 212.95, 250.7, 306.0, 510.55, 562.25, 1176.05, 1889.85] * h / cm / eV * c_0,
		[102.25, 172.05, 242.65, 265.55, 303.45, 520.8, 646.8, 794.1500000000001, 1032.6999999999998, 1344.35, 1658.15, 3551.35] * h / cm / eV * c_0,
		[129.65, 155.55, 189.8, 227.5, 2073.25] * h / cm / eV * c_0
	],
	:echem_params => [
		[-0.000286600929 / (μA/cm^2)^2, 0.0297720125 / (μA/cm^2)],
		[-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)],
		[-0.000189106972 / (μA/cm^2)^2,-0.00942574086 / (μA/cm^2)],
	],
	
)

# ╔═╡ ea5c99f8-cbd5-4542-ab17-f2582b49e565
const transition_states = DataFrame(
	:name => ["COOH-H₂O-e⁻_t"],
	:energy => [0.95] .* eV,
	:frequencies => [Float64[]],
	:echem_params => [[-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)]],
)

# ╔═╡ 29b9eacb-ef85-400d-8301-2e59d541b70f
const surfaces = DataFrame(
	:name => ["_t"],
	:energy => [0.0] .* eV,
)

# ╔═╡ 7787f977-66fb-4f0f-8ef2-41c78a518933
const fictious = DataFrame(
	:name => ["OH⁻", "e⁻"], 
	:energy => [0.0, 0.0], 
	:frequencies => [Float64[], Float64[]]
)

# ╔═╡ ce88f49d-2ea6-40fd-8e21-9fcecc938e02
const thermo_corrections = let
	ngases 		= nrow(gases)
	nadsorbates = nrow(adsorbates)
	ntransition = nrow(transition_states)
	nsurfaces 	= nrow(surfaces)
	nfictious 	= nrow(fictious)
	ntotal 		= ngases + nadsorbates + ntransition + nsurfaces + nfictious
	Dict(zip(
		[gases.name; adsorbates.name; transition_states.name; surfaces.name; fictious.name], zeros(ntotal)
	))
end

# ╔═╡ 7d4c60a7-8a5b-47f1-a91e-bb1b2c76c05d
const energies = let
	dicts = []
	for df in [gases, adsorbates, transition_states, fictious]
		push!(dicts, Dict(name => energy for (; name, energy) in eachrow(df)))
	end
	merge(dicts...)
end

# ╔═╡ 26cf49cc-7a4e-4232-85b2-b370db46378e
const translation = Dict(
		"H2_g" => "H₂",
		"H2O_g" => "H₂O",
		"CO2_g" => "CO₂",
		"CO_g" => "CO",
		"CO2_t" => "CO₂_t",
		"COOH_t" => "COOH_t",
		"CO_t" => "CO_t",
		"COOH-H2O-ele_t" => "COOH-H₂O-e⁻_t",
		"OH_g" => "OH⁻",
		"ele_g" => "e⁻",
		"t" => "_t",
)

# ╔═╡ ba894fa1-fba8-40db-9459-7ec78462a246
md"""
## CatMAP
"""

# ╔═╡ 1f020e47-e4cb-4202-80f7-39421e4ce471
function instantiate_catmap_template(template_file_path, ϕ_we, ϕ, local_pH)
    
    input_instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

	replacements = Dict(
		r"descriptor_ranges.?=.*" =>"descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]",
		r"voltage_diff_drop.?=.*" => "voltage_diff_drop = $ϕ",
		r"pH.?=.*" => "pH = $local_pH",
		r"Upzc.?=.*" => "Upzc = $ϕ_pzc",
		r"\nsigma_input.?=.*" => "\nsigma_input = $((ϕ_we - ϕ - ϕ_pzc) * C_gap/(μF/cm^2))",
	)
    input_instance_string = replace(input_instance_string, replacements...)

    input_instance_file_name = "CO2R.mkm"
    open(input_instance_file_name, "w") do input_instance_file
        write(input_instance_file, input_instance_string)
    end

    return input_instance_file_name
end

# ╔═╡ 9d4a2505-9c24-4ed2-bf02-f6598f296f11
begin
	template_file_path = abspath("../catmap_CO2R_template.mkm")
	isfile(template_file_path)
end

# ╔═╡ 6120ec04-80c8-423c-a92a-44583ba9bac7
catmap_df = let
py"""
from catmap import ReactionModel

def catmap_kinetic_model(setup_file):

	# ReactionModel is the main class that is initialized with a setup-file
	model = ReactionModel(setup_file=setup_file)

	# some solver parameters have to be set manually (?!)
	import mpmath as mp
	model.solver._mpfloat = mp.mpf
	model.solver._math = mp
	model.solver._matrix = mp.matrix

	model.solver.compile() # compiles all templates, here (rate_constants) are needed

	# define descriptor grid
	## model.resolution can either be a list of integers or an integer, if it's only
	## an integer, then duplicate the same resolution for all descriptors to normalize
	## model.resolution to a list
	if isinstance(model.resolution, int):
		model.resolution = [model.resolution] * len(model.descriptor_names)

	ndescriptor_tuples = 1
	for res in model.resolution:
		ndescriptor_tuples *= res

	## Note: define linspace to keep dependencies minimal
	def linspace(start, stop, nitems):
		if start >= stop:
			return [start]
		delta = (stop - start)/(nitems-1)
		r = []
		for i in range(nitems):
			r.append(start + i * delta)
		return r

	## list of the ranges of the descriptors
	descriptor_ranges = [linspace(_range[0], _range[1], res) for _range, res in zip(model.descriptor_ranges, model.resolution)]

	## The combinations of the descriptors are ordered in lexicographical order
	## This function returns the position of the combination at position idx in the
	## descriptor grid
	def get_idx_tuple(idx, lengths):
		idx_tuple = []
		rest = idx
		for length in lengths:
			idx_tuple.append(rest % length)
			rest = rest // length
		return idx_tuple
	
	descriptor_grid = [0] * ndescriptor_tuples
	rxn_params_grid = [0] * ndescriptor_tuples
	
	for i in range(ndescriptor_tuples):		
		descriptor_values = [
			descriptor_range[idx] for idx, descriptor_range in zip(get_idx_tuple(i, model.resolution), descriptor_ranges)
		]
		descriptor_grid[i] = descriptor_values
		rxn_params_grid[i] = {
			k:v for k, v in zip(
				model.adsorbate_names + model.transition_state_names,
				model.scaler.get_rxn_parameters(descriptor_values)
			)
		}
		for k, v in zip(model.gas_names, model.solver._gas_energies):
			rxn_params_grid[i][k] = v
		for k, v in zip(model.site_names, model.solver._site_energies):
			rxn_params_grid[i][k] = v
	return rxn_params_grid[0], model.species_definitions
"""
	df = nothing
	for ϕ in range(ϕ_we, 0.0, 11)
		σ = (ϕ_we - ϕ - ϕ_pzc) * C_gap
		params = Dict("ϕ_we" => ϕ_we, "ϕ" => ϕ, "σ" => σ)
		
		catmap_instance_path = instantiate_catmap_template(template_file_path, ϕ_we, ϕ, local_pH)
		catmap_instance_path = abspath(joinpath(".", catmap_instance_path))
		energies, specd = py"catmap_kinetic_model"(catmap_instance_path)
		energies = convert(Dict{String, Float64}, energies)
		
		if isnothing(df)
			df = DataFrame(params..., energies...)
		else
			push!(df, merge(params, energies))
		end
	end
	rename!(df, translation)
	df
end

# ╔═╡ 8a87380d-20e6-4c27-950b-25003b98bf30
md"""
# Explicit
"""

# ╔═╡ 30c7d656-e197-495e-bd71-9c71214d1857
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
	get_thermal_correction_adsorbate = py"get_thermal_correction_adsorbate"
	get_thermal_correction_ideal_gas = py"get_thermal_correction_ideal_gas"
end

# ╔═╡ 0b3c9e82-c87b-4f3c-b504-af10d355c0e2
function add_simple_electrochem_corrections!(electro_corrections, ϕ_we, ϕ)
	β = 0.50 # hack
	electro_corrections["e⁻"] -= (ϕ_we - ϕ) * eV
	# Frumking correct: ϕ_we - u[iϕ]
	electro_corrections["COOH-H₂O-e⁻_t"] += (-(ϕ_we - ϕ) + β * (ϕ_we - ϕ - ϕ_pzc)) * eV
	nothing
end

# ╔═╡ dc768d05-de5a-4c44-b30a-7aea60fb1926
function add_hbond_corrections!(electro_corrections)
	## hbond_electrochemical
	electro_corrections["COOH_t"] += -0.25 * eV
	electro_corrections["CO₂_t"] += 0.0 * eV
	electro_corrections["CO_t"] += -0.1 * eV
end

# ╔═╡ 119b4542-41e6-49fc-a281-a0c2c1b5981f
md"""
The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$
"""

# ╔═╡ 36a30e9a-6fce-40d3-a36a-aaef3d8b7f2a
function add_hbond_surface_charge_corrections!(electro_corrections, σ)
	## hbond_surface_charge_density
	for species in [adsorbates, transition_states]
		for (; name, echem_params) in eachrow(species)
			electro_corrections[name] += echem_params' * [σ^2, σ] * eV
		end
	end
end

# ╔═╡ 5ccf3fa5-df4d-47e5-a56d-82ba4b099b92
function add_pH_correction!(electro_corrections, thermo_corrections, local_pH)
	## _get_echem_corrections
	gases_by_name = groupby(gases, :name)
	G_H2O       = energies["H₂O"] + thermo_corrections["H₂O"]
	G_H2        = energies["H₂"] + thermo_corrections["H₂"]
	G_H         = 0.5 * G_H2 - .0592 * local_pH / 298.14 * T * eV
	G_OH        = G_H2O - G_H
	electro_corrections["OH⁻"] += G_OH
end

# ╔═╡ e257853d-2bb2-42f0-b1bd-eedb0736807a
function add_ideal_gas_correction!(thermo_corrections)
	for (; name,frequencies,symmetrynumber,geometry,spin,components) in eachrow(gases)
		if name ∉ ["OH⁻", "e⁻"]
			thermo_corrections[name] += get_thermal_correction_ideal_gas(
				T, frequencies, symmetrynumber, geometry, spin, components
			) * eV
		end
	end
	nothing
end

# ╔═╡ cea31912-c4ff-424f-9000-c537dc8faea4
function add_harmonic_adsorbate_correction!(thermo_corrections)
	for (; name, frequencies) in eachrow(adsorbates)
		thermo_corrections[name] += get_thermal_correction_adsorbate(T, frequencies) * eV
		#println("$sp: $(E_f[sp] / eV)")
	end
	thermo_corrections["COOH-H₂O-e⁻_t"] += (thermo_corrections["COOH_t"] + thermo_corrections["CO_t"]) / 2
	nothing
end

# ╔═╡ b69719cc-9378-4eb5-a6e4-56bd4f250cc5
begin
	for sp in keys(thermo_corrections)
		thermo_corrections[sp] = 0.0
	end
	add_ideal_gas_correction!(thermo_corrections)
	add_harmonic_adsorbate_correction!(thermo_corrections)
	
	function compute_energies_explicit!(Gf, σ, ϕ_we, ϕ, local_pH)

		add_simple_electrochem_corrections!(Gf, ϕ_we, ϕ)
		add_hbond_corrections!(Gf)
		add_hbond_surface_charge_corrections!(Gf, σ)
		add_pH_correction!(Gf, thermo_corrections, local_pH)
	
		for species in [gases, adsorbates, transition_states, surfaces, fictious]
			for (; name, energy) in eachrow(species)
				Gf[name] += energy + thermo_corrections[name]
			end
		end	
		nothing
	end
end

# ╔═╡ 77b9a080-7a3b-4835-9d69-afa13eff3757
md"""
#### Explicit Rate Calculations
"""

# ╔═╡ 90d20692-ce49-419e-8377-9e2b72d176a6
explicit_df = let
	Gf = Dict(zip(
			keys(thermo_corrections), zeros(Float64, length(thermo_corrections))
		))
	ns = names(catmap_df)
	df = DataFrame(fill(Float64[], (length(ns),)), ns)
	
	for ϕ in range(ϕ_we, 0.0, 11)
		for k in keys(Gf)
			Gf[k] = 0.0
		end
		σ = (ϕ_we - ϕ - ϕ_pzc) * C_gap
		params = Dict("ϕ" => ϕ, "ϕ_we" => ϕ_we, "σ" => σ)
		compute_energies_explicit!(Gf, σ, ϕ_we, ϕ, local_pH)
		Gf["g"] = 0.0
		for (k, v) in pairs(Gf)
			Gf[k] = v / eV
		end
		push!(df, merge(params, Gf))
	end
	df
end

# ╔═╡ 3215092e-a913-4c73-ac01-fc725beccb05
catmap_df .- explicit_df

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LessUnitful = "f29f6376-6e90-4d80-80c9-fb8ec61203d5"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

[compat]
DataFrames = "~1.6.1"
LessUnitful = "~0.6.1"
PyCall = "~1.96.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0-beta1"
manifest_format = "2.0"
project_hash = "33e2c182cb2ae1d2c697d94218cd0571206bd9f4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "8c86e48c0db1564a1d49548d3515ced5d604c408"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.9.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LessUnitful]]
deps = ["PhysicalConstants", "Unitful"]
git-tree-sha1 = "7d7b92d12f4748bc4ed196414f04499fb141a365"
uuid = "f29f6376-6e90-4d80-80c9-fb8ec61203d5"
version = "0.6.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.0.1+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "51d946d38d62709d6a2d37ea9bcc30c80c686801"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.9.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "cd4da9d1890bc2204b08fe95ebafa55e9366ae4e"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "43d304ac6f0354755f1d60730ece8c499980f7ba"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "ff42754a57bb0d6dcfe302fd0d4272853190421f"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.19"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "607c142139151faa591b5e80d8055a15e487095b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.16.3"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═20f09fe8-3d2f-11ee-3749-550c6dbe9f0e
# ╟─44a506d1-b0e8-4bb0-8265-fac97271efff
# ╠═3afa3da7-d14b-43a9-bd43-c0ce0550086f
# ╟─1dc9ffb6-02e2-4278-99c8-d1495818256a
# ╠═b24253cc-e044-48ad-8d1c-f12a6765c398
# ╟─15b52261-ba3e-42ba-96a9-c772626e5496
# ╟─b66f7fe9-009d-4287-9d89-a36dbde1e434
# ╟─ea5c99f8-cbd5-4542-ab17-f2582b49e565
# ╟─29b9eacb-ef85-400d-8301-2e59d541b70f
# ╟─7787f977-66fb-4f0f-8ef2-41c78a518933
# ╠═ce88f49d-2ea6-40fd-8e21-9fcecc938e02
# ╠═7d4c60a7-8a5b-47f1-a91e-bb1b2c76c05d
# ╟─26cf49cc-7a4e-4232-85b2-b370db46378e
# ╟─ba894fa1-fba8-40db-9459-7ec78462a246
# ╠═1f020e47-e4cb-4202-80f7-39421e4ce471
# ╠═9d4a2505-9c24-4ed2-bf02-f6598f296f11
# ╠═6120ec04-80c8-423c-a92a-44583ba9bac7
# ╟─8a87380d-20e6-4c27-950b-25003b98bf30
# ╟─30c7d656-e197-495e-bd71-9c71214d1857
# ╠═0b3c9e82-c87b-4f3c-b504-af10d355c0e2
# ╠═dc768d05-de5a-4c44-b30a-7aea60fb1926
# ╟─119b4542-41e6-49fc-a281-a0c2c1b5981f
# ╠═36a30e9a-6fce-40d3-a36a-aaef3d8b7f2a
# ╠═5ccf3fa5-df4d-47e5-a56d-82ba4b099b92
# ╠═e257853d-2bb2-42f0-b1bd-eedb0736807a
# ╠═cea31912-c4ff-424f-9000-c537dc8faea4
# ╠═b69719cc-9378-4eb5-a6e4-56bd4f250cc5
# ╟─77b9a080-7a3b-4835-9d69-afa13eff3757
# ╠═90d20692-ce49-419e-8377-9e2b72d176a6
# ╠═3215092e-a913-4c73-ac01-fc725beccb05
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
