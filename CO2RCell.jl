module Example111_CO2RCell
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors 
using StaticArrays
using InteractiveUtils
using ForwardDiff

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using PyCall

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

function main(;nref=0,
              voltages=(-1.5:0.1:-0.6)*ufac"V",
              molarities=[0.001,0.01,0.1,1],
              scheme=:μex,
              xmax=1,
              κ=10.0,
              Plotter=PyPlot,
              new=false,
              kwargs...)

    @local_phconstants N_A e R ε_0 k_B
    F = N_A*e
    c_0 = 2.99792458e8
    @local_unitfactors cm μF mol dm s mA A nm bar eV μA


    
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose="e",
              reltol=1.0e-8,
              tol_mono=1.0e-10)

    kwargs=merge(defaults, kwargs) 

    hmin=1.0e-1*ufac"μm"*2.0^(-nref)
    hmax=1.0*ufac"μm"*2.0^(-nref)
    L=80.0 * ufac"μm"
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)

    # sigmas = []
    # energies = []

    # ϕs = []
    # rs = []

    T   = 273.15 + 25 * ufac"K"
    pH  = 6.8


    # Henry constants
    Hcp_CO  = 9.7e-6 * ufac"mol/(m^3 * Pa)"
    Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"

    # kinetic model
    # 'CO2_g + 2*_t <-> CO2*_t',	                  #1
    # 'CO2*_t + H2O_g + ele_g <-> COOH*_t + OH_g',  #2
    # 'COOH*_t + H2O_g + ele_g <-> COOH-H2O-ele_t <-> CO*_t + H2O_g + OH_g + *_t; beta=0.5', #3
    # 'CO*_t <-> CO_g + *_t',	                      #4


    bulk_species        = ["H2_g", "CO2_g", "CO_g", "OH_g", "H2O_g", "ele_g"]
    surface_species     = ["CO2*_t", "COOH*_t", "CO*_t"]
    transition_state    = ["COOH-H2O-ele_t"]


    ## formation energies at 0 K from DFT calculations
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

    electro_correction_params = Dict(
        "CO2*_t"  => [-0.000286600929 / (μA/cm^2)^2, 0.0297720125 / (μA/cm^2)],
        "COOH*_t" => [-9.0295682e-05 / (μA/cm^2)^2, 0.00226896383 / (μA/cm^2)],
        "CO*_t"   => [-0.000189106972 / (μA/cm^2)^2,-0.00942574086 / (μA/cm^2)]
    )

    ## thermodynamics corrections
    thermo_corrections = Dict(zip(keys(E_f), zeros(length(E_f))))

    ideal_gas_params = Dict(
        "H2O_g" => (2, "nonlinear", 0, "H2O"),
        "CO_g" => (1, "linear", 0, "CO"),
        "CO2_g" => (2, "linear", 0, "CO2"),
        "H2_g" => (2, "linear", 0, "H2")
    )
    ### ideal gases
    for sp in bulk_species
        if sp ∉ ["OH_g", "ele_g"]
            thermo_corrections[sp] += py"get_thermal_correction_ideal_gas"(T, frequencies[sp], ideal_gas_params[sp]...) * eV
        end
    end

    ### harmonic adsorbates
    for sp ∈ surface_species
        thermo_corrections[sp] += py"get_thermal_correction_adsorbate"(T, frequencies[sp]) * eV
        println("$sp: $(E_f[sp] / eV)")
    end

    ### BEP scaling for transition states
    thermo_corrections["COOH-H2O-ele_t"] += (thermo_corrections["COOH*_t"] + thermo_corrections["CO*_t"]) / 2

    ## _get_echem_corrections
    G_H2O       = E_f["H2O_g"] + thermo_corrections["H2O_g"]
    G_H2        = E_f["H2_g"] + thermo_corrections["H2_g"]
    G_H         = 0.5 * G_H2 - .0592 * pH / 298.14 * T * eV
    G_OH        = G_H2O - G_H
    thermo_corrections["OH_g"] += G_OH


    C_gap = 20 * ufac"μF/cm^2"
    ϕ_pzc = 0.16 * ufac"V"
    
    ikplus      = 1
    ihco3       = 2
    ico3        = 3
    ico2        = 4
    ico         = 5
    iohminus    = 6
    ihplus      = 7
    ico_ad      = 8
    ico2_ad     = 9
    icooh_ad    = 10


    function halfcellbc(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
        (;nc,na,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,T, RT, ε)=data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)


        if bnode.region==Γ_we

            #boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            
            # Robin b.c. for the Poisson equation

            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)

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
            sigma = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
            #sigma = C_gap *(ϕ_we - ϕ_pzc)
            for sp in surface_species
                electro_corrections[sp] += electro_correction_params[sp]' * [sigma^2, sigma] * eV
            end


            G_f = Dict(zip([bulk_species; surface_species; transition_state], zeros(Tval, nc + na + length(transition_state))))
            for sp in [bulk_species; surface_species; transition_state]
                G_f[sp] += E_f[sp] + thermo_corrections[sp] + electro_corrections[sp]
            end


            # rate constants
            G_IS = zeros(Tval, 4)
            G_FS = zeros(Tval,4)
            G_TS = zeros(Tval, 4)

            kf = zeros(Tval, 4)
            kr = zeros(Tval, 4)

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

            # rates
            # 'CO2_g + 2*_t <-> CO2*_t',	                  #1
            # 'CO2*_t + H2O_g + ele_g <-> COOH*_t + OH_g',  #2
            # 'COOH*_t + H2O_g + ele_g <-> COOH-H2O-ele_t <-> CO*_t + H2O_g + OH_g + *_t; beta=0.5', #3
            # 'CO*_t <-> CO_g + *_t',	                      #4
            S       = 9.61e-5 / N_A * (1.0e10)^2 * ufac"mol/m^2"
            θ_free  = (1- u[ico2_ad] - u[ico_ad] - u[icooh_ad])

            rates = zeros(Tval, 4)

            rates[1] = kf[1] * (u[ico2] / Hcp_CO2 / bar) * θ_free^2 - kr[1] * (u[ico2_ad] )
            rates[2] = kf[2] * (u[ico2_ad] ) * 1.0 * 1.0 - kr[2] * (u[icooh_ad] ) * (u[iohminus])
            rates[3] = kf[3] * (u[icooh_ad] ) * 1.0 * 1.0 - kr[3] * (u[ico_ad] ) * 1.0 * (u[iohminus]) * θ_free
            rates[4] = kf[4] * (u[ico_ad] )  - kr[4] * (u[ico] / Hcp_CO / bar) * θ_free

            #println("rate constants: $(ForwardDiff.value.(kf)) and $(ForwardDiff.value.(kr))")
            println("rates: $(ForwardDiff.value.(rates))")

            # bulk species
            f[ico] += -rates[4] * S
            f[ico2] += rates[1] * S
            f[iohminus] += -rates[2] * S - rates[3] * S
            
            # surface species
            f[ico2_ad] += -rates[1] + rates[2]
            f[ico_ad] += -rates[3] + rates[4]
            f[icooh_ad] += -rates[2] + rates[3]

        end
        nothing
    end
    
    
    celldata=ElectrolyteData(;nc=7,
                             na=3,
                             z=[1,-1,-2,0,0,-1,1],
                             D=[1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s", # from Ringe paper
                             T=T,
                             eneutral=false,
                             κ=fill(κ,7),
                             Γ_we=1,
                             Γ_bulk=2,
                             scheme)

    (;iϕ::Int,ip::Int)=celldata
    
    celldata.c_bulk[ikplus]         = 0.1 * mol/dm^3
    celldata.c_bulk[ihco3]          = (0.1 - 9.53936e-8) * mol/dm^3
    celldata.c_bulk[ico3]           = 9.53936e-8 * mol/dm^3
    celldata.c_bulk[ico2]           = 0.033 * mol/dm^3
    celldata.c_bulk[ico]            = 0.0 * mol/dm^3
    celldata.c_bulk[iohminus]       = 10^(pH - 14) * mol/dm^3
    celldata.c_bulk[ihplus]         = 10^(-pH) * mol/dm^3

    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-10)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)
    
        
    volts, currs, sols = ivsweep(cell; voltages, ispec=iohminus, kwargs...)
    vis=GridVisualizer(;Plotter, layout=(1,1))
    scalarplot!(vis[1,1], volts, currs*ufac"cm^2/mA",color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP",clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(mA/cm^2)", yscale=:log)
    #scalarplot!(vis[2,1], sigmas, energies, color="black",clear=true,xlabel="σ/(μC/cm^s)",ylabel="ΔE/eV")
    #scalarplot!(vis[2,1], ϕs, rs, xlimits=(-1.5,-0.6), yscale=:log, xlabel="Δϕ/V", ylabel="c(CO2)/M")
    for (volt, curr) in zip(volts, currs)
        println("$volt,$curr")
    end
    return reveal(vis)
end

end

#=
```@example Example111_CO2RCell
Example111_CO2RCell.main()
```
=#
