module CO2RCell_wo_microkinetics
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors 
using StaticArrays
using InteractiveUtils

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using PyCall

function main(;nref=0,
              voltages=(-1.5:0.1:-0.6)*ufac"V",
              scheme=:μex,
              κ=10.0,
              Plotter=PyPlot,
              kwargs...)

    @local_phconstants N_A e R ε_0 c_0 h
    @local_unitfactors cm μF mol dm s mA A nm bar eV cm m K μA Pa


    
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


    sigmas = []
    energies = []

    ϕs = []
    rs = []

    T   = 273.15 + 25 * ufac"K"
    pH  = 6.8

    #Hcp_CO  = 9.7e-6 * mol/(m^3 * Pa)
    Hcp_CO2 = 3.3e-4 * mol/(m^3 * Pa)

    E_ads_CO2 = 0.657600203 * eV
    frequencies = [136.85, 183.6, 212.95, 250.7, 306.0, 510.55, 562.25, 1176.05, 1889.85] * h * c_0 / cm / eV

    py"""
    from ase.thermochemistry import HarmonicThermo

    def get_thermal_correction_adsorbate(T, frequencies):
        thermo = HarmonicThermo(frequencies)
        return thermo.get_helmholtz_energy(T, verbose=False)
    """
    harmonic_adsorbate_correction   = py"get_thermal_correction_adsorbate"(T, frequencies) * eV

    C_gap = 20 * μF/cm^2
    ϕ_pzc = 0.16 * ufac"V"
    
    ikplus      = 1
    ihco3       = 2
    ico3        = 3
    ico2        = 4
    ico         = 5
    iohminus    = 6
    ihplus      = 7


    ϕ_curr = Ref(100.0)
    function halfcellbc(f,u,bnode,data)
        (;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,RT, ε)=data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)


        if bnode.region==Γ_we

            #boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            
            # Robin b.c. for the Poisson equation

            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)

            # Flux conditions for CO2 and CO
            prefactor   = 1.0e13
            γ_CO2       = 1.0
            a_CO2       = γ_CO2 * u[ico2] / Hcp_CO2 / bar # * ufac"dm^3/mol" #
            θ_free      = 0.9999
            S           = 9.61e-5 / N_A * (1.0e10)^2 * mol/m^2
            
            sigma                           = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
            electrochemical_correction      = [-0.000286600929, 0.0297720125]' * [(sigma / ufac"μA/cm^2")^2 , (sigma / ufac"μA/cm^2")] * eV
            
            ΔG  = (E_ads_CO2 + harmonic_adsorbate_correction + electrochemical_correction) * N_A            
            r   = prefactor * a_CO2 * θ_free * exp(-ΔG / RT - 2.3*pH) * S
            
            if ϕ_we == ϕ_curr[]
                pop!(sigmas)
                pop!(energies)
                pop!(rs)
            else
                ϕ_curr[] = ϕ_we
                push!(ϕs, ϕ_we)
            end

            push!(sigmas, sigma)
            push!(energies, ΔG / (eV * N_A))
            #push!(rs, r)
            push!(rs, u[ico2] * dm^3/mol)

            f[ico2] = r
            f[ico]  = -r
            f[iohminus] = -2*r

        end
        nothing
    end
    
    
    celldata=ElectrolyteData(;nc=7,
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
    
    ## Calculate current density-voltage curve
    volts, currs, sols = ivsweep(cell; voltages, ispec=iohminus, kwargs...)
    vis=GridVisualizer(;Plotter, layout=(2,1))
    scalarplot!(vis[1,1], volts, currs*cm^2/mA,color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP",clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(mA/cm^2)", yscale=:log)
    #scalarplot!(vis[2,1], sigmas, energies, color="black",clear=true,xlabel="σ/(μC/cm^s)",ylabel="ΔE/eV")
    scalarplot!(vis[2,1], ϕs, rs, xlimits=(-1.5,-0.6), yscale=:log, xlabel="Δϕ/V", ylabel="c(CO2)/M")
    for (volt,curr) in zip(volts,currs)
        println("$volt, $curr")
    end
    return reveal(vis)
end

end

#=
```@example CO2RCell_wo_microkinetics
CO2RCell_wo_microkinetics.main()
```
=#

