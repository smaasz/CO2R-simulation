module Example112_CO2RCell
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors 
using StaticArrays
using InteractiveUtils
using ForwardDiff
using CatmapInterface

mutable struct CatmapData
    rate_constant_fns
    compute_rates
end

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


function main(;nref=0,
              voltages=(-1.5:0.1:-0.6)*ufac"V",
              scheme=:μex,
              κ=10.0,
              Plotter=PyPlot,
              kwargs...)

    @local_phconstants N_A e R ε_0 k_B c_0
    F = N_A*e
    @local_unitfactors cm μF mol dm s mA A nm bar eV μA μm


    
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose="e",
              reltol=1.0e-8,
              tol_mono=1.0e-10)

    kwargs=merge(defaults, kwargs) 

    hmin    = 1.0e-1*μm*2.0^(-nref)
    hmax    = 1.0*μm*2.0^(-nref)
    L       = 80.0 * μm
    X       = geomspace(0,L,hmin,hmax)
    grid    = simplexgrid(X)

    # environment parameters
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

    catmap_setup_file   = "catmap_CO2R_template.mkm"
    catmap_energy_file  = "catmap_CO2R_energies.txt"

    instance_input_file_name                = instantiate_catmap_template(catmap_setup_file, 0.0)
    ~, rate_constant_fns, ~, compute_rates  = get_catmap_output(instance_input_file_name, catmap_energy_file, true)
    catmap_data                             = CatmapData(rate_constant_fns, compute_rates)  

    function halfcellbc(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
        (;nc,na,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,T, RT, ε)=data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)


        if bnode.region==Γ_we

            #boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            
            # Robin b.c. for the Poisson equation

            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)


            sigma = C_gap / ufac"μF/cm^2" * (ϕ_we - u[iϕ] - ϕ_pzc)

            println([ForwardDiff.value(kf(sigma)) for kf in catmap_data.rate_constant_fns[1][1]])

            rates = catmap_data.compute_rates(catmap_data.rate_constant_fns[1], [u[ico2_ad], u[icooh_ad], u[ico_ad]], [u[ico2] / Hcp_CO2 / bar, u[ico] / Hcp_CO / bar, 1.0, 1.0, 1.0, 1.0], sigma)

            println("$(ForwardDiff.value.(rates))")

            S       = 9.61e-5 / N_A * (1.0e10)^2 * ufac"mol/m^2"

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
    
    function pre(sol, t)
        
        instance_input_file_name = instantiate_catmap_template(catmap_setup_file, -t)
        ~, rate_constant_fns, ~, compute_rates  = get_catmap_output(instance_input_file_name, catmap_energy_file, true)
        catmap_data.rate_constant_fns           = rate_constant_fns
        catmap_data.compute_rates               = compute_rates

        nothing
    end

        
    tsol, j_we, j_bulk  = ivsweep(cell; voltages, more_pre=pre, kwargs...)
    vis                 = GridVisualizer(;Plotter, layout=(1,1))
    currs               =  [F * j[iohminus] for j in j_we] * ufac"cm^2/mA"
    scalarplot!(vis[1,1], tsol.t, currs, color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP",clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(mA/cm^2)", yscale=:log)
    #scalarplot!(vis[2,1], sigmas, energies, color="black",clear=true,xlabel="σ/(μC/cm^s)",ylabel="ΔE/eV")
    #scalarplot!(vis[2,1], ϕs, rs, xlimits=(-1.5,-0.6), yscale=:log, xlabel="Δϕ/V", ylabel="c(CO2)/M")
    for (volt, curr) in zip(tsol.t, currs)
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
