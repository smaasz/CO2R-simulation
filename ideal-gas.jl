using PyCall
using LessUnitful
using LinearAlgebra

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

py"""
from ase.collections import g2

def moleculedata(formula):
    a = g2[formula]
    return a.numbers, a.get_masses(), a.get_positions()
"""

@enum MoleculeGeometry monoatomic linear nonlinear

@kwdef struct IdealGas
    elements::Vector{Int}
    masses::Vector{Float64}
    positions::Matrix{Float64}
    symmetrynumber::Int32
    frequencies::Vector{Float64}
    spin::Float64
    geometry::MoleculeGeometry
    temperature::Float64
end

function totalmass(idealgas::IdealGas)
    sum(idealgas.masses)
end

function momentsofinertia(idealgas::IdealGas)
    m = idealgas.masses
    x = idealgas.positions
    M = sum(m)
    
    cm = sum(x .* m, dims=1) / M
    x .-= cm

    I = x' * (x .* m) - diagm(diag((x .* m) * x'))
    eigvals(I)
end

function entropy(idealgas::IdealGas)
    @local_phconstants k_B h R ħ
    @local_unitfactors bar eV

    N = length(idealgas.elements)
    M = totalmass(idealgas)
    T = idealgas.temperature
    σ = idealgas.symmetrynumber
    geometry = idealgas.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 6)]
    end
    S = idealgas.spin
    (IA, IB, IC) = momentsofinertia(idealgas)
    if geometry == linear
        @assert isapprox(IB, 0.0, atol = 1e-6) && isapprox(-IA, IC, rtol = 1e-6)
    end
    P° = 1 * bar

    translational = k_B * (log((2 * π * M * k_B * T / h^2)^(3/2) * k_B * T / P°) + 5 / 2)
    rotational = 0.0
    if geometry == linear
        rotational +=  k_B * (log(2 * IC * k_B * T / (ħ^2 * σ)) + 1)
    elseif geometry == nonlinear
        rotational += k_B * (log(sqrt(π * IA * IB * IC) / σ * (2 * k_B * T / ħ^2)^(3/2)) + 3/2)
    end
    vibrational = k_B * sum(h * ω / (k_B * T * (exp.(h * ω / (k_B * T)) .- 1)) .- log.(1 .- exp.(-h * ω / (k_B * T))))
    electronic = k_B * log(2 * S + 1)

    #println("$(translational/eV), $(rotational/eV), $(vibrational/eV), $(electronic/eV)")
    translational + rotational + vibrational + electronic
end

function enthalpy(idealgas::IdealGas)
    @local_phconstants k_B h R
    @local_unitfactors eV

    N = length(idealgas.elements)
    T = idealgas.temperature
    geometry = idealgas.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 6)]
    end

    zpe = h * sum(ω) / 2
    translational = 3 / 2 * k_B * T
    rotational = 0
    if geometry == linear
        rotational += k_B * T
    elseif  geometry == nonlinear
        rotational += 3 / 2 * k_B * T
    end
    vibrational = h * sum(ω / (exp.(h * ω / (k_B * T)) .- 1))
    corr = k_B * T

    #println("$zpe, $translational, $rotational, $vibrational, $corr")
    zpe + translational + rotational + vibrational + corr
end

function test_CO2()
    @local_phconstants c_0 h
    @local_unitfactors cm eV u Å

    (elements, masses, positions) = py"moleculedata"("CO2")
    masses *= u
    positions *= Å
    symmetrynumber = 2
    frequencies = [24.1, 70.7, 635.8, 640.5, 1312.2, 2361.2] * c_0 / cm # * ufac"h / cm / eV" * c_0
    spin = 0
    geometry = linear
    temperature = 298.14

    co2 = IdealGas(elements, masses, positions, symmetrynumber, frequencies, spin, geometry, temperature)
    H = enthalpy(co2)
    S = entropy(co2)

    G = H - temperature * S

    G_py = py"get_thermal_correction_ideal_gas"(temperature, frequencies * h / eV, symmetrynumber, "linear", spin, "CO2")

    G, G_py * eV
end
