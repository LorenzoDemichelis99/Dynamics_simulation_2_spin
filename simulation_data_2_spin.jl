# this file is used to define the struct which contains all the important information about the simulation of the dynamics of the s-spin model

# The important quantities for the simulations are:
# G: underlying graph of the system
# J: matrix of the coupling constants
# D: parameter of the noise term
# T: number of time steps of the simulation
# Δ: time steps
# x: matrix containing the trajectories generated in the simulation
# λ: vector containing the trajectory of the Lagrange multiplier

mutable struct simulation_data_2_spin{F<:AbstractFloat,H<:Function}
    G::SimpleGraph
    J::Matrix{F}
    g::Vector{H}
    T::Int
    Δ::F
    x::Matrix{F}
    λ::Vector{F}
end

# this is the constructor of the struct defined above
function simulation_data_2_spin(G::SimpleGraph, T::Int, Δ::F, gx::Function, jx::Function, x0_init::Function) where F<:AbstractFloat
    # extraction of the number of vertices
    n = nv(G)

    # allocation of the required memory
    J = Matrix{F}(undef, (n,n))
    g = Vector{Function}(undef, n)
    x = Matrix{F}(undef, (n,T))
    λ = Vector{F}(undef, T)

    # initialization of the required quantities
    coupling!(G, J, jx)
    gterm!(g, Δ, gx)
    init_x!(x, x0_init)

    # initialization of the struct
    simulation_data_2_spin(G, J, g, T, Δ, x, λ)
end

# this is the function that builds the matrix of couplings
function coupling!(G::SimpleGraph, J::Matrix{F}, jx::Function) where F <: AbstractFloat
    J .= 0.0
    @inbounds for i in 1:size(J, 1)
        J[i,i] = 0.0
        for j in i+1:size(J, 2)
            if has_edge(G,i,j)
                J[i,j] = J[i,j] |> jx
                J[j,i] = copy(J[i,j])
            else
                J[i,j] = 0.0
            end
        end
    end
end

# this is the function that builds the noise terms
function gterm!(g::Vector{H}, Δ::F, gx::Function) where F <: AbstractFloat where H <: Function
    @inbounds @Threads.threads for i in 1:length(g)
        g[i] = x -> sqrt(2 * (x |> gx) * Δ) * randn()
    end
end

# this is the function that initializes the trajectories of the degrees of freedom
function init_x!(x::Matrix{F}, x0_init::Function) where F <: AbstractFloat
    @inbounds @Threads.threads for i in 1:size(x,1)
        x[i,1] = x0_init(x[i,1])
    end
end