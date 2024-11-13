# this file implements the simulation of the dynamics

# The function implementing the simulation takes the following parameters:
# G: the underlying graph of the dynamics
# T: the number of time steps
# Δ: the time step

function simulation_2_spin(G::SimpleGraph, T::Int64, Δ::F, gx::Function, jx::Function, x0_init::Function, noise::Bool) where F <: AbstractFloat
    # initialization
    sim = simulation_data_2_spin(G, T, Δ, gx, jx, x0_init)
    N = nv(G)
    if noise == true
        current_η = Vector{Float64}(undef, N)
        previous_η = Vector{Float64}(undef, N)
    end

    # setting the progress meter
    p = Progress(T; showspeed=true)

    # computation of the initial value of the Lagrange multiplier
    if noise == true
        current_η[:] = noise_values(sim.g, sim.x[:,1])
        λ_temp_1 = 0.0
        λ_temp_2 = 0.0
        for i in 1:N
            λ_temp_2 += sim.x[i,1] * current_η[i]
            for k in 1:N
                λ_temp_1 += sim.J[i,k] * sim.x[k,1] * sim.x[i,1]
            end
        end
        sim.λ[1] = (1/N) * (λ_temp_1 + λ_temp_2)
        previous_η[:] = copy(current_η)
    else
        sim.λ[1] = 0.0
        for i in 1:N
            for k in 1:N
                sim.λ[1] += sim.J[i,k] * sim.x[i,1] * sim.x[k,1]
            end
        end
        sim.λ[1] /= N
    end

    # simulation of the dynamics
    @inbounds for t in 2:T
        if noise == true
            # generation of the new values of the trajectories 
            for i in 1:N
                sim.x[i,t] = sim.x[i,t-1] - (sim.λ[t-1] * sim.x[i,t-1] * Δ) + previous_η[i]
                for k in 1:N
                    sim.x[i,t] += Δ * sim.J[i,k] * sim.x[k,t-1]
                end
            end

            # generation of a new instance of the noise
            current_η[:] = noise_values(sim.g, sim.x[:,1])

            # generation of the new value of the Lagrange multiplier
            λ_temp_1 = 0.0
            λ_temp_2 = 0.0
            for i in 1:N
                λ_temp_2 += sim.x[i,t] * current_η[i]
                for k in 1:N
                    λ_temp_1 += sim.J[i,k] * sim.x[k,t] * sim.x[i,t]
                end
            end
            sim.λ[t] = (1/N) * (λ_temp_1 + λ_temp_2)

            previous_η[:] = copy(current_η)
        else
            # generation of the new values of the trajectories
            for i in 1:N
                sim.x[i,t] = sim.x[i,t-1] - (sim.λ[t-1] * sim.x[i,t-1] * Δ) 
                for k in 1:N
                    sim.x[i,t] += Δ * sim.J[i,k] * sim.x[k,t-1]
                end
            end

            # generation of the new value of the Lagrange multiplier 
            sim.λ[t] = 0.0
            for i in 1:N
                for k in 1:N
                    sim.λ[t] += sim.J[i,k] * sim.x[k,t-1] * sim.x[i,t-1]
                end
            end
            sim.λ[t] /= N
        end
        next!(p)
    end

    return sim
end