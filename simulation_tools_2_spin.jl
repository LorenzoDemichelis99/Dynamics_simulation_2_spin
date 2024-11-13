# this file contatins all the small functions required to implement the simulation of the dynamics

# this function generates a realization of the noises on all the nodes in the graph
function noise_values(g::Vector{H}, x::Vector{F}) where F <: AbstractFloat where H <: Function
    return x .|> g
end

# this function plots some trajectories of the dynamics together with the average trajectory
function plot_trajectories(sim_res::simulation_data_2_spin, Nsamples::Int64)
    # allocating the resources
    indices = Vector{Int64}(undef, Nsamples)
    avg_trajectory = Vector{Float64}(undef, sim_res.T)

    # computing the average trajectory
    avg_trajectory[:] = reshape(mean(sim_res.x;dims=1), sim_res.T)

    # sampling the indices of the trajectories
    indices[:] = sample(1:nv(sim_res.G), Nsamples; replace=false)

    # plotting the trajectories
    for i in 1:Nsamples
        if i == 1
            plot(1:sim_res.T, sim_res.x[indices[i],:], color="blue", label=false, xlabel=L"T", ylabel=L"\langle x(t) \rangle", title=L"Plot \,\, of \,\, the \,\, trajectories", size=(1200,600), linewidth=1.0, margin=5mm)
        else
            plot!(1:sim_res.T, sim_res.x[indices[i],:], color="blue", label=false, linewidth=1.0)
        end
    end

    # plotting the average trajectory
    plot!(1:sim_res.T, avg_trajectory, color="red", label="average trajectory", linewidth=1.2)
end

# this function plots the Lagrange multiplier of the dynamics
function plot_lgm(sim_res::simulation_data_2_spin)
    plot(1:sim_res.T, sim_res.λ, size=(1200,600), color="orange", xlabel=L"T", ylabel=L"\lambda(t)", linewidth=1.2, title=L"Plot \, \, of \, \, the \, \, Lagrange \, \, multiplier", label=false, margin=5mm)
end

# this function computes some meaningful statistical quantities of the dynamics
function stats_analysis(sim_res::simulation_data_2_spin)
    # allocation of the resources
    μ = Vector{Float64}(undef, sim_res.T)
    C = Matrix{Float64}(undef, (sim_res.T, sim_res.T))
    ρ = Matrix{Float64}(undef, (sim_res.T, sim_res.T))

    # computation of the mean
    μ = reshape(mean(sim_res.x;dims=1), sim_res.T)

    # computation of the covariance matrix
    C = cov(sim_res.x;dims=1)

    # computation of the correlation matrix
    ρ = cor(sim_res.x;dims=1)

    return μ, C, ρ
end

# this function plots the mean function μ(t), the covariance function C(t-t₀,t₀) and the correlation function ρ(t-t₀,t₀) in the interval [t₀,tₙ]
function plot_stats(μ::Vector{Float64}, C::Matrix{Float64}, ρ::Matrix{Float64}, t₀::Int64, tₙ::Int64)
    plot(t₀:tₙ, μ[t₀:tₙ], layout=(3,1), subplot=1, xlabel=L"T", ylabel=L"\mu(t)", size=(1200,800), color="green", margin=5mm, label=false)
    plot!(t₀:tₙ, C[t₀:tₙ,t₀], subplot=2, xlabel=L"T", ylabel=L"C(t-t_{0},t_{0})", color="red", label=false, margin=5mm)
    plot!(t₀:tₙ, ρ[t₀:tₙ,t₀], subplot=3, xlabel=L"T", ylabel=L"\rho(t-t_{0},t_{0})", color="orange", label=false, margin=5mm)
end

# this function saves the trajectories obtained from the simulation
function save_trajectories(sim_res::simulation_data_2_spin, folder_path::String)
    # creating the path for the file
    filename = joinpath(folder_path, "x.csv")

    # saving the trajectories in the desired path
    f = open(filename, "w")
    writedlm(f, sim_res.x)
    close(f)
end

# this function saves the results of the statistical analysis
function save_stats(μ::Vector{Float64}, C::Matrix{Float64}, ρ::Matrix{Float64}, folder_path::String)
    # creating the paths for the files
    filename_μ = joinpath(folder_path, "μ.csv")
    filename_C = joinpath(folder_path, "C.csv")
    filename_ρ = joinpath(folder_path, "ρ.csv")

    # saving the results of the statistical analysis in the desired path
    f = open(filename_μ, "w")
    writedlm(f, μ)
    close(f)

    f = open(filename_C, "w")
    writedlm(f, C)
    close(f)

    f = open(filename_ρ, "w")
    writedlm(f, ρ)
    close(f)
end