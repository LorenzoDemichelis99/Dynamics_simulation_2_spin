# Simulation of the dynsmics of the 2-spin model with additive noise on a generic graph 

## Structure of the files

This repository contains a Jupyter notebook and some Julia files to implement the simulation of the 2-spin model on a generic graph, with the possibility of introducing also an additive noise. 
The file structure is the following: simulation_data_2_spin.jl contains the code that defines the mutable struct used to store all the important information about the simulation and the code, together with the constructor of the struct and the functions that initialize the relevant quantities of the simulation; simulation_2_spin.jl contains the function implementing the actual simulation; simulation_tools_2_spin.jl contains all the utility functions required to plot and save the results of the simulation.


## Dynamics implemented on the graph

Let $G=(V,E)$ be the graph on which the dynamics is defined, with $V$ being the node set and $E$ being the edge set; the evolution of the degree of freedom $x_{i}(t)$, with $i \in V$, is determined by the following stochastic differential equation:
```math
\begin{equation}
    \frac{dx_{i}}{dt} = f[x_{i}(t)] + \sum_{j \in \partial i} J_{ij} x_{j}(t) + \eta_{i}(t)
\end{equation}
```
with $f[x_{i}(t)]$ being the local term of the dynamics, $\partial i$ being the neighborhood of node $i$, $J_{ij}$ being the coupling constant between $x_{i}(t)$ and $x_{j}(t)$ and $\eta_{i}(t)$ being the noise, which has the following properties:
```math
\begin{align*}
    & \langle \eta_{i}(t) \rangle = 0 \\
    & \langle \eta_{i}(t) \eta_{j}(t') \rangle = 2 D \delta_{i,j} \delta(t-t')
\end{align*}
```
