##Code to look at simple C-R model and how adding constant Pred driven mortality impacts dynamics
##Initiated: Oct 9, 2025
##Contributors: Marie Gutgesell

using Parameters: @with_kw, @unpack ##imports Parameters package that provides convenient macros for working with keyword arugemnts, parameter structs and unpacking variables 
using Plots
using DifferentialEquations
using Symbolics
using ForwardDiff
using LinearAlgebra
using NLsolve
using DataFrames
using Interact
using Statistics


##functional responses between resources and consumer 
function f_R1C(u, p, t)
    @unpack aR_C, hR_C= p  
    R1,  C = u
    return  aR_C * R1 / (1 + aR_C * hR_C * R1)
end

#Set up paramters for first C-R model
@with_kw mutable struct ModelPar_test_1
  
    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 2.0
    K = 5.0
    aR_C = 0.8  ##attack rate  of consumer on R
    e = 0.7   ##energy conversion 
    mC = 0.5 ##C mortality rate
    hR_C = 1.0 ##handling time of C on R
   
    
   ##temporal variation in R parameters
   l1 = 0.1  #magnitude of variation in K of R1 
   pf = 5.0 ##period of fluctuation 
     e1::Function = t -> sin(2π / pf * t)
   
    
    ##C functional responses
    f_r1c::Function = f_R1C  

    ##Constant predator driven mortality
    Z = 0.0
    Y = 0.0
end



##Model 
function model_unforced!(du, u, p ,t)
    @unpack r, K,aR_C, e, mC, hR_C, l1,  e1, Z, Y  = p
   R1, C = u 

## consumer functional responses
  f_r1c = p.f_r1c(u, p, t)

   
   ##ODEs
  # du[1] = r * R1 * (1 - R1 / (K - l1*(e1(t) - 0.5))) - C * f_r1c #temporally forced R 
   du[1] = r * R1 * (1 - R1 / K) - C * f_r1c - Y * R1
   du[2] = e * C * (f_r1c)  - mC * C - Z*C
 
   return du
 end 

 # Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_unforced(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_unforced!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end


###Set up functions for eigenvalue analysis
##trying my own equilibrium function, that follows structure similar to KC 
##other functions for eigenvalue analysis - based on KC code
find_eq(u, p) = nlsolve((du, u) -> model_unforced!(du, u, p, zero(u)), u).zero
cmat(u, p) = ForwardDiff.jacobian(x -> rhs_unforced(x, p), u)


"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
λ1_stability(M) = maximum(real.(eigvals(M)))

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
function λ1_stability_imag(M) 
    ev = eigvals(M)
    imag.(ev)[findall(real.(ev) .== maximum(real.(ev)))[1]]
end 


"""M is the community matrix, we can be calculated with `cmat(u, p)`
Note: `\nu` is the what to input `ν` which looks a bit too much like `v` for my taste
"""
ν_stability(M) = λ1_stability((M + M') / 2)



##trying function to calculate equilibrium that i can then loop over for values of K etc. -using KC approaches
 # Same equilibrium before and after
function equilibrium_unforced(p, t; frac_window=0.1, n_window=1000)
    u0 = [3.0, 1.5]                      # initial condition
    tspan = (0.0, t)

    # simulate dynamics to approach equilibrium
    prob = ODEProblem(model_unforced!, u0, tspan, deepcopy(p))
    sol  = solve(prob, reltol=1e-8, abstol=1e-8)

    # equilibrium via root-finding using terminal state as initial guess
    u_approx = sol(t)                     # state at time t
    eq = nlsolve((du, u) -> model_unforced!(du, u, deepcopy(p), 0.0), u_approx).zero
    ##can ask eq to only detect coexistence equilibrium, and feed that into M #####
    # compute stability metrics from community matrix
    M       = cmat(eq, p)
    λ1      = λ1_stability(M)
    λ1_imag = λ1_stability_imag(M)
    react   = ν_stability(M)

    # ----- summary stats "at equilibrium" -----
    # evaluate solution over a window at the end of the simulation
    t_start = max(0.0, t * (1 - frac_window))
    t_grid  = range(t_start, t; length=n_window)
    sol_grid = sol(t_grid)                # each row = state variable, each column = time sample

    mean_state = vec(mean(sol_grid; dims=2))
    sd_state   = vec(std(sol_grid; dims=2))
    cv_state   = sd_state ./ mean_state   # may produce Inf/NaN if mean≈0 (expected)
    min_state  = vec(minimum(sol_grid; dims=2))
    max_state  = vec(maximum(sol_grid; dims=2))

    return (
        eq        = eq,
        λ1        = λ1,
        λ1_imag   = λ1_imag,
        react     = react,
        mean      = mean_state,
        sd        = sd_state,
        cv        = cv_state,
        min       = min_state,
        max       = max_state
    )
end

##########################################################################################
 
##Solve ODE 
##set initial condition
u0 = [3.0, 1.5]
tspan = (0.0, 500.0)
##set Parameters
p = ModelPar_test_1(K = 4.0, aR_C = 1.7, Z = 0.0, Y = 0.0)

##Define the ODE problem
prob_1 = ODEProblem(rhs_unforced, u0, tspan, p)
sol_1 = solve(prob_1)

##plot timeseries
plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution")

##trade-off b/w K and Z -- increasing K needed to maintain increased consistent mortality on C

####Conduct stability analysis and look at bifurcations 
results_K_all = []
for K in 3.0:0.1:8.0
    p = ModelPar_test_1(K=K)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_K_all, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all)

p = ModelPar_test_1(Z = 1.0)
t = 100
k_sol = equilibrium_unforced(p, t)

plot(df_eq.K, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.K, R1_min, label = "R1 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, R1_max, label = "R1 max")

plot(df_eq.K, C_min, col = "red", label = "C min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, C_max, col = "red", label = "C max")

##over range of aR_C
results_a_all = []
for aR_C in 0.5:0.1:1.5
    p = ModelPar_test_1(aR_C=aR_C)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_a_all, (; aR_C=aR_C, eq_data...))
end
df_eq = DataFrame(results_a_all)

plot(df_eq.aR_C, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.aR_C, R1_min, label = "R1 min", xlabel = "aR_C", ylabel ="min/max")
plot!(df_eq.aR_C, R1_max, label = "R1 max")

plot(df_eq.aR_C, C_min, col = "red", label = "C min", xlabel = "aR_C", ylabel ="min/max")
plot!(df_eq.aR_C, C_max, col = "red", label = "C max")





##over range of e
results_e_all = []
for e in 0.6:0.1:1.0
    p = ModelPar_test_1(e=e)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_e_all, (; e=e, eq_data...))
end
df_eq = DataFrame(results_e_all)

plot(df_eq.e, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.e, R1_min, label = "R1 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, R1_max, label = "R1 max")

plot(df_eq.e, C_min, col = "red", label = "C min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, C_max, col = "red", label = "C max")

##over range of m
results_m_all = []
for mC in 0.4:0.1:1.5
    p = ModelPar_test_1(mC=mC)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_m_all, (; mC=mC, eq_data...))
end
df_eq = DataFrame(results_m_all)

plot(df_eq.mC, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.mC, R1_min, label = "R1 min", xlabel = "m", ylabel ="min/max")
plot!(df_eq.mC, R1_max, label = "R1 max")

plot(df_eq.mC, C_min, col = "red", label = "C min", xlabel = "m", ylabel ="min/max")
plot!(df_eq.mC, C_max, col = "red", label = "C max")

##over range of Z
results_Z_all = []
for Z in 0.0:0.1:2.0
    p = ModelPar_test_1(Z=Z, K = 7.0)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_Z_all, (; Z=Z, eq_data...))
end
df_eq = DataFrame(results_Z_all)

plot(df_eq.Z, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.Z, R1_min, label = "R1 min", xlabel = "Z", ylabel ="min/max")
plot!(df_eq.Z, R1_max, label = "R1 max")

plot(df_eq.Z, C_min, col = "red", label = "C min", xlabel = "Z", ylabel ="min/max")
plot!(df_eq.Z, C_max, col = "red", label = "C max")

##over range of Y
results_Y_all = []
for Y in 0.0:0.1:1.5
    p = ModelPar_test_1(Y=Y, K = 5.0)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_Y_all, (; Y=Y, eq_data...))
end

df_eq = DataFrame(results_Y_all)

plot(df_eq.Y, df_eq.λ1)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
C_cv = [row.cv[2] for row in eachrow(df_eq)]
C_max = [row.max[2] for row in eachrow(df_eq)]
C_min = [row.min[2] for row in eachrow(df_eq)]


##Looking at bifurcations
plot(df_eq.Y, R1_min, label = "R1 min", xlabel = "Y", ylabel ="min/max")
plot!(df_eq.Y, R1_max, label = "R1 max")

plot(df_eq.Y, C_min, col = "red", label = "C min", xlabel = "Y", ylabel ="min/max")
plot!(df_eq.Y, C_max, col = "red", label = "C max")



##Looking at both K and Z
results_Z_K_all_uf = []
for K in 3.0:0.5:8.0,  Z in 0.0:0.1:1.0
    p = ModelPar_test_1(Z = Z,  K = K, aR_C = 1.5)
    t = 100
    eq_data = equilibrium_unforced(p, t)
    push!(results_Z_K_all_uf, (; Z = Z, K = K, eq_data...))
end
df_eq = DataFrame(results_Z_K_all_uf)
print(df_eq)
# Get unique values
Z_vals = unique(df_eq.Z)
K_vals = unique(df_eq.K)

# Sort them to be safe
sort!(Z_vals)
sort!(K_vals)

# Create matrix for λmax
λ1_mat = [df_eq[(df_eq.Z .== Z) .& (df_eq.K .== K), :λ1][1] for Z in Z_vals, K in K_vals]

# Plot heatmap
heatmap(Z_vals,K_vals, λ1_mat;
        xlabel = "Consumer constant harvest (Z)",
        ylabel = "K",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

surface(Z_vals, K_vals, λ1_mat;
        xlabel = "Consumer constant harvest (Z)",
        ylabel = "K",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)



#Looking at Y and Z        
results_Z_Y_all_uf = []
for Z in 0.0:0.1:0.75,  Y in 0.0:0.1:0.75
    p = ModelPar_test_1(Z = Z, Y = Y, K = 7.0)
    t = 100
    eq_data = equilibrium_unforced(p, t)
    push!(results_Z_Y_all_uf, (; Z = Z, Y = Y, eq_data...))
end
df_eq = DataFrame(results_Z_Y_all_uf)

# Get unique values
Z_vals = unique(df_eq.Z)
Y_vals = unique(df_eq.Y)

# Sort them to be safe
sort!(Z_vals)
sort!(Y_vals)

# Create matrix for λmax
λ1_mat = [df_eq[(df_eq.Z .== Z) .& (df_eq.Y .== Y), :λ1][1] for Z in Z_vals, Y in Y_vals]

# Plot heatmap
heatmap(Z_vals, Y_vals, λ1_mat;
        xlabel = "Consumer constant harvest (Z)",
        ylabel = "Resource constant harvest (Y)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

surface(Z_vals, Y_vals, λ1_mat;
        xlabel = "Consumer constant harvest (Z)",
        ylabel = "Resource constant harvest (Y)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis) 


##Looking at both K and Z
Ks    = range(3.0, stop=8.0, length=11)
aR_Cs = range(0.1, stop=3.0, length=11)
results_a_K_all_uf = []
for K in Ks, aR_C in aR_Cs
    p = ModelPar_test_1(K = K, aR_C = aR_C)
    t = 100
    eq_data = equilibrium_unforced(p, t)
    push!(results_a_K_all_uf, (;  K ,aR_C, eq_data...))
end
df_eq = DataFrame(results_a_K_all_uf)
print(df_eq)
# Get unique values
a_vals = unique(df_eq.aR_C)
K_vals = unique(df_eq.K)

# Sort them to be safe
sort!(a_vals)
sort!(K_vals)

# Create matrix for λmax
λ1_mat = [df_eq[(df_eq.aR_C .== aR_C) .& (df_eq.K .== K), :λ1][1] for aR_C in a_vals, K in K_vals]

# Plot heatmap
heatmap(a_vals,K_vals, λ1_mat;
        xlabel = "aR_C",
        ylabel = "K",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

surface(a_vals, K_vals, λ1_mat;
        xlabel = "aRC",
        ylabel = "K",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)


#####Plotting Isoclines -- NOT WORKING YET 
# Consumer isocline (vertical line): R = R*
function Rstar(p::ModelPar_test_1)
    num = (p.mC + p.Z)
    den = p.aR_C * (p.e - p.hR_C*(p.mC + p.Z))
    if den <= 0
        return NaN      # no positive R* if energetics violated
    end
    return num / den
end

# Resource isocline: C(R) = [r(1 - R/K) - Y] * (1 + a h R) / a
function C_on_R_isocline(R, p::ModelPar_test_1)
    term = p.r * (1 - R/p.K) - p.Y
    return term * (1 + p.aR_C * p.hR_C * R) / p.aR_C
end

# Effective carrying capacity after Y
Keff(p::ModelPar_test_1) = p.K * (1 - p.Y / p.r)

# ----------------------
# ODEs (for vector field/trajectory)
#
# Plotting helper
# ----------------------
function plot_isoclines(p::ModelPar_test_1; show_vectorfield=true, show_trajectory=true)
    # R grid restricted to feasible region (where C(R) >= 0 and Keff > 0)
    Ke = Keff(p)
    if !(Ke > 0)
        @warn "Keff <= 0 (Y ≥ r). Resource collapses; no positive isocline segment."
        Ke = max(1e-6, p.K)  # still make a small range to avoid errors
    end
    Rmin, Rmax = 1e-6, max(Ke, 1e-3)
    Rgrid = range(Rmin, Rmax; length=400)

    # Resource isocline curve (clip negatives)
    C_R = [max(C_on_R_isocline(R, p), NaN) for R in Rgrid]  # NaN hides negative parts in Plots

    # Consumer isocline (vertical line)
    Rs = Rstar(p)

    # Pick a reasonable C-axis max for plotting
    Cmax_guess = maximum(skipmissing(C_R))
    Cmax = isfinite(Cmax_guess) ? 1.1*Cmax_guess : 1.0

    plt = plot(
        xlabel = "R",
        ylabel = "C",
        legend = :topright,
        framestyle = :box,
        size = (750, 500)
    )

    # Plot R isocline (hump-shaped curve)
    plot!(plt, Rgrid, C_R, lw=2, label="R-isocline:  Ṙ = 0")

    # Plot vertical C-isocline if it exists
    if isfinite(Rs)
        plot!(plt, [Rs, Rs], [0, max(Cmax, 0.1)], lw=2, ls=:dash, label="C-isocline:  Ċ = 0")
    else
        @warn "No positive R* (consumer isocline) — violates e > h(m+Z)"
    end

    # Optional: vector field (arrows)
    if show_vectorfield
        nx, ny = 22, 22
        Rv = range(Rmin, Rmax; length=nx)
        Cv = range(0.0, max(Cmax, 0.1); length=ny)
        RV = repeat(collect(Rv), inner=ny)
        CV = repeat(collect(Cv), outer=nx)

        dR = similar(RV)
        dC = similar(CV)
        for i in eachindex(RV)
            du = zeros(2)
            model_unforced!(du, (RV[i], CV[i]), p, 0.0)
            dR[i] = du[1]
            dC[i] = du[2]
        end
        # normalize arrows for visibility
        len = sqrt.(dR.^2 .+ dC.^2) .+ 1e-12
        uR = 0.06 .* dR ./ len .* (Rmax - Rmin)
        uC = 0.06 .* dC ./ len .* max(Cmax, 0.1)
        quiver!(plt, RV, CV, quiver=(uR, uC), alpha=0.5, label=false)
    end

    # Optional: sample trajectory
    if show_trajectory && isfinite(Rs)
        u0 = [max(Rs*0.7, 1e-3), max(0.5*max(Cmax, 0.2), 1e-6)]
        prob = ODEProblem(model_unforced!, u0, (0.0, 200.0), p)
        sol = solve(prob; reltol=1e-9, abstol=1e-12)
        plot!(plt, sol[1,:], sol[2,:], lw=2, c=:black, label="trajectory")
        scatter!(plt, [sol[1,end]], [sol[2,end]], ms=5, c=:black, label=false)
    end

    return plt
end

# ----------------------
# Run it
# ----------------------
p = ModelPar_test_1()          # tweak parameters here
plt = plot_isoclines(p; show_vectorfield=true, show_trajectory=true)
display(plt)