##Code to work out how to generate asynchrony in R1 and R2 
##Initiated: July 3, 2025
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

##Just using simple C, R1, R2 model
##habitat preference - so preference for habitat changes with shifting densities of C across patches (but preference w itself is fixed, i.e., when prey are at equal densities, this is pref for patch1, influences speed of switching across patches)
function hab_pref(u, p, t)
    R1, R2, C = u
    return p.w * R1/(p.w * R1 + (1 - p.w)*R2)  
end

##functional responses between resources and consumer 
function f_R1C(u, p, t)
    @unpack aR_C, hR_C= p  
    R1, R2, C = u
    W1 = p.W(u, p, t)

    return W1 * aR_C * R1 / (1 + aR_C * hR_C * R1)
end

function f_R2C(u, p, t)
    @unpack aR_C, hR_C = p  
    R1, R2, C= u  
    W1 = p.W(u, p, t)
    return (1- W1) * aR_C * R2 / (1 + aR_C * hR_C * R2)
end


#Set up paramters 
@with_kw mutable struct ModelPar_test
  
    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 2.0
    K = 5.0
    aR_C = 0.8  ##attack rate  of consumer on R
    e = 0.7   ##energy conversion 
    mC = 0.5 ##C mortality rate
    hR_C = 1.0 ##handling time of C on R
   
    #habitat preference: 
    w = 0.2
    ##Density dependent habitat preference function (simplifying for now to remove foraging scale)
    W::Function = hab_pref
    
   ##temporal variation in R parameters
   l1 = 0.1  #magnitude of variation in K of R1 
   l2 = 0.1  #magnitude of variation in K of R2
   pf = 5.0 ##period of fluctuation 
   D = 0.5 ##phase delay between K1 and K2 (0.5 = perfectly asynchronous)
    e1::Function = t -> sin(2π / pf * t)
    e2::Function = t -> sin(2π / pf * (t - D *pf))
    
    ##C functional responses
    f_r1c::Function = f_R1C
    f_r2c::Function = f_R2C
    
end


##Model 
function model_forced!(du, u, p ,t)
    @unpack r, K,aR_C, e, mC, hR_C, l1, l2, e1, e2 = p
   R1, R2, C = u 

## consumer functional responses
  f_r1c = p.f_r1c(u, p, t)
  f_r2c = p.f_r2c(u, p, t)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / (K - l1*(e1(t) - 0.5))) - C * f_r1c
   du[2] = r * R2 * (1 - R2 / (K - l2*(e2(t) - 0.5)))  - C * f_r2c
   du[3] = e * C * (f_r1c + f_r2c)  - mC * C
 
   return du
 end 

 # Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_forced(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_forced!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end


##Solve ODE 
##set initial condition
u0 = [1.0, 1.0, 1.5]
tspan = (0.0, 500.0)
##set Parameters
p = ModelPar_test(K = 5.0)

##Define the ODE problem
prob_1 = ODEProblem(rhs_forced, u0, tspan, p)
sol_1 = solve(prob_1)

##plot timeseries
plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution")

##trying my own equilibrium function, that follows structure similar to KC 
##other functions for eigenvalue analysis - based on KC code
find_eq(u, p) = nlsolve((du, u) -> model_forced!(du, u, p, zero(u)), u).zero
cmat(u, p) = ForwardDiff.jacobian(x -> rhs_forced(x, p), u)


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
function equilibrium_forced(p)
    u0 = [1.0,1.0, 1.5] ##initial condition
    
    t_warmup = 300.0 ##run long enough to reach the limit cycle 
    t_eval = 500.0 ##window to evaluate system properties
    tspan = (0.0, t_eval)
    
    t_grid = range(t_warmup, t_eval, length = 1000) ##extract dynamics after settling 

    #extract solutions over limit cycle 
    prob = ODEProblem(model_forced!, u0, tspan, deepcopy(p)) ##only need to use deepcopy if you are changing parameters inside the function
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    sol_grid = sol(t_grid)

   #use ODE result as initial guess for equilibrium
    #u_approx = sol(t) ##returns full vector of state variables at time t
   # eq = nlsolve((du, u) -> model_forced!(du, u, deepcopy(p), t), u_approx).zero ##there is no fixed equilibrium, so don't use nlsolve 
    #eq = sol(t)

    ##Calculate mean state metrics over the limit cycle 
    mean_state = mean(sol_grid, dims = 2)
    range_state = maximum(sol_grid, dims = 2) - minimum(sol_grid, dims = 2)
    min_state = minimum(sol_grid, dims =2)
    ##Calculate jacobian at final state
    x_eval = sol(t_eval) ##extract solution at time = 500, when want to evaluate system
    J = ForwardDiff.jacobian(x -> rhs_forced(x, p, t_eval), x_eval)

    #M = cmat(eq, p)
    λ1 = λ1_stability(J)
    λ1_imag = λ1_stability_imag(J)
    react = ν_stability(J)

    return( mean_state = vec(mean_state),
        amplitude = vec(range_state),
        min = vec(min_state),
        λ1 = λ1,
        λ1_imag = λ1_imag,
        react = react,
        sol = sol)
end 


p = ModelPar_test()

K_results = equilibrium_forced(p)

results_K_all = []
for K in 0.1:0.1:10.0
    p = ModelPar_test(K=K)
    eq_data = equilibrium_forced(p)
    push!(results_K_all, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all)

plot(df_eq.K, df_eq.λ1)

R1_amp = [row.amplitude[1] for row in eachrow(df_eq)]
R2_amp = [row.amplitude[2] for row in eachrow(df_eq)]
C_amp = [row.amplitude[3] for row in eachrow(df_eq)]

plot(df_eq.K, R1_amp)
plot(df_eq.K, R2_amp)
plot(df_eq.K, C_amp)


results_w_all = []
for w in 0.5:0.1:1.0
    p = ModelPar_test(w=w)
    eq_data = equilibrium_forced(p)
    push!(results_w_all, (; w=w, eq_data...))
end
df_eq_w = DataFrame(results_w_all)

plot(df_eq_w.w, df_eq_w.λ1)


C_min = [row.min[3] for row in eachrow(df_eq_w)]
plot(df_eq_w.w, C_min)

##okay great, so trends are working the way i would expect based on McCann et al., 2005, figure 3b 

# Extract time and solution values
ts = sol_1.t
R1_vals = sol_1[1, :]
R2_vals = sol_1[2, :]


# Filter for t between 1 and 100
idx = findall(t -> t ≥ 1 && t ≤ 25, ts)

# Subset data
ts_sub = ts[idx]
R1_sub = R1_vals[idx]
R2_sub = R2_vals[idx]

# Plot
plot(ts_sub, R1_sub, label = "R1", xlabel = "Time", ylabel = "Abundance", lw = 2)
plot!(ts_sub, R2_sub, label = "R2", lw = 2, linestyle = :dash)
# Plot
plot(ts, R1_vals, label = "R1", xlabel = "Time", ylabel = "Abundance", lw = 2)
plot!(ts, R2_vals, label = "R2", lw = 2, linestyle = :dash)


##Looking into dynamics of just K for R1 and R2 to see if asynchrnous dynamics are working
K = 3.0
l1 = 0.5
l2 = 0.5
D = 0.5
pf = 1.0 ##period of fluctuation
t = 0:0.1:25
e1 = @. sin(2π / pf * t)
e2 = @. sin(2π / pf * (t - D *pf))


Keff_1 = @. K - l1 * (e1 - 0.5)  # assuming e1(t) = sin(...) form
Keff_2 = @. K - l2 * (e2 - 0.5)


plot(t, Keff_1, label = "Effective K for R1", lw = 2, xlabel = "Time", ylabel = "Keff", ylim = (1.5, 4.5))
plot!(t, Keff_2, label = "Effective K for R2", lw = 2, linestyle = :dash)




##Comparing stability responses to unforced model
##Model 
@with_kw mutable struct ModelPar_test_unforced
  
    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 3.0
    K = 1.0
    aR_C = 2.0  ##attack rate  of consumer on R
    e = 0.7   ##energy conversion 
    mC = 0.2 ##C mortality rate
    hR_C = 0.4 ##handling time of C on R
   
     ##C functional responses
    f_r1c::Function = f_R1C
    f_r2c::Function = f_R2C
    
end

function model_unforced!(du, u, p ,t)
    @unpack r, K,aR_C, e, mC, hR_C = p
   R1, C = u 

## consumer functional responses
  f_r1c = p.f_r1c(u, p, t)
 # f_r2c = p.f_r2c(u, p, t)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 /K) - C * f_r1c
  # du[2] = r * R2 * (1 - R2 /K)  - C * f_r2c
   du[2] = e * C * (f_r1c)  - mC * C
 
   return du
 end 

 # Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_unforced(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_unforced!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

##Solve ODE 
##set initial condition
u0 = [1.0, 1.5]
tspan = (0.0, 100.0)
##set Parameters
p = ModelPar_test_unforced(K = 1.3)

##Define the ODE problem
prob_2 = ODEProblem(rhs_unforced, u0, tspan, p)
sol_2 = solve(prob_2)

##plot timeseries
plot(sol_2, xlabel="Time", ylabel="Population", title="ODE Solution")


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
function equilibrium_unforced(p, t)
    u0 = [1.0, 1.5] ##initial condition
    tspan = (0.0, t)
    
    #simulate dynamics to approach equilibrium
    prob = ODEProblem(model_unforced!, u0, tspan, deepcopy(p)) ##only need to use deepcopy if you are changing parameters inside the function
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
   
   #use ODE result as initial guess for equilibrium
    u_approx = sol(t) ##returns full vector of state variables at time t
    eq = nlsolve((du, u) -> model_unforced!(du, u, deepcopy(p), 0.0), u_approx).zero

    ##Compute the community matrix and stability metrics 
    M = cmat(eq, p)
    λ1 = λ1_stability(M)
    λ1_imag = λ1_stability_imag(M)
    react = ν_stability(M)

    return(eq = eq, λ1 = λ1, λ1_imag =λ1_imag, react = react)
end 

p = ModelPar_test_unforced()
t = (100.00)
K_results = equilibrium_unforced(p, t)

results_K_all = []
for K in 0.1:0.1:2.0
    p = ModelPar_test_unforced(K=K)
    t = (100.0)
    eq_data = equilibrium_unforced(p, t)
    push!(results_K_all, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all)

plot(df_eq.K, df_eq.λ1)


test