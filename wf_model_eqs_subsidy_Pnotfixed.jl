##Code to set up equations for wild food model - WITH SUBSIDY, Gutgesell omnivory function, different paramters per trophic level, foraging scale removed to simplify
##Date Initiated: June 18, 2025, added temporal variance in R July 3, 2025
##Contributor(s): Marie K. Gutgesell

#using Pkg
#Pkg.add("Plots")
##Load Libraries needed (ODE for later)
using Parameters: @with_kw, @unpack ##imports Parameters package that provides convenient macros for working with keyword arugemnts, parameter structs and unpacking variables 
using DifferentialEquations
using Plots
using Symbolics
using ForwardDiff
using LinearAlgebra
using NLsolve
using DataFrames
using Interact
using Statistics

##Re-creating model from McCann et al., 2005, Ecology Letters and adding in preference for external subsidy to P 
#u = state variables where u[1] = R1, u[2] = R2, u[3] = C1, u[4] = C2, u[5] = P
#p = Parameters
#t = time

##Next steps:
##work through checks to make sure the math makes sense
##sensitivity analysis across parameter spaces  - does stability change as expected based on existing theory? 
##bifurcation analysis 

##introducing abiotic asynchrony into resources 
##making model stochastic 
##different parameters for r, K, a, e, m, h for different state variables 



##Habitat preference (coupling)
##habitat preference - so preference for habitat changes with shifting densities of C across patches (but preference w itself is fixed, i.e., when prey are at equal densities, this is pref for patch1, influences speed of switching across patches)
function hab_pref(u, p, t)
    R1, R2, C1, C2, P = u
    return p.w * C1/(p.w * C1 + (1 - p.w)*C2)  
end

##simplifying model to remove scale of foraging, so coupling just driven by habitat preference 
##So hab_pref = Si, and so 1-hab_pref is Sj


### Omnivory preference functions -- density dependent preference, where o measures the speed of switching, based on Gutgesell et al., 2022
function om_i_pref_active(u, p, t)
    R1, R2, C1, C2, P = u 
    return (p.o * R1) / (p.o * R1 + (1-p.o) * C1) 
end

function om_j_pref_active(u, p, t)
     R1, R2, C1, C2, P = u 
    return (p.o * R2) / (p.o * R2 + (1-p.o) * C2) 
end

function om_i_pref_fixed(u,p,t)
    return p.o
end

function om_j_pref_fixed(u,p,t)
    return p.o
end

##Grocery Preference -- also density dependent, and H measures speed of switching 
function sub_pref_func(u, p ,t)
    R1, R2, C1, C2, P = u 
    return(p.H * p.G) / (p.H * p.G + (1-p.H)*R1 + (1-p.H)*R2 + (1-p.H) * C1 + (1-p.H)*C2)
end
##is this the right way to have the denominator?  i think so yes
##is this the same as (1-p.H)*(R1+R2+C1+C2)? check math 

##functional responses between resources and consumer 
function f_R1C1(u, p, t)
    @unpack aR_C, hR_C, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    
    return aR_C * R1 / (1 + aR_C * hR_C * R1)
end

function f_R2C2(u, p, t)
    @unpack aR_C, hR_C, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    
    return aR_C * R2 / (1 + aR_C * hR_C * R2)
end



##functional response between resources and predator 
function f_R1P(u, p, t)
    @unpack aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    W1 = p.W(u, p, t) ##function that defines habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    
    numerator = (1-H1) * W1 * aR_P * o1 * R1
    denominator = 1 + (W1 * aR_P * hR_P * o1 * R1 + (1-W1)* aR_P * hR_P * o2 * R2 + W1 * aC_P * hC_P * (1-o1) * C1 + (1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
    return numerator / denominator
end

function f_R2P(u, p, t)
    @unpack aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * (1 - W1) * aR_P * o2 * R2
    denominator = 1 + (W1 * aR_P * hR_P * o1 * R1 + (1-W1)* aR_P * hR_P * o2 * R2 + W1 * aC_P * hC_P * (1-o1) * C1 + (1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
    return numerator / denominator
end

function f_C1P(u, p, t)
    @unpack aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * W1 * aC_P * (1 - o1) * C1
     denominator = 1 + (W1 * aR_P * hR_P * o1 * R1 + (1-W1)* aR_P * hR_P * o2 * R2 + W1 * aC_P * hC_P * (1-o1) * C1 + (1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
   return numerator / denominator
end

function f_C2P(u, p, t)
    @unpack aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * (1 - W1) * aC_P * (1 - o2) * C2
    denominator = 1 + (W1 * aR_P * hR_P * o1 * R1 + (1-W1)* aR_P * hR_P * o2 * R2 + W1 * aC_P * hC_P * (1-o1) * C1 + (1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
    return numerator / denominator
end

##trying out functional resposne for G.. type 2 functional response 
function f_GP(u, p, t)
    @unpack aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = H1 * aG_P * G #trying if i remove scaling / suppression of G by other foraging preferences, i think this makes biological sense (but keep 1-H1 in other FRs)
   denominator = 1 + (W1 * aR_P * hR_P * o1 * R1 + (1-W1)* aR_P * hR_P * o2 * R2 + W1 * aC_P * hC_P * (1-o1) * C1 + (1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
    return numerator / denominator
end


##linear (type 1) functional response 
function f_GP_2(u, p, t)
    @unpack a, h, H, G = p 
    R1, R2, C1, C2 = u 
    return G * H
end



#Set up paramters - starting with just looking at effect of coupling, no subsidy and no omnivory
@with_kw mutable struct ModelPar_passive
    w = 0.5   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
    H = 0.5 ##H = preference for groceries (G)

    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 2.0
    K = 3.0
    aR_C = 1.0  ##attack rate  of consumer on R
    aR_P = 0.2 ##attack rate of P on R 
    aC_P = 0.5  ##attack rate of P on C 
    aG_P = 0.5  ##attack rate of P on G
    e = 0.5   ##energy conversion 
    mC = 0.4 ##C mortality rate
    mP = 0.2   ## P mortality rate
    hR_C = 0.4 ##handling time of C on R
    hR_P = 0.6  ##handling time of P on R  
    hC_P = 0.3 ##handling time of P on C
    hG_P = 0.3 ##handling time of P on G

    ##Density dependent habitat preference function (simplifying for now to remove foraging scale)
    W::Function = hab_pref
    
    ##Omnivory preference function
    d_om_i::Function = om_i_pref_fixed ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = om_j_pref_fixed ##degree of omnivory in patch 2 
    
    ##Subsidy preference function
    sub_pref::Function = sub_pref_func

    ##Initial/constant value of groceries, so effectively instantly replenishes 
    G = 2.0
    ##Predator functional responses
    f_r1c1::Function = f_R1C1
    f_r2c2::Function = f_R2C2
    f_r1p::Function = f_R1P
    f_r2p::Function = f_R2P
    f_c1p::Function = f_C1P
    f_c2p::Function = f_C2P
    f_gp::Function = f_GP
    f_gp2::Function = f_GP_2

      ##temporal variation in R parameters
   l1 = 0.1  #magnitude of variation in K of R1 
   l2 = 0.1  #magnitude of variation in K of R2
   pf = 10.0 ##period of fluctuation 
   D = 0.5 ##phase delay between K1 and K2 (0.5 = perfectly asynchronous)
    e1::Function = t -> sin(2π / pf * t)
    e2::Function = t -> sin(2π / pf * (t - D *pf))
end

#Set up paramters - starting with just looking at effect of coupling, no subsidy and no omnivory
@with_kw mutable struct ModelPar_active
    w = 0.7   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.8 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
    H = 0.2 ##H = preference for groceries (G)

    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 2.0
    K = 3.0
    aR_C = 1.0  ##attack rate  of consumer on R
    aR_P = 0.2 ##attack rate of P on R 
    aC_P = 0.5  ##attack rate of P on C 
    aG_P = 0.5  ##attack rate of P on G
    e = 0.5   ##energy conversion 
    mC = 0.4 ##C mortality rate
    mP = 0.2   ## P mortality rate
    hR_C = 0.4 ##handling time of C on R
    hR_P = 0.6  ##handling time of P on R  
    hC_P = 0.3 ##handling time of P on C
    hG_P = 0.3 ##handling time of P on G

    ##Density dependent habitat preference function (simplifying for now to remove foraging scale)
    W::Function = hab_pref
    
    ##Omnivory preference function
    d_om_i::Function = om_i_pref_active ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = om_j_pref_active ##degree of omnivory in patch 2 
    
    ##Subsidy preference function
    sub_pref::Function = sub_pref_func

    ##Initial/constant value of groceries, so effectively instantly replenishes 
    G = 2.0
    ##Predator functional responses
    f_r1c1::Function = f_R1C1
    f_r2c2::Function = f_R2C2
    f_r1p::Function = f_R1P
    f_r2p::Function = f_R2P
    f_c1p::Function = f_C1P
    f_c2p::Function = f_C2P
    f_gp::Function = f_GP
    f_gp2::Function = f_GP_2

      ##temporal variation in R parameters
   l1 = 0.1  #magnitude of variation in K of R1 
   l2 = 0.1  #magnitude of variation in K of R2
   pf = 10.0 ##period of fluctuation 
   D = 0.5 ##phase delay between K1 and K2 (0.5 = perfectly asynchronous)
    e1::Function = t -> sin(2π / pf * t)
    e2::Function = t -> sin(2π / pf * (t - D *pf))
end

##Model 
function model_forced!(du, u, p ,t)
    @unpack r, K, aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G, e, mC, mP, H, l1, l2, e1, e2 = p
   R1, R2, C1, C2, P = u 
 
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
    f_gp =p.f_gp(u, p, t)
   # f_gp2 = p.f_gp2(u, p, t)
## consumer functional responses
  f_r1c1 = p.f_r1c1(u, p, t)
  f_r2c2 = p.f_r2c2(u, p, t)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / (K - l1*(e1(t) - 0.5))) - C1 * f_r1c1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / (K - l2*(e2(t) - 0.5))) - C2 * f_r2c2 - P * f_r2p
   du[3] = e * C1 * f_r1c1 - P * f_c1p - mC * C1
   du[4] = e * C2 * f_r2c2 - P * f_c2p - mC * C2 
   du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp - mP * P
  #  du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 

function model_unforced!(du, u, p ,t)
    @unpack r, K, aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G, e, mC, mP, H, l1, l2, e1, e2 = p
   R1, R2, C1, C2, P = u 
 
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
    f_gp =p.f_gp(u, p, t)
   # f_gp2 = p.f_gp2(u, p, t)
## consumer functional responses
  f_r1c1 = p.f_r1c1(u, p, t)
  f_r2c2 = p.f_r2c2(u, p, t)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_r1c1 - P * f_r1p ##if want to use equations w/o temporal forcing
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_r2c2 - P * f_r2p
   du[3] = e * C1 * f_r1c1 - P * f_c1p - mC * C1
   du[4] = e * C2 * f_r2c2 - P * f_c2p - mC * C2 
   du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp - mP * P
  #  du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 


##function to calculate total harvest for P, by summing functional response for P at each time step
function total_FR_into_P(u, p ,t)
    @unpack r, K, aR_P, aC_P,aG_P, hR_P, hC_P, hG_P, G, e, mC, mP, H = p
   R1, R2, C1, C2, P = u 
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
   f_gp = p.f_gp(u, p, t) 
   # f_gp2 = p.f_gp2(u, p, t)
    total_FR = f_r1p + f_r2p + f_c1p + f_c2p + f_gp
   return (total_FR, f_r1p, f_r2p, f_c1p, f_c2p, f_gp)
end 


# Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
function rhs_forced(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_forced!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

function rhs_unforced(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_unforced!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end


##other functions for eigenvalue analysis - based on KC code
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
function equilibrium_forced(p)
    u0 = [1.5, 1.5, 1.0, 1.0, 0.5] ##initial condition
    
    t_warmup = 300.0 ##run long enough to reach the limit cycle 
    t_eval = 500.0 ##window to evaluate system properties
    tspan = (0.0, t_eval)
    
    t_grid = range(t_warmup, t_eval, length = 1000) ##extract dynamics after settling 

    #extract solutions over limit cycle 
    prob = ODEProblem(model_forced!, u0, tspan, deepcopy(p)) ##only need to use deepcopy if you are changing parameters inside the function
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    sol_grid = sol(t_grid)
    
    ##Calculate mean state metrics over the limit cycle 
    mean_state = mean(sol_grid, dims = 2)
    range_state = maximum(sol_grid, dims = 2) - minimum(sol_grid, dims = 2)
    min_state = minimum(sol_grid, dims =2)

    ##Calculate CV for one period of K
    pf = p.pf
    t_cv = range(t_eval - pf, t_eval, length = 1000)
    sol_cv = sol(t_cv)

    mean_cv = mean(sol_cv, dims = 2)
    sd_cv = std(sol_cv, dims = 2)
    cv = sd_cv ./ mean_cv
    
    ##want to add loop here so can calculate eigenvalue at each point t over whole period, and integrate to get total eigenvalue (see Bieg et al., 2023)
   # Eigenvalue tracking over one period
t_eig = range(t_eval - pf, t_eval, length = 500)
λ_max_vals = zeros(length(t_eig))

for (i, t) in enumerate(t_eig)
    x = sol(t)
    J = ForwardDiff.jacobian(x -> rhs_forced(x, p, t), x)
    λ_max_vals[i] = maximum(real.(eigvals(J)))
end

# Integrate λ_max over the period (e.g., trapezoidal rule)
dt = step(t_eig)
λ_integrated = sum(λ_max_vals[2:end-1]) * dt + 0.5 * dt * (λ_max_vals[1] + λ_max_vals[end])



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
        cv = vec(cv),
        λ1 = λ1,
        λ1_imag = λ1_imag,
         λ_integrated = λ_integrated,
        react = react,
        sol = sol)
end 


function equilibrium_unforced(p, t)
    u0 = [1.5, 1.5, 1.0, 1.0, 0.5] ##initial condition
    tspan = (0.0, t)
    
    #simulate dynamics to approach equilibrium
    prob = ODEProblem(model_unforced!, u0, tspan, deepcopy(p)) ##only need to use deepcopy if you are changing parameters inside the function
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
   
   #use ODE result as initial guess for equilibrium
    u_approx = sol(t) ##returns full vector of state variables at time t
   eq = nlsolve((du, u) -> model_unforced!(du, u, deepcopy(p), 0.0), u_approx).zero
 cmat(u, p) = ForwardDiff.jacobian(x -> rhs_unforced(x, p), u)

    ##Compute the community matrix and stability metrics 
    M = cmat(eq, p)
    λ1 = λ1_stability(M)
    λ1_imag = λ1_stability_imag(M)
    react = ν_stability(M)

    return(eq = eq, λ1 = λ1, λ1_imag =λ1_imag, react = react)
end 




##STRUCTURE 1: Plotting dynamics, equilibrium, eigenvalue analysis 
##Solve ODE 
##set initial condition
u0 = [1.5, 1.5, 1.0, 1.0, 0.5]
tspan = (0.0, 1000.0)
##set Parameters
p = ModelPar_passive(w = 0.2, o = 0.3, H = 0.0, K = 3.0)

##Define the ODE problem
prob_1 = ODEProblem(rhs_forced, u0, tspan, p)
sol_1 = solve(prob_1)

##plot timeseries
plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution - Passive Omnivory")

for (u, t) in zip(sol_1.u, sol_1.t)
    W = hab_pref(u, p, t)
    o = om_i_pref_fixed(u, p, t)
    H = sub_pref_func(u, p, t)
    println("t=$t, W=$W, o=$o, H=$H")
end

##Look at predator consumption  
fr_vals_1 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_1.u, sol_1.t)]

# Separate each series into its own array
fr_total = [x[1] for x in fr_vals_1]
fr_R1 = [x[2] for x in fr_vals_1]
fr_R2 = [x[3] for x in fr_vals_1]
fr_C1 = [x[4] for x in fr_vals_1]
fr_C2 = [x[5] for x in fr_vals_1]
fr_G = [x[6] for x in fr_vals_1]
times = sol_1.t

   ##plot total P consumption 
plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Predator Consumption")  
plot!(times, fr_R1, label = "R1 → P")
plot!(times, fr_R2, label = "R2 → P")
plot!(times, fr_C1, label = "C1 → P")
plot!(times, fr_C2, label = "C2 → P")
plot!(times, fr_G, label = "G → P")



##Look at dynamics with active omnivory 
##set Parameters
u0 = [0.6, 0.8, 0.45, 0.65, 0.2]
tspan = (0.0, 500.0)
p = ModelPar_active(w = 0.2, o = 0.2, H = 1.0)
##Define the ODE problem
prob_2 = ODEProblem(rhs_forced, u0, tspan, p)
sol_2 = solve(prob_2)

##plot timeseries
plot(sol_2, xlabel="Time", ylabel="Population", title="ODE Solution - Active Omnivory")
  
##Look at predator consumption  
fr_vals_2 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)]

# Separate each series into its own array
fr_total = [x[1] for x in fr_vals_2]
fr_R1 = [x[2] for x in fr_vals_2]
fr_R2 = [x[3] for x in fr_vals_2]
fr_C1 = [x[4] for x in fr_vals_2]
fr_C2 = [x[5] for x in fr_vals_2]
fr_G = [x[6] for x in fr_vals_2]
times = sol_2.t

   ##plot total P consumption 
plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Predator Consumption")  
plot!(times, fr_R1, label = "R1 → P")
plot!(times, fr_R2, label = "R2 → P")
plot!(times, fr_C1, label = "C1 → P")
plot!(times, fr_C2, label = "C2 → P")
plot!(times, fr_G, label = "G → P")


 