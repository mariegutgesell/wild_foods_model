##Code to set up equations for wild food model - WITH SUBSIDY, Gutgesell omnivory function, different paramters per trophic level, foraging scale removed to simplify
##Date Initiated: June 18, 2025, added temporal variance in R July 3, 2025
##Contributor(s): Marie K. Gutgesell

#using Pkg
#Pkg.add("Plots")
#Pkg.add("QuadGK")
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
using QuadGK: quadgk
##Re-creating model from McCann et al., 2005, Ecology Letters and adding in preference for external subsidy to P 
#u = state variables where u[1] = R1, u[2] = R2, u[3] = C1, u[4] = C2, u[5] = P
#p = Parameters
#t = time

##Next steps:
##work through checks to make sure the math makes sense
##sensitivity analysis across parameter spaces  - does stability change as expected based on existing theory? 
##bifurcation analysis 
##making model stochastic 


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
    G = p.G_func(t)

    return(p.H * p.G) / (p.H * p.G + (1-p.H)*(R1 +R2 + C1 +C2))
end
##is this the right way to have the denominator?  i think so yes
##is this the same as (1-p.H)*(R1+R2+C1+C2)? check math 

##functional responses between resources and consumer 
function f_R1C1(u, p, t)
    @unpack aR_C, hR_C = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    G = p.G_func(t)
    return aR_C * R1 / (1 + aR_C * hR_C * R1)
end

function f_R2C2(u, p, t)
    @unpack aR_C, hR_C = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    G = p.G_func(t)

    return aR_C * R2 / (1 + aR_C * hR_C * R2)
end



##functional response between resources and predator 
function f_R1P(u, p, t)
    @unpack aR_P,aC_P, aG_P, hR_P, hC_P, hG_P = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2, P = u  ##defines state variables, G = groceries
    W1 = p.W(u, p, t) ##function that defines habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)
    
    numerator = (1-H1) * W1 * aR_P * o1 * R1
    denominator = 1 + ((1-H1)*W1 * aR_P * hR_P * o1 * R1 + (1-H1)*(1-W1)* aR_P * hR_P * o2 * R2 + (1-H1)*W1 * aC_P * hC_P * (1-o1) * C1 + (1-H1)*(1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
    return numerator / denominator
end

function f_R2P(u, p, t)
   @unpack aR_P,aC_P, aG_P, hR_P, hC_P, hG_P = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * (1 - W1) * aR_P * o2 * R2
   denominator = 1 + ((1-H1)*W1 * aR_P * hR_P * o1 * R1 + (1-H1)*(1-W1)* aR_P * hR_P * o2 * R2 + (1-H1)*W1 * aC_P * hC_P * (1-o1) * C1 + (1-H1)*(1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
        return numerator / denominator
end

function f_C1P(u, p, t)
   @unpack aR_P,aC_P, aG_P, hR_P, hC_P, hG_P = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * W1 * aC_P * (1 - o1) * C1
  denominator = 1 + ((1-H1)*W1 * aR_P * hR_P * o1 * R1 + (1-H1)*(1-W1)* aR_P * hR_P * o2 * R2 + (1-H1)*W1 * aC_P * hC_P * (1-o1) * C1 + (1-H1)*(1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
      return numerator / denominator
end

function f_C2P(u, p, t)
   @unpack aR_P,aC_P, aG_P, hR_P, hC_P, hG_P = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * (1 - W1) * aC_P * (1 - o2) * C2
    denominator = 1 + ((1-H1)*W1 * aR_P * hR_P * o1 * R1 + (1-H1)*(1-W1)* aR_P * hR_P * o2 * R2 + (1-H1)*W1 * aC_P * hC_P * (1-o1) * C1 + (1-H1)*(1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
        return numerator / denominator
end

##trying out functional resposne for G.. type 2 functional response 
function f_GP(u, p, t)
   @unpack aR_P,aC_P, aG_P, hR_P, hC_P, hG_P = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2, P = u  ##defines state variables
    W1 = p.W(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = H1 * aG_P * G #trying if i remove scaling / suppression of G by other foraging preferences, i think this makes biological sense (but keep 1-H1 in other FRs)
   denominator = 1 + ((1-H1)*W1 * aR_P * hR_P * o1 * R1 + (1-H1)*(1-W1)* aR_P * hR_P * o2 * R2 + (1-H1)*W1 * aC_P * hC_P * (1-o1) * C1 + (1-H1)*(1-W1) * aC_P * hC_P * (1-o2) * C2 + H1 * aG_P * hG_P * G)
        return numerator / denominator
end


##linear (type 1) functional response 
function f_GP_2(u, p, t)
    @unpack a, h, H = p 
    R1, R2, C1, C2 = u 
    G = p.G_func(t)
    return G * H
end

##set function so before pulse G = G_pre, and during pulse is G_pulse 
G_func(t) = (t < t_pulse || t ≥ t_recover) ? G_pre : G_pulse



#Set up paramters - starting with just looking at effect of coupling, no subsidy and no omnivory
@with_kw mutable struct ModelPar_passive
    w = 0.5   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
    H = 0.5 ##H = preference for groceries (G)

    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level 
    r = 1.0
    K = 3.25
    aR_C = 2.5  ##attack rate  of consumer on R
    aR_P = 2.5 ##attack rate of P on R 
    aC_P = 2.5  ##attack rate of P on C 
    aG_P = 2.5  ##attack rate of P on G
    e = 0.8   ##energy conversion 
    mC = 1.0 ##C mortality rate
    mP = 1.0   ## P mortality rate
    hR_C = 0.5 ##handling time of C on R
    hR_P = 0.5  ##handling time of P on R  
    hC_P = 0.5 ##handling time of P on C
    hG_P = 0.5 ##handling time of P on G

    ##Density dependent habitat preference function (simplifying for now to remove foraging scale)
    W::Function = hab_pref
    
    ##Omnivory preference function
    d_om_i::Function = om_i_pref_fixed ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = om_j_pref_fixed ##degree of omnivory in patch 2 
    
    ##Subsidy preference function
    sub_pref::Function = sub_pref_func

    ##Initial/constant value of groceries, so effectively instantly replenishes 
    G = 2.0
    G_base = 2.0 ##need base for pulse perturbation experiment
   
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

    ##time varying function of G
     G_func::Function = t -> 1.0

end

#Set up paramters - starting with just looking at effect of coupling, no subsidy and no omnivory
@with_kw mutable struct ModelPar_active
    w = 0.2   ##w = habitat preference for patch i 
    o = 0.2 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
    H = 0.1 ##H = preference for groceries (G)

    ##Model parameters, for now just keeping these parameters the same for each patch, different per trophic level - based on values in Fig 3, McCann et al., 2005
    r = 1.0
    K = 3.25
    aR_C = 2.5  ##attack rate  of consumer on R
    aR_P = 4.0 ##attack rate of P on R 
    aC_P = 3.4  ##attack rate of P on C 
    aG_P = 3.4  ##attack rate of P on G
    e = 0.8   ##energy conversion -
    mC = 1.0 ##C mortality rate
    mP = 0.45   ## P mortality rate
    hR_C = 0.4 ##handling time of C on R
    hR_P = 1.25  ##handling time of P on R  
    hC_P = 1.25 ##handling time of P on C
    hG_P = 1.25 ##handling time of P on G

    ##Density dependent habitat preference function (simplifying for now to remove foraging scale)
    W::Function = hab_pref
    
    ##Omnivory preference function
    d_om_i::Function = om_i_pref_active ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = om_j_pref_active ##degree of omnivory in patch 2 
    
    ##Subsidy preference function
    sub_pref::Function = sub_pref_func

    ##Initial/constant value of groceries, so effectively instantly replenishes 
    G = 2.0
    G_base = 2.0


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
   #l1 = 1.0  #magnitude of variation in K of R1 
   #l2 = 1.0  #magnitude of variation in K of R2
   l = 1.0
   pf = 10.0 ##period of fluctuation 
   D = 0.5 ##phase delay between K1 and K2 (0.5 = perfectly asynchronous)
    e1::Function = t -> sin(2π / pf * t)
    e2::Function = t -> sin(2π / pf * (t - D *pf))

    ##time varying function of G
     G_func::Function = t -> 1.0

     ##fixed P value
    # P_fixed = 0.25
end

##Model 
function model_forced!(du, u, p ,t)
    @unpack r, K,  aR_P, aC_P, aG_P, hR_P, hC_P, hG_P, e,  mC, H, l, e1, e2 = p
   R1, R2, C1, C2, P = u 
   G = p.G_func(t)
 
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
  # du[1] = r * R1 * (1 - R1 / (K - l*(e1(t) - 0.5))) - C1 * f_r1c1 - P * f_r1p ##subtracting 0.5 from e1 essentially increases mean K (Kmean = K + 0.5l1)
  # du[2] = r * R2 * (1 - R2 / (K - l*(e2(t) - 0.5))) - C2 * f_r2c2 - P * f_r2p
    du[1] = r * R1 * (1 - R1 / (K - l*(e1(t)))) - C1 * f_r1c1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / (K - l*(e2(t)))) - C2 * f_r2c2 - P * f_r2p
   du[3] = e * C1 * f_r1c1 - P * f_c1p - mC * C1
   du[4] = e * C2 * f_r2c2 - P * f_c2p - mC * C2 
   #du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp - mP * P
    du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 

 ##can adjust model and do forced/unforced just by changing parameters l1 andl2, will simplify code so dont need to do all forced/unforced i think.... just set up equilibirum functions to do both, maybe have if statements 
function model_unforced!(du, u, p ,t)
   @unpack r, K,  aR_P, aC_P, aG_P, hR_P, hC_P, hG_P, e,  mC, H = p
    R1, R2, C1, C2, P = u 
     G = p.G_func(t)
 
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
  # du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp - mP * P
   du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 


##function to calculate total harvest for P, by summing functional response for P at each time step
function total_FR_into_P(u, p ,t)
  @unpack r, K,  aR_P, aC_P, aG_P, hR_P, hC_P, hG_P, e, mC, H = p
    R1, R2, C1, C2, P = u 
   G = p.G_func(t)
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

##Note: next if using NLSolve need to do wrapper function - solving for equilibrium only in a subset of the state variables
function wrapped_model_forced!(residuals, x, p, P_fixed)
    R1, R2, C1, C2 = x
    T = eltype(x)
    P = convert(T, P_fixed)
    u = T[R1, R2, C1, C2, P]  # use T[...] to create uniform type array

    du = zeros(T,5)
    model_forced!(du, u, p, 0.0)

    residuals[1:4] .= du[1:4]
    return residuals
end
##vector-returning version, to allow for jacobian calculation
function wrapped_model_vec_forced(x, p, P_fixed)
    T = eltype(x)
    du = zeros(T, 4)
    wrapped_model_forced!(du, x, p, P_fixed)
    return du
end

##Note: next if using NLSolve need to do wrapper function 
function wrapped_model_unforced!(residuals, x, p, P_fixed)
    R1, R2, C1, C2 = x
    T = eltype(x)
    P = convert(T, P_fixed)
    u = T[R1, R2, C1, C2, P]  # use T[...] to create uniform type array

    du = zeros(T,5)
    model_unforced!(du, u, p, 0.0)

    residuals[1:4] .= du[1:4]
    return residuals
end
##vector-returning version, to allow for jacobian calculation
function wrapped_model_vec_unforced(x, p, P_fixed)
    T = eltype(x)
    du = zeros(T, 4)
    wrapped_model_unforced!(du, x, p, P_fixed)
    return du
end

#The wrapped functions (e.g., wrapped_model_unforced!) are only needed when using solvers like nlsolve that require:
#   Equal-length input/output vectors (i.e., length(x) == length(f(x))).
 #   A lower-dimensional system (in your case, solving only for [R1, R2, C1, C2] with P held constant).
  #  Residual format rather than time-dependent ODEs.
#In contrast, ODEProblem:
 #   Solves full time-dependent systems over all state variables.
  #  Accepts the full system (i.e., model_unforced!) and will handle all five state variables just fine.

function find_eq_forced(u, p)
    x0 = u[1:4]
    P_fixed = u[5]

    result = nlsolve((du, x) -> wrapped_model_forced!(du, x, deepcopy(p), P_fixed), x0)

    if result.f_converged
        return vcat(result.zero, P_fixed)
    else
        error("Equilibrium solver did not converge.")
    end
end


cmat_forced(u4, p, Pstar) = ForwardDiff.jacobian(x -> wrapped_model_vec_forced(x, p, Pstar), u4)

##other functions for eigenvalue analysis - based on KC code
"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
#λ1_stability(M) = maximum(real.(eigvals(M)))

function λ1_stability(M)
    if any(isnan.(M)) || any(isinf.(M))
        return NaN
    end
    λs = eigvals(M)
    return maximum(real.(λs))
end

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
function λ1_stability_imag(M) 
    ev = eigvals(M)
    imag.(ev)[findall(real.(ev) .== maximum(real.(ev)))[1]]
end 

"""M is the community matrix, we can be calculated with `cmat(u, p)`
Note: `\nu` is the what to input `ν` which looks a bit too much like `v` for my taste
"""
ν_stability(M) = λ1_stability((M + M') / 2)

# overshoot and oscillation range

abs_sol(sol, t, eq) = abs.(sol(t) .- eq)

function overshoot(sol, eq, spc, t_beg, t_end)
    return quadgk(t -> abs_sol(sol, t, eq)[spc], t_beg, t_end)[1]
end


# Calculate max-min metric
function min_max(sol, spc, t_beg, t_end, len = 100000)
    return maximum(sol(range(t_beg, t_end, length = len))[spc, :]) - 
    minimum(sol(range(t_beg, t_end, length = len))[spc, :])
end

# find first time equilibrium is hit
function find_times_hit_equil_press(res)
    eq = res[1, end], res[2, end], res[3, end], res[4, end]
    times = zeros(4)
    for spc in 1:4
        for i in 20:length(res)
            # cannot be too strict here otherwise the value of the 
            # first ht time varies a lort which will have serious 
            # impact on min and max (overshoot too) leading to major oscillations lenth of the ts must be high enough too.
            if isapprox(res[spc, i], eq[spc], atol = 0.01)
                times[spc] = res.t[i]
                break
            end
        end
    end
    return times
end


##trying function to calculate equilibrium that i can then loop over for values of K etc. -using KC approaches
function equilibrium_forced(p, P0)
    u0 = [1.5, 1.5, 1.0, 1.0, P0] ##initial condition
    
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
    sd_state = std(sol_grid, dims = 2; corrected = false)
    range_state = maximum(sol_grid, dims = 2) - minimum(sol_grid, dims = 2)
    min_state = minimum(sol_grid, dims =2)
    max_state = maximum(sol_grid, dims = 2)
    cv = sd_state ./mean_state
   
    ##want to add loop here so can calculate eigenvalue at each point t over whole period, and integrate to get total eigenvalue (see Bieg et al., 2023)
   # Eigenvalue tracking over one period
   pf = p.pf
t_eig = range(t_eval - pf, t_eval, length = 500)
λ_max_vals = zeros(length(t_eig))

for (i, t) in enumerate(t_eig)
    x = sol(t)[1:4] ##exclude P from eigenvalue analysis, don't want it in the Jacobian
    J = ForwardDiff.jacobian(x -> wrapped_model_vec_forced(x, p, t), x)
    λ_max_vals[i] = maximum(real.(eigvals(J)))
end

println("Min λₘₐₓ = ", minimum(λ_max_vals))
println("Max λₘₐₓ = ", maximum(λ_max_vals))
println("Mean λₘₐₓ = ", mean(λ_max_vals))

# Integrate λ_max over the period (e.g., trapezoidal rule)
dt = step(t_eig)
λ_integrated = sum(λ_max_vals[2:end-1]) * dt + 0.5 * dt * (λ_max_vals[1] + λ_max_vals[end])

    ##Calculate jacobian at final state
    x_eval = sol(t_eval)[1:4] ##extract solution at time = 500, when want to evaluate system
    J = ForwardDiff.jacobian(x -> wrapped_model_vec_forced(x, p, t_eval), x_eval)

    #M = cmat(eq, p)
    λ1 = λ1_stability(J)
    λ1_imag = λ1_stability_imag(J)
    react = ν_stability(J)

    return( mean = vec(mean_state),
        sd = vec(sd_state),
        amplitude = vec(range_state),
        min = vec(min_state),
        max = vec(max_state),
        cv = vec(cv),
        λ1 = λ1,
        λ1_imag = λ1_imag,
         λ_integrated = λ_integrated,
        react = react,
        sol = sol)
end 

# ---- CV helper ----
cv(x) = mean(x) == 0 ? NaN : (std(x) / mean(x))

##equilibrium forced without eigenvalue calculation and integration
function equilibrium_forced_2(p, P0; t_warmup = 300.0, t_eval = 500.0, ngrid = 1000, reltol = 1e-8, abstol = 1e-8)
    u0 = @views [1.5, 1.5, 1.0, 1.0, P0] ##initial condition
    tspan = (0.0, t_eval)
    t_grid = range(t_warmup, t_eval, length = ngrid) ##extract dynamics after settling 

    #extract solutions over limit cycle 
    prob = ODEProblem(model_forced!, u0, tspan,p) ##only need to use deepcopy if you are changing parameters inside the function
    sol = solve(prob, Tsit5(); reltol = reltol, abstol = abstol, saveat = t_grid, save_everystep = false, dense = false)
    U = Array(sol)
    
        # Robust CV helper for a vector time series
    coeffvar(v) = begin
        μ = mean(v)
        if !isfinite(μ) || abs(μ) ≤ 1e-10
            NaN
        else
            std(v; corrected=false) / μ
        end
    end

        # State stats over time window
    mean_state = dropdims(mean(U; dims=2), dims=2)
    sd_state   = dropdims(std(U;  dims=2, corrected=false), dims=2)
    min_state  = dropdims(minimum(U; dims=2), dims=2)
    max_state  = dropdims(maximum(U; dims=2), dims=2)
    cv_state   = sd_state ./ mean_state
    cv_state[.!isfinite.(cv_state)] .= NaN   # guard

    # Flux time series into P at each saved time
    # Iterate over columns of U, each is a 5-vector state at time t_grid[j]
    nT = length(t_grid)
    fr_total = Vector{Float64}(undef, nT)
    fr_R1    = similar(fr_total); fr_R2 = similar(fr_total)
    fr_C1    = similar(fr_total); fr_C2 = similar(fr_total)
    fr_G     = similar(fr_total)

    @inbounds for j in 1:nT
        u = @view U[:, j]
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t_grid[j])
        fr_total[j] = tot; fr_R1[j] = r1; fr_R2[j] = r2
        fr_C1[j]    = c1;  fr_C2[j] = c2; fr_G[j]  = g
    end

    return (
        mean      = mean_state,            # 5-vector (R1,R2,C1,C2,P)
        sd        = sd_state,
        amplitude = max_state .- min_state,
        min       = min_state,
        max       = max_state,
        cv        = cv_state,              # state CVs (not a function name!)
        cv_total  = coeffvar(fr_total),    # flux CVs
        cv_R1     = coeffvar(fr_R1),
        cv_R2     = coeffvar(fr_R2),
        cv_C1     = coeffvar(fr_C1),
        cv_C2     = coeffvar(fr_C2),
        cv_G      = coeffvar(fr_G)
    )
end 
#p = ModelPar_active()
#P0 = 0.25
#test = equilibrium_forced_2(p, P0; t_warmup = 300.0, t_eval = 500.0, ngrid = 1000, reltol = 1e-8, abstol = 1e-8)


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



function equilibrium_unforced(p, P0)
    u0 = [1.5, 1.5, 1.0, 1.0, P0] ##initial condition
    t_warmup = 300.0 ##run long enough to reach the equilibrium/limit cycle 
    t_eval = 500.0 ##window to evaluate system properties
    #tspan = (0.0, t_eval)
    #t_grid = range(t_warmup, t_eval, length = 1000) ##extract dynamics after settling 

Δt = (t_eval - t_warmup) / 999           # 1000 samples



prob = ODEProblem(model_unforced!, u0, (0.0, t_eval), deepcopy(p))
sol  = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8,
           saveat=t_warmup:Δt:t_eval, save_everystep=false, dense=false,
            save_idxs=1:5)              # only save the states you need
##look up if function has an option for precision or number of digits used in the computation 
    
function robust_stats(U; clip_negatives=true, tol_abs=1e-4, tol_rel=1e-4, cv_for_absent=NaN)
    U2 = clip_negatives ? max.(U, 0.0) : U            # states×times
    μ  = dropdims(mean(U2; dims=2), dims=2)
    σ  = dropdims(std(U2;  dims=2), dims=2)
    occ = dropdims(maximum(U2; dims=2), dims=2)
    min = dropdims(minimum(U2; dims=2), dims=2)
    cv = similar(μ)
    for i in eachindex(μ)
        tol = max(tol_abs, tol_rel * occ[i])          # relative+absolute floor
        if occ[i] ≤ tol || abs(μ[i]) ≤ tol            # effectively extinct/zero
            cv[i] = cv_for_absent                     # NaN or 0.0
        else
            cv[i] = σ[i] / μ[i]
        end
    end
    return (mean=μ, sd=σ, max=occ, cv=cv, min = min)
end
U = Array(sol)
stats = robust_stats(U; clip_negatives=true, cv_for_absent=NaN)

mean_state = stats.mean
sd_state   = stats.sd
min_state = stats.min
max_state  = stats.max
cv         = stats.cv
   
   #use ODE result as initial guess for equilibrium
    u_approx = sol(t_eval) ##returns full vector of state variables at time t
    Pstar = u_approx[5]
    eq = nlsolve((du, u) -> wrapped_model_unforced!(du, u, deepcopy(p), Pstar), u_approx[1:4]).zero
    cmat(u, p) = ForwardDiff.jacobian(x -> wrapped_model_vec_unforced(x, p, Pstar), u)


##trying to see if can keep function running and see where getting inf/NAs in matrix
    # safe eigenvalue calculation
    λ1 = try
        M = cmat(eq, p)
        if any(!isfinite, M)
            @warn "NaN/Inf in community matrix" p=p eq=eq
            NaN
        else
            maximum(real.(eigvals(M)))
        end
    catch e
        @warn "λ1 calculation failed" p=p eq=eq exception=(e, catch_backtrace())
        NaN
    end

    λ1_imag = try
        M = cmat(eq, p)
        if any(!isfinite, M)
            NaN
        else
            vals = eigvals(M)
            imag(vals[argmax(real.(vals))])
        end
    catch
        NaN
    end

    react = try
        M = cmat(eq, p)
        if any(!isfinite, M)
            NaN
        else
            ν_stability(M)
        end
    catch
        NaN
    end

    return(eq = eq, λ1 = λ1, λ1_imag =λ1_imag, react = react, cv=cv, min = min_state, max = max_state,  mean = mean_state, sd = sd_state)
end 

###Functions to calculate CV of total harvest
# ---- CV helper ----
cv(x) = mean(x) == 0 ? NaN : (std(x) / mean(x))

# ---- run, sample post-warmup, compute FR_into_P series + CVs ----
function fr_cv_unforced(p; u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    tspan = (0.0, t_eval)
    prob  = ODEProblem(model_unforced!, u0, tspan, deepcopy(p))
    sol   = solve(prob; reltol=1e-8, abstol=1e-8)

    t_grid  = range(t_warmup, t_eval; length=ngrid)
    us      = sol.(t_grid)

    # time series of fluxes into P
    fr_total = Float64[]; fr_R1 = Float64[]; fr_R2 = Float64[]
    fr_C1 = Float64[];    fr_C2 = Float64[]; fr_G  = Float64[]

    for (u, t) in zip(us, t_grid)
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t)
        push!(fr_total, tot); push!(fr_R1, r1); push!(fr_R2, r2)
        push!(fr_C1, c1);     push!(fr_C2, c2); push!(fr_G,  g)
    end

    return (; 
        cv_total = cv(fr_total),
        cv_R1    = cv(fr_R1),
        cv_R2    = cv(fr_R2),
        cv_C1    = cv(fr_C1),
        cv_C2    = cv(fr_C2),
        cv_G     = cv(fr_G)
    )
end

function fr_cv_forced(p; u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    tspan = (0.0, t_eval)
    prob  = ODEProblem(model_forced!, u0, tspan, deepcopy(p))
    sol   = solve(prob; reltol=1e-8, abstol=1e-8)

    t_grid  = range(t_warmup, t_eval; length=ngrid)
    us      = sol.(t_grid)

    # time series of fluxes into P
    fr_total = Float64[]; fr_R1 = Float64[]; fr_R2 = Float64[]
    fr_C1 = Float64[];    fr_C2 = Float64[]; fr_G  = Float64[]

    for (u, t) in zip(us, t_grid)
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t)
        push!(fr_total, tot); push!(fr_R1, r1); push!(fr_R2, r2)
        push!(fr_C1, c1);     push!(fr_C2, c2); push!(fr_G,  g)
    end

    return (; 
        cv_total = cv(fr_total),
        cv_R1    = cv(fr_R1),
        cv_R2    = cv(fr_R2),
        cv_C1    = cv(fr_C1),
        cv_C2    = cv(fr_C2),
        cv_G     = cv(fr_G),
        mean_total = mean(fr_total),
        mean_R1 = mean(fr_R1),
        mean_R2 = mean(fr_R2),
        mean_C1 = mean(fr_C1),
        mean_C2 = mean(fr_C2),
        mean_G = mean(fr_G),
        sd_total = std(fr_total),
        sd_R1 = std(fr_R1),
        sd_R2 = std(fr_R2),
        sd_C1 = std(fr_C1),
        sd_C2 = std(fr_C2),
        sd_G = std(fr_G),
    )
end

##STRUCTURE 1: Plotting dynamics, equilibrium, eigenvalue analysis 
##Solve ODE 
##set initial condition
u0 = [1.5, 1.5, 1.0, 1.0, 0.25]
#u0 = [0.6, 0.8, 0.45, 0.61, 0.2]

tspan = (0.0, 500.0)
G_pre = 1.0 
G_pulse = G_pre
t_pulse = 200.0 ##time when disturbance occurs, want to be once model at equilibirum
t_recover = 250.0 ##time when decline in resources ends 
 
##set Parameters
#p = ModelPar_active(o = 0.1, w = 0.5, H = 0.0, K = 3.05)
p = ModelPar_active(o = 0.0, w = 0.0, H = 0.0, r = 1.2999999999999996, K = 1.9, aR_C = 1.8013606526895989, aR_P = 1.1777706183189485, aC_P = 1.7486164433891764, aG_P = 1.6654077546996051, hR_P = 1.7239456424939552, hR_C = 0.47568284600108846, hC_P = 1.0250642119658746, hG_P = 1.395392617078673, e = 0.9249999999999999, mC = 0.6462558154613841, G = 3.714285714285714)
##Define the ODE problem
prob_1 = ODEProblem(rhs_unforced, u0, tspan, p)
sol_1 = solve(prob_1)

##also this ODE solver is working, why in function am i then getting NAs/Infs in matrix? 

##plot timeseries
plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution - Active Omnivory")

plot(sol_1, tspan=(400, 500),
     xlabel="Time", ylabel="Population",
     title="ODE Solution - Active Omnivory")

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
   colors = [:black, :darkgreen, :salmon, :lightgreen, :pink, :blue]

plot(times, fr_total, label = false, xlabel = "Time", ylabel = "Total Harvest", ylims = (0, 1.0), color = colors[1], linewidth = 2.5)  
plot!(times, fr_R1, label = "R1 → P", color = colors[2], linewidth = 2.5)
plot!(times, fr_R2, label = "R2 → P", color = colors[3], linewidth = 2.5)
plot!(times, fr_C1, label = "C1 → P", color = colors[4], linewidth = 2.5)
plot!(times, fr_C2, label = "C2 → P", color = colors[5], linewidth = 2.5)
plot!(times, fr_G, label = "G → P", color = colors[6], linewidth = 2.5)

#plot without legend
plot(times, fr_total, xlabel = "Time", ylabel = "Predator Consumption", ylims = (0, 1.0), color = colors[1], linewidth = 2.5, legend = false)  
plot!(times, fr_R1,  color = colors[2], linewidth = 2.5)
plot!(times, fr_R2,  color = colors[3], linewidth = 2.5)
plot!(times, fr_C1,  color = colors[4], linewidth = 2.5)
plot!(times, fr_C2, color = colors[5], linewidth = 2.5)
plot!(times, fr_G, color = colors[6], linewidth = 2.5)
plot!(times, fr_total,color = colors[1], linewidth = 2.5)  



##Plot consumption dynamics after reaching equilibrium
# Find index where time >= 250
start_idx = findfirst(t -> t ≥ 250, times)

# Subset time and consumption arrays
times_sub_4 = times[start_idx:end]
fr_total_sub = fr_total[start_idx:end]
fr_R1_sub = fr_R1[start_idx:end]
fr_R2_sub = fr_R2[start_idx:end]
fr_C1_sub = fr_C1[start_idx:end]
fr_C2_sub = fr_C2[start_idx:end]
fr_G_sub = fr_G[start_idx:end]

# Define custom colors (optional - adjust as needed)
colors = [:black, :darkgreen, :salmon, :lightgreen, :pink, :red]

# Plot
plot(times_sub_4, fr_total_sub, label = "total → P", xlabel = "Time", ylabel = "Community Consumption", ylims = (0.0, 1.25), color = colors[1], linewidth = 2.5)
plot!(times_sub_4, fr_R1_sub, label = "R1 → P", color = colors[2],linewidth = 2.5)
plot!(times_sub_4, fr_R2_sub, label = "R2 → P", color = colors[3],linewidth = 2.5)
plot!(times_sub_4, fr_C1_sub, label = "C1 → P", color = colors[4],linewidth = 2.5)
plot!(times_sub_4, fr_C2_sub, label = "C2 → P", color = colors[5],linewidth = 2.5)
plot!(times_sub_4, fr_G_sub, label = "G → P", color = colors[6],linewidth = 2.5)

















##Look at dynamics with active omnivory 
##set Parameters
#P_fixed = 0.25
#u0 = [1.5, 1.5, 1.0, 1.0, P_fixed]
tspan = (0.0, 100.0)
p = ModelPar_active(w = 0.1, o = 0.1, H = 0.9, K =7.5, pf = 10.0, D = 0.5)

##Define the ODE problem
prob_2 = ODEProblem(rhs_forced, u0, tspan, p)
sol_2 = solve(prob_2)

##plot timeseries
plot(sol_2, xlabel="Time", ylabel="Population", title="ODE Solution - Active Omnivory")
  
##unforced
prob_3 = ODEProblem(rhs_unforced, u0, tspan, p)
sol_3 = solve(prob_3)

##plot timeseries
plot(sol_3, xlabel="Time", ylabel="Population", title="ODE Solution - Active Omnivory")
  
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


 