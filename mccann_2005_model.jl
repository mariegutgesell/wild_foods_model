##Reproducing model from McCann et al., 2005
##want to ensure that I can reproduce results as expected, before adding in subsidy and holding P constant

##Date Initiated: June 18, 2025
##Contributor(s): Marie K. Gutgesell

using Pkg
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

##Habitat preference (coupling)
##habitat preference - so preference for habitat changes with shifting densities of C across patches (but preference w itself is fixed, i.e., when prey are at equal densities, this is pref for patch1, influences speed of switching across patches)
function hab_pref(u, p, t)
    R1, R2, C1, C2 = u
    return p.w * C1/(p.w * C1 + (1 - p.w)*C2)  
end

##Scale of foraging 
function forage_scale(u, p, t)
    W = p.h_pref(u, p, t)
    return p.Q + (1 - p.Q) * W
end 

##So hab_pref = Wi, and so 1-hab_pref is Wj, and since Si = Wi when Q = 0, then Sj = 1 - Si (or 1 - Wj)
#When Q = 1, functional response returns to classic multi-species functional response where consumer percieves prey as well-mixed
#When Q = 0, consumer must make a choice between foraging in different habitats 
##We will start with case where Sc = Sh (i.e., Q= 0), as this most reflects humans making a choice between foraging in aquatic vs. terrestrial habitats

# Degree of omnivory functions -- where o represents preference, based on McCann et al., 2005 
function deg_om_i(u, p, t) ##constant speed of switching, where degree of omnivory shifts in response to changing resource/consumer densities in given patch i -  and degree of omnivory used in functional response (can consider later adding in dynamic omnivory)
   R1, C1 = u
    return (p.o * p.e * p.aP * R1) / (p.o * p.e * p.aP * R1 + (1 - p.o)*p.e * p.aP * C1)
end


function deg_om_j(u, p ,t)
    R2, C2 = u
    return (p.o * p.e * p.aP * R2) / (p.o * p.e * p.aP * R2 + (1 - p.o)*p.e * p.aP * C2)
end




##functional response between resources and predator
function f_R1P(u, p, t)
    @unpack aR_P, hR_P, aC_P, hC_P= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2 = u  ##defines state variables, G = groceries
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 

    numerator = S1 * aR_P * o1 * R1
    denominator = 1 + (S1 * aR_P * hR_P * o1 * R1 + (1-S1)* aR_P * hR_P * o2 * R2 + S1 * aC_P * hC_P * (1-o1) * C1 + (1-S1) * aC_P * hC_P * (1-o2) * C2)
    return numerator / denominator
end

function f_R2P(u, p, t)
    @unpack aR_P, hR_P, aC_P, hC_P= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 

    numerator = (1 - S1) * aR_P * o2 * R2
   denominator = 1 + (S1 * aR_P * hR_P * o1 * R1 + (1-S1)* aR_P * hR_P * o2 * R2 + S1 * aC_P * hC_P * (1-o1) * C1 + (1-S1) * aC_P * hC_P * (1-o2) * C2)
     return numerator / denominator
end

function f_C1P(u, p, t)
    @unpack aR_P, hR_P, aC_P, hC_P= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 

    numerator = S1 * aC_P * (1 - o1) * C1
   denominator = 1 + (S1 * aR_P * hR_P * o1 * R1 + (1-S1)* aR_P * hR_P * o2 * R2 + S1 * aC_P * hC_P * (1-o1) * C1 + (1-S1) * aC_P * hC_P * (1-o2) * C2)
     return numerator / denominator
end

function f_C2P(u, p, t)
    @unpack aR_P, hR_P, aC_P, hC_P= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 

    numerator = (1 - S1) * aC_P * (1 - o2) * C2
   denominator = 1 + (S1 * aR_P * hR_P * o1 * R1 + (1-S1)* aR_P * hR_P * o2 * R2 + S1 * aC_P * hC_P * (1-o1) * C1 + (1-S1) * aC_P * hC_P * (1-o2) * C2)
     return numerator / denominator
end


#Set up paramters 
@with_kw mutable struct ModelPar
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 0.7   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
            ##when o = 0, no omnivory (P only feeds on C)
    ##Model parameters, for now just keeping these parameters the same (could make unique ones based on patch/TL)
    r = 1.0
    K = 2.25
    aR_C = 2.5  ##attack rate  of consumer on R
    aR_P = 3.4 ##attack rate of P on R 
    aC_P =   ##attack rate of P on C 
    e = 0.8   ##energy conversion 
    mC = 1.0 ##C mortality rate
    mP = 0.45   ## P mortality rate
    hR_C = 0.4 ##handling time of C on R
    hR_P = 1.25  ##handling time of P on R  
    hC_P =  ##handling time of P on C
    ##Habitat preference function 
    h_pref::Function = hab_pref  

    ##Foraging scale
    S::Function = forage_scale
    
    ##Omnivory preference function
    d_om_i::Function = deg_om_i ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = deg_om_j ##degree of omnivory in patch 2 


    ##Predator functional responses
    f_r1p::Function = f_R1P
    f_r2p::Function = f_R2P
    f_c1p::Function = f_C1P
    f_c2p::Function = f_C2P

end

##McCann et al., 2005 model w/ omnivory
function model_2!(du, u, p ,t)
    @unpack r, K, aC, aP,  e, hC, hP, mC, mP = p
   R1, R2, C1, C2, P = u 
   #S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   #o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   #o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
  
## consumer functional responses
    f_RC1 = aC * R1 / (1 + aC * hC * R1)
   f_RC2 = aC * R2 / (1 + aC * hC * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * f_r2p
   du[3] = e * C1 * f_RC1 - P * f_c1p - mC * C1
   du[4] = e * C2 * f_RC2 - P * f_c2p - mC * C2 
   du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P  - mP * P
   # du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

   return du
 end 


# Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_2(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_2!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

##STRUCTURE 1: Plotting dynamics, equilibrium, eigenvalue analysis 
##Solve ODE 
##set initial condition
u0 = [0.6, 0.8, 0.45, 0.61, 0.2]
tspan = (0.0, 1000.0)
##set Parameters
p = ModelPar()

##Define the ODE problem
prob_1 = ODEProblem(rhs_2, u0, tspan, p)
sol_1 = solve(prob_1)

##plot timeseries
   plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution")
  