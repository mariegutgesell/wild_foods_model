##Code to unpack and look at functional responses in model, with preference for patch and omnivory 
##Date Initiated: April 1, 2025
##Contributor(s): Marie K. Gutgesell

#using Pkg
#Pkg.add("Plots")
##Load Libraries needed (ODE for later)
using Parameters: @with_kw, @unpack ##imports Parameters package that provides convenient macros for working with keyword arugemnts, parameter structs and unpacking variables 
using DifferentialEquations
using Plots
#using Symbolics
using ForwardDiff
using LinearAlgebra
using NLsolve
#using DataFrames
#using Interact

##Re-creating model from McCann et al., 2005, Ecology Letters
#u = state variables where u[1] = R1, u[2] = R2, u[3] = C1, u[4] = C2, u[5] = P
#p = Parameters
#t = 

##Habitat preference (coupling)
##habitat preference - so preference for habitat changes with shifting densities of C across patches (but preference w itself is fixed)
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

# Omnivory preference functions
function deg_om_i(u, p, t) ##constant preference, where degree of omnivory shifts in response to changing resource/consumer densities in given patch i -  and degree of omnivory used in functional response (can consider later adding in dynamic omnivory)
   R1, C1 = u
    return (p.o * p.e * p.a * R1) / (p.o * p.e * p.a * R1 + (1 - p.o)*p.e * p.a * C1)
end


function deg_om_j(u, p ,t)
    R2, C2 = u
    return (p.o * p.e * p.a * R2) / (p.o * p.e * p.a * R2 + (1 - p.o)*p.e * p.a * C2)
end





##functional response between resources and predator
function f_R1P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    numerator = S1 * a * o1 * R1
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
    return numerator / denominator
end

function f_R2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    numerator = (1 - S1) * a * o2 * R2
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
    return numerator / denominator
end

function f_C1P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    numerator = S1 * a * (1 - o1) * C1
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
    return numerator / denominator
end

function f_C2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    numerator = (1 - S1) * a * (1 - o2) * C2
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
    return numerator / denominator
end

#Set up paramters 
@with_kw mutable struct ModelPar
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 0.8   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.2  ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)

    ##Model parameters, for now just keeping these parameters the same (could make unique ones based on patch/TL)
    r = 1.0
    K = 3.25
    a = 2.5  ##attack rate 
    e = 0.8   ##energy conversion 
    m = 1.0
    h = 0.5
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

##Plotting functional responses across range of R/C
p = ModelPar()

##Plot P functional response across range of R1, holding all other variables constant
R1_vals = range(0, 50, length = 100)
response_vals = [f_R1P((R1, 10.0, 5.0, 5.0), p, 0.0) for R1 in R1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R1_vals, response_vals,
xlabel = "R1 Density", ylabel = "Predator Consumption Rate of R1", title = "Functional Response to R1")



##Plot P functional response across range of R2, holding all other variables constant
R2_vals = range(0, 100, length = 100)
response_vals_R2 = [f_R2P((10.0, R2, 5.0, 5.0), p, 0.0) for R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R2_vals, response_vals_R2,
xlabel = "R2 Density", ylabel = "Predator Consumption Rate of R2", title = "Functional Response to R2")


##Plot P functional response across range of R1, holding all other variables constant
C1_vals = range(0, 50, length = 100)
response_vals_C1 = [f_C1P((10.0, 10.0, C1, 5.0), p, 0.0) for C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C1_vals, response_vals_C1,
xlabel = "C1 Density", ylabel = "Predator Consumption Rate of C1", title = "Functional Response to C1")



##Plot P functional response across range of R1, holding all other variables constant
C2_vals = range(0, 100, length = 100)
response_vals_C2 = [f_C2P((10.0, 10.0, 5.0, C2), p, 0.0) for C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C2_vals, response_vals_C2,
xlabel = "C2 Density", ylabel = "Predator Consumption Rate of C2", title = "Functional Response to C2")

##now getting basically type 2 / type 3 functional response... is that right? yes i think so.. 


##Looking into S1, o1 and o2, to see if they are responding how I think they should be across C/R densities
C1_vals = range(0.1, 100, length = 100)
S1_vals_C1 = [p.S((10.0, 10.0, C1, 5.0), p, 0.0) for C1 in C1_vals]
plot(C1_vals, S1_vals_C1, xlabel = "C1 Density", ylabel = "S1 (foraging in Patch 1)", title = "Foraging Preference (S1) vs. C1 Density")

C2_vals = range(0.1, 100, length = 100)
S1_vals_C2 = [p.S((10.0, 10.0, 5.0, C2), p, 0.0) for C2 in C2_vals]
plot(C2_vals, S1_vals_C2, xlabel = "C2 Density", ylabel = "S1 (foraging in Patch 1)", title = "Foraging Preference (S1) vs. C2 Density")


R1_vals = range(0.1, 100, length = 100)
S1_vals_R1 = [p.S((R1, 10.0, 5.0, 5.0), p, 0.0) for R1 in R1_vals]
plot(R1_vals, S1_vals_R1, xlabel = "R1 Density", ylabel = "S1 (foraging in Patch 1)", title = "Foraging Preference (S1) vs. R1 Density")
##so habitat preference is driven by differences in C densities, so P needs to make a choice, and R density has no effect on foraging 


##Looking at o1
o1_vals_R1 = [p.d_om_i((R1, 10.0), p, 0.0) for R1 in R1_vals]
plot(R1_vals, o1_vals_R1, xlabel = "R1 Density", ylabel = "Degree of Omnivory in Patch i")
##increasing o (higher preference for R) makes degree of omnivory increase rapidly, so working as expected

o1_vals_C1 = [p.d_om_i((10.0, C1), p, 0.0) for C1 in C1_vals]
plot(C1_vals, o1_vals_C1, xlabel = "C1 Density", ylabel = "Degree of Omnivory in Patch i")



o1_vals_R2 = [p.d_om_i((10.0, 5.0), p, 0.0) for R2 in R2_vals]
plot(R2_vals, o1_vals_R2, xlabel = "R2 Density", ylabel = "Degree of Omnivory in Patch i")



##testing model by calling above functional responses ... 
##trying model that holds dP/dt constant ... and maybe calculates changes in harvest? sum of functional responses.. 
function model_2!(du, u, p ,t)
    @unpack r, K, a, e, h, m = p
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
    f_RC1 = a * R1 / (1 + a * h * R1)
   f_RC2 = a * R2 / (1 + a * h * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * f_r2p
   du[3] = e * C1 * f_RC1 - P * f_c1p - m * C1
   du[4] = e * C2 * f_RC2 - P * f_c2p - m * C2 
   du[5] = 0

   return du
 end 


  # Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_2(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_2!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

##set initial condition
u0 = [0.6, 0.8, 0.45, 0.61, 0.2]
tspan = (0.0, 1000.0)

##set Parameters
p = ModelPar()


##Define the ODE problem
prob_2 = ODEProblem(rhs_2, u0, tspan, p)
sol_2 = solve(prob_2)


   ##plot timeseries
   plot(sol_2, xlabel="Time", ylabel="Population", title="ODE Solution")
  

    ##function to calculate total harvest for P, by summing functional response for P at each time step
function total_FR_into_P(u, p ,t)
    @unpack a, h = p
   R1, R2, C1, C2, P = u 
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)

    total_FR = f_r1p + f_r2p + f_c1p + f_c2p
   return (total_FR, f_r1p, f_r2p, f_c1p, f_c2p)
end 

fr_vals_2 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)]

# Separate each series into its own array
fr_total = [x[1] for x in fr_vals_2]
fr_R1 = [x[2] for x in fr_vals_2]
fr_R2 = [x[3] for x in fr_vals_2]
fr_C1 = [x[4] for x in fr_vals_2]
fr_C2 = [x[5] for x in fr_vals_2]
times = sol_2.t

   ##plot total P consumption 
   plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Predator Consumption")  
   plot!(times, fr_R1, label = "R1 → P")
plot!(times, fr_R2, label = "R2 → P")
plot!(times, fr_C1, label = "C1 → P")
plot!(times, fr_C2, label = "C2 → P")


##lets test pushing this to git 