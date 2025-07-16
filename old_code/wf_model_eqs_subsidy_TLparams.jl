##Code to set up equations for wild food model - WITH SUBSIDY
##Date Initiated: April 8, 2025
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

#To look into:
##McCann 2005 or Gutgesell 2022 omnivory approach 

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


### Omnivory preference functions -- density dependent preference, where o measures the speed of switching, based on Gutgesell et al., 2022
function om_i_pref_active(u, p, t)
    R1, C1 = u
    return (p.o * R1) / (p.o * R1 + (1-p.o) * C1) 
end

function om_j_pref_active(u, p, t)
    R2, C2 = u
    return (p.o * R2) / (p.o * R2 + (1-p.o) * C2) 
end

function om_i_pref_fixed(u,p,t)
    return p.o_fixed
end

function om_j_pref_fixed(u,p,t)
    return p.o_fixed
end
##for now, using degree of omnivory as in McCann et al., 2005, think about using omnivory preference later ... come back to 

##Grocery Preference -- also density dependent, and H measures speed of switching 
function sub_pref_func(u, p ,t)
    R1, R2, C1, C2 = u 
    return(p.H * p.G) / (p.H * p.G + (1-p.H)*R1 + (1-p.H)*R2 + (1-p.H) * C1 + (1-p.H)*C2)
end

##functional responses between resources and consumer 
function f_R1C1(u, p, t)
    @unpack aRC, hRC, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2 = u  ##defines state variables, G = groceries
    
    return aRC * R1 / (1 + aRC * hRC * R1)
end

function f_R2C2(u, p, t)
    @unpack aRC, hRC, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2 = u  ##defines state variables, G = groceries
    
    return aRC * R2 / (1 + aRC * hRC * R2)
end

##functional response between resources and predator
function f_R1P(u, p, t)
    @unpack  aRP, aCP, aGP, hRP, hCP, hGP, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2 = u  ##defines state variables, G = groceries
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    
    numerator = (1-H1) * S1 * aRP * o1 * R1
    denominator = 1 + (S1 * aRP * hRP * o1 * R1 + (1-S1)* aRP * hRP * o2 * R2 + S1 * aCP * hCP * (1-o1) * C1 + (1-S1) * aCP * hCP * (1-o2) * C2 + H1 * aGP * hGP * G)
    return numerator / denominator
end

function f_R2P(u, p, t)
    @unpack  aRP, aCP, aGP, hRP, hCP, hGP, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * (1 - S1) * aRP * o2 * R2
    denominator = 1 + (S1 * aRP * hRP * o1 * R1 + (1-S1)* aRP * hRP * o2 * R2 + S1 * aCP * hCP * (1-o1) * C1 + (1-S1) * aCP * hCP * (1-o2) * C2 + H1 * aGP * hGP * G)
     return numerator / denominator
end

function f_C1P(u, p, t)
    @unpack  aRP, aCP, aGP, hRP, hCP, hGP, G = p##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * S1 * aCP * (1 - o1) * C1
    denominator = 1 + (S1 * aRP * hRP * o1 * R1 + (1-S1)* aRP * hRP * o2 * R2 + S1 * aCP * hCP * (1-o1) * C1 + (1-S1) * aCP * hCP * (1-o2) * C2 + H1 * aGP * hGP * G)
    return numerator / denominator
end

function f_C2P(u, p, t)
    @unpack  aRP, aCP, aGP, hRP, hCP, hGP, G = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = (1-H1) * (1 - S1) * aCP * (1 - o2) * C2
    denominator = 1 + (S1 * aRP * hRP * o1 * R1 + (1-S1)* aRP * hRP * o2 * R2 + S1 * aCP * hCP * (1-o1) * C1 + (1-S1) * aCP * hCP * (1-o2) * C2 + H1 * aGP * hGP * G)
    return numerator / denominator
end

##trying out functional resposne for G.. type 2 functional response 
function f_GP(u, p, t)
    @unpack  aRP, aCP, aGP, hRP, hCP, hGP, G = p  ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)

    numerator = H1 * aGP * G #trying if i remove scaling / suppression of G by other foraging preferences, i think this makes biological sense (but keep 1-H1 in other FRs)
    denominator = 1 + (S1 * aRP * hRP * o1 * R1 + (1-S1)* aRP * hRP * o2 * R2 + S1 * aCP * hCP * (1-o1) * C1 + (1-S1) * aCP * hCP * (1-o2) * C2 + H1 * aGP * hGP * G)
    return numerator / denominator
end


##linear (type 1) functional response 
function f_GP_2(u, p, t)
    @unpack a, h, H, G = p 
    R1, R2, C1, C2 = u 
    return G * H
end


##second form of parameters, where omnivory is done using Gutgesell et al., active omnivory, and different parameters for each trophic level, but same in each patch  
@with_kw mutable struct ModelPar2 
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 1.0   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.0 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)
    o_fixed = 0.1
    H = 0.5 ##H = preference for groceries (G)

    ##Model parameters, for now just keeping these parameters the same in patch, unique across TL (values based on Gutgesell et al., 2022)
    r = 2.0
    K = 3.0
    aRC = 1.0  ##attack rate of C on R 
    aRP = 0.2 ##attack rate of P on R 
    aCP = 0.5 ##attack rate of P on C 
    aGP = 0.5 ##attack rate of P on G (groceries) 
    e = 0.8   ##energy conversion - same for all (as in McCann et al., 2005)
    mC = 0.4
    mP = 0.2
    hRC = 0.4 ##handling time of C consuming R 
    hRP = 0.6 ##handling time of P consuming R
    hCP = 0.3 ##handling time of P consuming C 
    hGP = 0.3 ##handling time of P consuming G (groceries) - just keeping same as C 
    ##Habitat preference function 
    h_pref::Function = hab_pref  

    ##Foraging scale
    S::Function = forage_scale
    
    ##Omnivory preference function
    d_om_i::Function = om_i_pref_active ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = om_j_pref_active ##degree of omnivory in patch 2 
    
    ##Subsidy preference function
    sub_pref::Function = sub_pref_func

    ##Initial/constant value of groceries, so effectively instantly replenishes 
    G = 2.0
    ##Predator functional responses
    f_r1p::Function = f_R1P
    f_r2p::Function = f_R2P
    f_c1p::Function = f_C1P
    f_c2p::Function = f_C2P
    f_gp::Function = f_GP
    f_gp2::Function = f_GP_2
end

##Model that holds P constant 
function model_2!(du, u, p ,t)
    @unpack r, K, aRC, aRP, aCP, aGP, hRC, hRP, hCP, hGP, e, mC, mP, H, G = p
   R1, R2, C1, C2, P = u 
   #S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   #o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   #o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   
   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
     f_gp =p.f_gp(u, p, t)
   # f_gp2 = p.f_gp2(u, p, t)
## consumer functional responses
    f_RC1 = aRC * R1 / (1 + aRC * hRC * R1)
   f_RC2 = aRC * R2 / (1 + aRC * hRC * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * f_r2p
   du[3] = e * C1 * f_RC1 - P * f_c1p - mC * C1
   du[4] = e * C2 * f_RC2 - P * f_c2p - mC * C2 
  # du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp2 - m * P
    du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 


# Utilities for doing eigenvalue analysis
   ##this is a right-hand-side function that is used to solve ODEs 
   function rhs_2(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model_2!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end


##Note: next if using NLSolve need to do wrapper function (not working yet)
function wrapped_model!(residuals, x, p, P_fixed)
    R1, R2, C1, C2 = x
    T = eltype(x)
    P = convert(T, P_fixed)
    u = T[R1, R2, C1, C2, P]  # use T[...] to create uniform type array

    du = zeros(T,5)
    model_2!(du, u, p, 0.0)

    residuals[1:4] .= du[1:4]
    return residuals
end
##vector-returning version, to allow for jacobian calculation
function wrapped_model_vec(x, p, P_fixed)
    T = eltype(x)
    du = zeros(T, 4)
    wrapped_model!(du, x, p, P_fixed)
    return du
end


##function to calculate total harvest for P, by summing functional response for P at each time step
function total_FR_into_P(u, p ,t)
 @unpack aRC, aRP, aCP, aGP, hRC, hRP, hCP, hGP, e, mC, mP, H, G = p
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


##STRUCTURE 1: Plotting dynamics, equilibrium, eigenvalue analysis 
##Solve ODE 
##set initial condition
u0 = [2.0, 2.0, 1.0, 1.0, 0.5]
tspan = (0.0, 1000.0)
##set Parameters
p = ModelPar2()

##Define the ODE problem
prob_1 = ODEProblem(rhs_2, u0, tspan, p)
sol_1 = solve(prob_1)

##plot timeseries
   plot(sol_1, xlabel="Time", ylabel="Population", title="ODE Solution")
  
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


##Calculating local equilibrium and stability (eigenvalues)
# Initial guess for dynamic variables: R1, R2, C1, C2
x0 = [1.0, 1.0, 0.5, 0.5]
P_fixed = 0.5
p = ModelPar2()

##Local stability analysis
#1) Calculate equilibrium values for system of equation 
result = nlsolve(x -> wrapped_model_vec(x, p, P_fixed), x0)

# Extract equilibrium values
x_star = result.zero
R1_star, R2_star, C1_star, C2_star = x_star
 

#2) Calculate numerical jacobian
J_num = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, P_fixed), x0)
println("Jacobian Matrix:")
display(J_num)

#3) Calculate the eigenvalues 
eigvals(J_num)


##Trying to find multiple equilibriums, as likely more than 1 
#Pkg.add("IterTools")
using IterTools

R_vals = [0.1, 0.5, 1.0, 2.0]
C_vals = [0.1, 0.5, 1.0]

# All combinations of R1, R2, C1, C2
guesses = collect(product(R_vals, R_vals, C_vals, C_vals))

P_fixed = 0.5
tolerance = 1e-4
equilibria_1 = []

for x0 in guesses
    result = nlsolve(x -> wrapped_model_vec(x, p, P_fixed), collect(x0))
    if result.f_converged
        x_star = result.zero
        # Check if it's a new equilibrium
        if all(norm(x_star .- eq) > tolerance for eq in equilibria_1)
            push!(equilibria_1, x_star)
        end
    end
end


max_eigenvalues_1 = Float64[]  # to store max real parts

println("Number of distinct equilibria found: ", length(equilibria_1))

for (i, eq) in enumerate(equilibria_1)
    J = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, P_fixed), eq)
    λ = eigvals(J)
    λ_real = real.(λ)
    λ_max = maximum(λ_real)

    is_stable = all(λ_real .< 0)
    stability = is_stable ? "Stable" : "Unstable or saddle"
    
    println("\nEquilibrium $i:")
    println("  State: ", round.(eq; digits=4))
    println("  Stability: ", stability)
    println("  Max real eigenvalue: ", round(λ_max; digits=6))
    println("  All eigenvalues: ", round.(λ; digits=4))

    push!(max_eigenvalues_1, λ_max)
end

##For any equilibirum where the state value is negative, is biologically impossible (can't have densities <0)
##Focus only on stable equilibrium? 

##Later: look at isoclines of pairs to dig into potenial eqiulibrium points 

##
