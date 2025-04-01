##Code to set up and explore equations for wild food model
##Date Initiated: October 3, 2024
##Contributor(s): Marie K. Gutgesell

using Pkg
Pkg.add("Plots")
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

#Set up paramters 
@with_kw mutable struct ModelPar
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 0.95   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.95  ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)

    ##Model parameters, for now just keeping these parameters the same (could make unique ones based on patch/TL)
    r = 1.0
    K = 3.25
    a = 0.3  ##attack rate 
    e = 0.8   ##energy conversion 
    m = 1.0
    h = 1.0
    ##Habitat preference function 
    h_pref::Function = hab_pref  

    ##Foraging scale
    S::Function = forage_scale
    
    ##Omnivory preference function
    d_om_i::Function = deg_om_i ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function = deg_om_j ##degree of omnivory in patch 2 
    
end

 

##trying model that holds dP/dt constant ... and maybe calculates changes in harvest? sum of functional responses.. 
function model_2!(du, u, p ,t)
    @unpack r, K, a, e, h, m = p
   R1, R2, C1, C2, P = u 
   S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   
   ##compute these once for all... essentially same as in functional responses above 
   num_R1P = S1 * a * o1 * R1
   num_R2P = (1 - S1) * a * o2 * R2
   num_C1P = S1 * a * (1 - o1) * C1
   num_C2P = (1 - S1) * a * (1 - o2) * C2
   denom_RCP =  1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
   f_RC1 = a * R1 / (1 + a * h * R1)
   f_RC2 = a * R2 / (1 + a * h * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * num_R1P / denom_RCP
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * num_R2P / denom_RCP
   du[3] = e * C1 * f_RC1 - P * num_C1P / denom_RCP - m * C1
   du[4] = e * C2 * f_RC2 - P * num_C2P / denom_RCP - m * C2 
   du[5] = 0

   return du
 end 
   
 ##function to calculate total harvest for P, by summing functional response for P at each time step
function total_FR_into_P(u, p ,t)
    @unpack a, h = p
   R1, R2, C1, C2, P = u 
   S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   
   ##compute these once for all... essentially same as in functional responses above 
   num_R1P = S1 * a * o1 * R1
   num_R2P = (1 - S1) * a * o2 * R2
   num_C1P = S1 * a * (1 - o1) * C1
   num_C2P = (1 - S1) * a * (1 - o2) * C2
   denom_RCP =  1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)

   ##Individual functional responses 
    fr_R1P = P*num_R1P/denom_RCP
    fr_R2P = P*num_R2P/denom_RCP
    fr_C1P = P*num_C1P/denom_RCP
    fr_C2P = P*num_C2P/denom_RCP

    total_FR = fr_R1P + fr_R2P + fr_C1P + fr_C2P
   return (total_FR, fr_R1P, fr_R2P, fr_C1P, fr_C2P)
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
   
   fr_vals_2 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)]

   # Separate each series into its own array
   fr_total = [x[1] for x in fr_vals_2]
fr_R1 = [x[2] for x in fr_vals_2]
fr_R2 = [x[3] for x in fr_vals_2]
fr_C1 = [x[4] for x in fr_vals_2]
fr_C2 = [x[5] for x in fr_vals_2]
times = sol_2.t

   ##plot timeseries
   plot(sol_2, xlabel="Time", ylabel="Population", title="ODE Solution")
   
   ##plot total P consumption 
   plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Functional Response")  
   plot!(times, fr_R1, label = "R1 → P")
plot!(times, fr_R2, label = "R2 → P")
plot!(times, fr_C1, label = "C1 → P")
plot!(times, fr_C2, label = "C2 → P")


##Look at functional responses here... as weird that populations are just crashing 
C1_vals = range(0.01, 5.0, length = 100)
fr_vals_c1 = [P * (S1 * a * (1 - o1) * C1)/1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2) for C1 in C1_vals]




  ##seing if i can track parameters in the model over time to see how the model is working
  S_vals = [p.S(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)]   
  o1_vals = [p.d_om_i(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)] 
  o2_vals = [p.d_om_(u, p, t) for (u, t) in zip(sol_2.u, sol_2.t)] 
##Next steps:
##work through checks to make sure the math makes sense
##write out and think about the math to get to the jacobian and the eigenvalues to interpret stability outcomes
##sensitivity analysis across parameter spaces ... 
##bifurcations
##stability metrics -- local and non-local 

##how to add in preference for subsidy
##calculating total harvest -- essentially as sum of P*functional responses ... 

 ##thinking about parameters, and then how to set up experiments to test applied questions 
 ##introducing abiotic asynchrony into resources 
 ##making model stochastic 


 ##much later next steps: thinking about perturbation experiments.. 





 ###CODE WITH P NOT HELD constant


 function model!(du, u, p ,t)
    @unpack r, K, a, e, h, m = p
   R1, R2, C1, C2, P = u 
   S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   
   ##compute these once for all... essentially same as in functional responses above 
   num_R1P = S1 * a * o1 * R1
   num_R2P = (1 - S1) * a * o2 * R2
   num_C1P = S1 * a * (1 - o1) * C1
   num_C2P = (1 - S1) * a * (1 - o2) * C2
   denom_RCP =  1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2)
   f_RC1 = a * R1 / (1 + a * h * R1)
   f_RC2 = a * R2 / (1 + a * h * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * num_R1P / denom_RCP
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * num_R2P / denom_RCP
   du[3] = e * C1 * f_RC1 - P * num_C1P / denom_RCP - m * C1
   du[4] = e * C2 * f_RC2 - P * num_C2P / denom_RCP - m * C2 
   du[5] = e * P * (num_R1P + num_R2P) / denom_RCP + e * P * (num_C1P + num_C2P) / denom_RCP - m * P ##not sure if this is the right way to add those functional responses together? 
   
   return du
   end



# Utilities for doing eigenvalue analysis
##this is a right-hand-side function that is used to solve ODEs 
function rhs(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

##set initial condition
u0 = [0.66, 0.66, 0.27, 0.27, 0.15]
tspan = (0.0, 1000.0)

##set Parameters
p = ModelPar()


##Define the ODE problem
prob = ODEProblem(rhs, u0, tspan, p)
sol = solve(prob)

fr_vals = [total_FR_into_P(u, p, t) for (u, t) in zip(sol.u, sol.t)]
##plot timeseries
plot(sol, xlabel="Time", ylabel="Population", title="ODE Solution")


##calculate numerical jacobian
J_num = ForwardDiff.jacobian(x -> rhs(x, p), u0)
println("Jacobian Matrix:")
display(J_num)
##find the equilibrium for this system of equations 
eq = nlsolve((du, u) -> model!(du, u, p, 0.0), u0).zero

# Display equilibrium
println("Equilibrium point (R*, C*): ", eq)

