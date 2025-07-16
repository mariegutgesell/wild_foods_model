##Pulse experiment on G 
##Date Initiated: May 5, 2025
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

# Degree of omnivory functions -- where o represents preference, based on McCann et al., 2005 
function deg_om_i(u, p, t) ##constant speed of switching, where degree of omnivory shifts in response to changing resource/consumer densities in given patch i -  and degree of omnivory used in functional response (can consider later adding in dynamic omnivory)
   R1, C1 = u
    return (p.o * p.e * p.a * R1) / (p.o * p.e * p.a * R1 + (1 - p.o)*p.e * p.a * C1)
end


function deg_om_j(u, p ,t)
    R2, C2 = u
    return (p.o * p.e * p.a * R2) / (p.o * p.e * p.a * R2 + (1 - p.o)*p.e * p.a * C2)
end


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
    G = p.G_func(t)
    return(p.H * G) / (p.H * G + (1-p.H)*R1 + (1-p.H)*R2 + (1-p.H) * C1 + (1-p.H)*C2)
end


##functional response between resources and predator
function f_R1P(u, p, t)
    @unpack a, h= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    ##H is preference for groceries  -holding constant 
    R1, R2, C1, C2 = u  ##defines state variables, G = groceries
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * S1 * a * o1 * R1
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2 + H1 * a * h * G)
    return numerator / denominator
end

function f_R2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * (1 - S1) * a * o2 * R2
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2 + H1 * a * h * G)
    return numerator / denominator
end

function f_C1P(u, p, t)
    @unpack a, h= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * S1 * a * (1 - o1) * C1
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2 + H1 * a * h * G)
    return numerator / denominator
end

function f_C2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = (1-H1) * (1 - S1) * a * (1 - o2) * C2
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2 + H1 * a * h * G)
    return numerator / denominator
end

##trying out functional resposne for G.. type 2 functional response 
function f_GP(u, p, t)
    @unpack a, h= p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    H1 = p.sub_pref(u, p, t)
    G = p.G_func(t)

    numerator = H1 * a * G #trying if i remove scaling / suppression of G by other foraging preferences, i think this makes biological sense (but keep 1-H1 in other FRs)
    denominator = 1 + (S1 * a * h * o1 * R1 + (1-S1)* a * h * o2 * R2 + S1 * a * h * (1-o1) * C1 + (1-S1) * a * h * (1-o2) * C2 + H1 * a * h * G)
    return numerator / denominator
end


##linear (type 1) functional response 
function f_GP_2(u, p, t)
    @unpack a, h, H = p 
    R1, R2, C1, C2 = u 
    G = p.G_func(t)
    return G * H
end

##Perturbation experiment 
##Define time varying G(t)
G_pre = 2.0
G_pulse = 0.25*G_pre
t_pulse = 500.0 ##time when disturbance occurs, want to be once model at equilibirum
t_recover = 550.0 ##time when decline in resources ends 
##set function so before pulse G = G_pre, and during pulse is G_pulse 
G_func(t) = (t < t_pulse || t ≥ t_recover) ? G_pre : G_pulse

#Set up paramters - STRUCTURE 1: Diverse harvest, across multiple habitats 
@with_kw mutable struct ModelPar2_pulse
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 0.6   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.2 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)

    H = 0.5 ##H = preference for groceries (G)

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

    ##time varying function of G
    G_func::Function = t -> 1.0
end

##STRUCTURE 2: low diversity, no coupling, no omnivory
@with_kw mutable struct ModelPar4_pulse
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 1.0   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.0 ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)

    H = 0.5 ##H = preference for groceries (G)

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

     ##time varying function of G
     G_func::Function = t -> 1.0
end


##set up model with G varying w/ time 
function model_2_pulse!(du, u, p ,t)
    @unpack r, K, a, e, h, m, H = p
   R1, R2, C1, C2, P = u 
   #S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
   #o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
   #o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
   G = p.G_func(t) 

   ##predator functional responses
   f_r1p = p.f_r1p(u, p, t)
   f_r2p = p.f_r2p(u, p, t)
   f_c1p = p.f_c1p(u, p, t)
   f_c2p = p.f_c2p(u, p, t)
     f_gp =p.f_gp(u, p, t)
   # f_gp2 = p.f_gp2(u, p, t)
## consumer functional responses
    f_RC1 = a * R1 / (1 + a * h * R1)
   f_RC2 = a * R2 / (1 + a * h * R2)
   
   ##ODEs
   du[1] = r * R1 * (1 - R1 / K) - C1 * f_RC1 - P * f_r1p
   du[2] = r * R2 * (1 - R2 / K) - C2 * f_RC2 - P * f_r2p
   du[3] = e * C1 * f_RC1 - P * f_c1p - m * C1
   du[4] = e * C2 * f_RC2 - P * f_c2p - m * C2 
  # du[5] = e * P * f_r1p + e * P * f_r2p + e * P * f_c1p + e * P * f_c2p + e * P * f_gp2 - m * P
    du[5] = 0 ##can do it this way here if using differential equaitons, but to find equilibrium using NLSOlve will need to use a wrapper function that says not to solve for P and G 

  ##don't treat G as state variable, is a constant, non-dynamic variable that instantly replenishes, have fixed, constant as parameter in functional responses 
   return du
 end 

    ##this is a right-hand-side function that is used to solve ODEs 
    function rhs_pulse(u, p, t=0.0)
        du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
        model_2_pulse!(du, u, p, t) 
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
    @unpack a, h, H = p
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




 ##try plotting for structure 1 
##Solve ODE 
##set initial condition
u0 = [0.6, 0.8, 0.45, 0.61, 0.2]
tspan = (0.0, 1000.0)
##set Parameters
p = ModelPar2_pulse(G_func = G_func)

##Define the ODE problem
prob_3 = ODEProblem(rhs_pulse, u0, tspan, p)
sol_3 = solve(prob_3)

##plot timeseries
   plot(sol_3, xlabel="Time", ylabel="Population", title="ODE Solution w/ Pulse Perturbation on G")
  
##Look at predator consumption  
fr_vals_3 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_3.u, sol_3.t)]

# Separate each series into its own array
fr_total = [x[1] for x in fr_vals_3]
fr_R1 = [x[2] for x in fr_vals_3]
fr_R2 = [x[3] for x in fr_vals_3]
fr_C1 = [x[4] for x in fr_vals_3]
fr_C2 = [x[5] for x in fr_vals_3]
fr_G = [x[6] for x in fr_vals_3]
times = sol_3.t

# Find index where time >= 250
start_idx = findfirst(t -> t ≥ 250, times)

# Subset time and consumption arrays
times_sub = times[start_idx:end]
fr_total_sub = fr_total[start_idx:end]
fr_R1_sub = fr_R1[start_idx:end]
fr_R2_sub = fr_R2[start_idx:end]
fr_C1_sub = fr_C1[start_idx:end]
fr_C2_sub = fr_C2[start_idx:end]
fr_G_sub = fr_G[start_idx:end]

# Define custom colors (optional - adjust as needed)
colors = [:black, :darkgreen, :salmon, :lightgreen, :pink, :red]

# Plot
plot(times_sub, fr_total_sub, label = "total → P", xlabel = "Time", ylabel = "Community Consumption", ylims = (0.0, 1.25), color = colors[1], linewidth = 2.5)
plot!(times_sub, fr_R1_sub, label = "R1 → P", color = colors[2], linewidth = 2.5)
plot!(times_sub, fr_R2_sub, label = "R2 → P", color = colors[3], linewidth = 2.5)
plot!(times_sub, fr_C1_sub, label = "C1 → P", color = colors[4], linewidth = 2.5)
plot!(times_sub, fr_C2_sub, label = "C2 → P", color = colors[5], linewidth = 2.5)
plot!(times_sub, fr_G_sub, label = "G → P", color = colors[6], linewidth = 2.5)


   ##plot total P consumption 
   plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Community Consumption",ylims=(0.0,1.25))  
   
   plot!(times, fr_R1, label = "R1 → P")
plot!(times, fr_R2, label = "R2 → P")
plot!(times, fr_C1, label = "C1 → P")
plot!(times, fr_C2, label = "C2 → P")
plot!(times, fr_G, label = "G → P")



 ##try plotting for structure 2 
##Solve ODE 
##set initial condition
u0 = [0.6, 0.8, 0.45, 0.61, 0.2]
tspan = (0.0, 1000.0)
##set Parameters
p = ModelPar4_pulse(G_func = G_func)

##Define the ODE problem
prob_4 = ODEProblem(rhs_pulse, u0, tspan, p)
sol_4 = solve(prob_4)

##plot timeseries
   plot(sol_4, xlabel="Time", ylabel="Population", title="ODE Solution w/ Pulse Perturbation on G")
  
##Look at predator consumption  
fr_vals_4 = [total_FR_into_P(u, p, t) for (u, t) in zip(sol_4.u, sol_4.t)]

# Separate each series into its own array
fr_total_4 = [x[1] for x in fr_vals_4]
fr_R1_4 = [x[2] for x in fr_vals_4]
fr_R2_4 = [x[3] for x in fr_vals_4]
fr_C1_4 = [x[4] for x in fr_vals_4]
fr_C2_4 = [x[5] for x in fr_vals_4]
fr_G_4 = [x[6] for x in fr_vals_4]
times_4 = sol_4.t

   ##plot total P consumption 
   plot(times_4, fr_total_4, label = "total → P", xlabel = "Time", ylabel = "Community Consumption",
   ylims=(0.0,1.25))  
   plot!(times_4, fr_R1_4, label = "R1 → P")
plot!(times_4, fr_R2_4, label = "R2 → P")
plot!(times_4, fr_C1_4, label = "C1 → P")
plot!(times_4, fr_C2_4, label = "C2 → P")
plot!(times_4, fr_G_4, label = "G → P")


# Find index where time >= 250
start_idx = findfirst(t -> t ≥ 250, times_4)

# Subset time and consumption arrays
times_sub_4 = times_4[start_idx:end]
fr_total_sub = fr_total_4[start_idx:end]
fr_R1_sub = fr_R1_4[start_idx:end]
fr_R2_sub = fr_R2_4[start_idx:end]
fr_C1_sub = fr_C1_4[start_idx:end]
fr_C2_sub = fr_C2_4[start_idx:end]
fr_G_sub = fr_G_4[start_idx:end]

# Define custom colors (optional - adjust as needed)
colors = [:black, :darkgreen, :salmon, :lightgreen, :pink, :red]

# Plot
plot(times_sub_4, fr_total_sub, label = "total → P", xlabel = "Time", ylabel = "Community Consumption", ylims = (0.0, 1.25), color = colors[1], linewidth = 2.5)
plot!(times_sub_4, fr_R1_sub, label = "R1 → P", color = colors[2],linewidth = 2.5)
plot!(times_sub_4, fr_R2_sub, label = "R2 → P", color = colors[3],linewidth = 2.5)
plot!(times_sub_4, fr_C1_sub, label = "C1 → P", color = colors[4],linewidth = 2.5)
plot!(times_sub_4, fr_C2_sub, label = "C2 → P", color = colors[5],linewidth = 2.5)
plot!(times_sub_4, fr_G_sub, label = "G → P", color = colors[6],linewidth = 2.5)


##Plotting both results from both model types
plot(times, fr_total, label = "total → P", xlabel = "Time", ylabel = "Predator Consumption")  
 