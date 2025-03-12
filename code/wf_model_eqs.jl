##Code to set up and explore equations for wild food model
##Date Initiated: October 3, 2024
##Contributor(s): Marie K. Gutgesell

using Pkg
Pkg.add("Plots")
##Load Libraries needed (ODE for later)
using Parameters: @with_kw, @unpack ##imports Parameters package that provides convenient macros for working with keyword arugemnts, parameter structs and unpacking variables 
using DifferentialEquations
using Plots

##Re-creating model from McCann et al., 2005, Ecology Letters
#u = state variables where u[1] = R1, u[2] = R2, u[3] = C1, u[4] = C2, u[5] = P
#p = Parameters
#t = 

##Habitat preference (coupling)
##habitat preference - so preference for habitat changes with shifting densities of C across patches (but preference w itself is fixed)
function hab_pref(u, p, t)
    C1, C2 = u
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
    return p.o * (p.e * p.a * R1 / (p.o * p.e * p.a * R1) + (1 - o)*p.e * p.a * C1)
end
##know the R/C isnt right... 
##hmm how do i do this so the state variables change depending on if i or j patch? - come back to this
##

function deg_om_j(u, p ,t)
    R2, C2 = u
    return p.o * (p.e * p.a * R2 / (p.o * p.e * p.a * R2) + (1 - o)*p.e * p.a * C2)
end

#Set up paramters 
@with_kw mutable struct ModelPar
    ##spatial scale of foraging 
    Sc = 1  ##Sc = spatial scale of consumer foraging
    Sh = 1  ##Sh = spatial scale of resources 
    Q = 0 ##for now starting with Q = 0 (Q = (Sc-Sh)/Sh) - #Q is the area of consumers current habitat where consumer percieves prey as well mixed; so 1-Q is area that requres consumers to switch patches in order to gain access to resources in that portion of the patch
    w = 0.85   ##w = habitat preference for patch i 
        ##do we need separate w for each C/R? wR1, wR2, wC1, wC2 ? 
    o = 0.5  ##o = omnivory preference, preference for either consumer or resource, starting with fixed omnivory preference (essentially passive case)


    ##Model parameters, for now just keeping these parameters the same (could make unique ones based on patch/TL)
    r = 1.0
    K = 3.25
    a = 2.5  ##attack rate 
    e = 0.8   ##energy conversion 
    m = 1.0
    h = 0.5
    ##Habitat preference function 
    h_pref::Function hab_pref  

    ##Foraging scale
    S::Function forage_scale
    
    ##Omnivory preference function
    d_om_i::Function deg_om_i ##if d_om = 0, have food chain (no omnivory) in patch 1
    d_om_j::Function deg_om_j ##degree of omnivory in patch 2 
    
   


end



##functional response between resources and predator
##hm separate functional response for R1 and R2... i think so... 
function f_R1P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    return (S1 * a * o1 * R1)  / 1 + ((S1 * a * h * o1 * R1) + ((1-S1)* a * h * o2 * R2)) + ((S1 * a * h * (1-o1) * C1) + ((1-S1) * a * h * (1-o2) * C2))
end

function f_R2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    return ((1 - S1) * a * o2 * R2)  / 1 + ((S1 * a * h * o1 * R1) + ((1-S1)* a * h * o2 * R2)) + ((S1 * a * h * (1-o1) * C1) + ((1-S1) * a * h * (1-o2) * C2))
end

function f_C1P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    return (S1 * a * (1 - o1) * C1)  / 1 + ((S1 * a * h * o1 * R1) + ((1-S1)* a * h * o2 * R2)) + ((S1 * a * h * (1-o1) * C1) + ((1-S1) * a * h * (1-o2) * C2))
end

function f_C2P(u, p, t)
    @unpack a, h = p ##note: may need to put the p directly in the equation, not sure if will work if calling p for the functions below
    R1, R2, C1, C2 = u  ##defines state variables
    S1 = p.S(u, p, t) ##function that defines foraging scale and habitat preference (since Q = 0, Si = Wi, so Sj = 1-Si)
    o1 = p.d_om_i(u, p, t) ##function that defines degree of omnivory in patch 1
    o2 = p.d_om_j(u, p, t) ##function that defines degree of omnivory in patch 2 
    
    return ((1 - S1) * a * (1 - o2) * C2)  / 1 + ((S1 * a * h * o1 * R1) + ((1-S1)* a * h * o2 * R2)) + ((S1 * a * h * (1-o1) * C1) + ((1-S1) * a * h * (1-o2) * C2))
end


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
denom_RCP =  1 + ((S1 * a * h * o1 * R1) + ((1-S1)* a * h * o2 * R2)) + ((S1 * a * h * (1-o1) * C1) + ((1-S1) * a * h * (1-o2) * C2))
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


