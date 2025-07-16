##Code to look at functional responses 
##Date Initiated: April 8, 2025
##Contributor(s): Marie K. Gutgesell

##source model - choose which based on which you want to investigate FR
#include("wf_model_eqs_subsidy.jl")

#include("wf_model_eqs_subsidy_TLparams.jl")

include("wf_model_eqs_subsidy_2.jl") ##model equations with P NOT held constant, unique parameters per trophic level, active and passive omnivory parameter structures
##Plotting functional responses across range of R/C
p = ModelPar_active()

##Plot P functional response across range of R1, holding all other variables constant
R1_vals = range(0, 50, length = 100)
response_vals = [f_R1P((R1, 10.0, 5.0, 5.0, 3.0), p, 0.0) for R1 in R1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R1_vals, response_vals,
xlabel = "R1 Density", ylabel = "Predator Consumption Rate of R1", title = "Functional Response to R1")

##Plot C functional response across range of R1, holding all other variables constant
response_vals_C1 = [f_R1C1((R1, 10.0, 5.0, 5.0, 3.0), p, 0.0) for R1 in R1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R1_vals, response_vals_C1,
xlabel = "R1 Density", ylabel = "C1 Consumption Rate of R1", title = "C1 Functional Response to R1")


##Plot P functional response across range of R2, holding all other variables constant
R2_vals = range(0, 100, length = 100)
response_vals_R2 = [f_R2P((10.0, R2, 5.0, 5.0, 3.0), p, 0.0) for R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R2_vals, response_vals_R2,
xlabel = "R2 Density", ylabel = "Predator Consumption Rate of R2", title = "Functional Response to R2")


##Plot P functional response across range of R1, holding all other variables constant
C1_vals = range(0, 50, length = 100)
response_vals_C1 = [f_C1P((10.0, 10.0, C1, 5.0, 3.0), p, 0.0) for C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C1_vals, response_vals_C1,
xlabel = "C1 Density", ylabel = "Predator Consumption Rate of C1", title = "Functional Response to C1")



##Plot P functional response across range of C2, holding all other variables constant
C2_vals = range(0, 100, length = 100)
response_vals_C2 = [f_C2P((10.0, 10.0, 5.0, C2, 3.0), p, 0.0) for C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C2_vals, response_vals_C2,
xlabel = "C2 Density", ylabel = "Predator Consumption Rate of C2", title = "Functional Response to C2")

##Plot P functional response across range of G, holding all other variables constant
G_vals = range(0, 100, length = 100)
response_vals_G = []

for G in G_vals
    p.G = G  # update the parameter G for this iteration
    push!(response_vals_G, f_GP((10.0, 10.0, 5.0, 5.0, 3.0), p, 0.0))
end

plot(G_vals,response_vals_G,
xlabel = "G Density", ylabel = "Predator Consumption Rate of G", title = "Functional Response to G")
##now getting basically type 2 / type 3 functional response... is that right? yes i think so.. 


##Looking into S1, o1 and o2, to see if they are responding how I think they should be across C/R densities
C1_vals = range(0.1, 100, length = 100)
W1_vals_C1 = [p.W((10.0, 10.0, C1, 5.0, 3.0), p, 0.0) for C1 in C1_vals]
plot(C1_vals, W1_vals_C1, xlabel = "C1 Density", ylabel = "W1 (foraging in Patch 1)", title = "Habitat Preference (W1) vs. C1 Density")

C2_vals = range(0.1, 100, length = 100)
W1_vals_C2 = [p.W((10.0, 10.0, 5.0, C2, 3.0), p, 0.0) for C2 in C2_vals]
plot(C2_vals, W1_vals_C2, xlabel = "C2 Density", ylabel = "W1 (foraging in Patch 1)", title = "Foraging Preference (W1) vs. C2 Density")


R1_vals = range(0.1, 100, length = 100)
W1_vals_R1 = [p.W((R1, 10.0, 5.0, 5.0, 3.0), p, 0.0) for R1 in R1_vals]
plot(R1_vals, W1_vals_R1, xlabel = "R1 Density", ylabel = "S1 (foraging in Patch 1)", title = "Foraging Preference (S1) vs. R1 Density")
##so habitat preference is driven by differences in C densities, so P needs to make a choice, and R density has no effect on foraging 


##Looking at o1
o1_vals_R1 = [p.d_om_i((R1, 10.0), p, 0.0) for R1 in R1_vals]
plot(R1_vals, o1_vals_R1, xlabel = "R1 Density", ylabel = "Degree of Omnivory in Patch i")
##increasing o (higher preference/speed of switching for R) makes degree of omnivory increase rapidly, so working as expected - density dependent 
##note: if degree of omnivory 

o1_vals_C1 = [p.d_om_i((10.0, C1), p, 0.0) for C1 in C1_vals]
plot(C1_vals, o1_vals_C1, xlabel = "C1 Density", ylabel = "Degree of Omnivory in Patch i")



o1_vals_R2 = [p.d_om_i((10.0, 5.0), p, 0.0) for R2 in R2_vals]
plot(R2_vals, o1_vals_R2, xlabel = "R2 Density", ylabel = "Degree of Omnivory in Patch i")

##Looking at H values (preference for groceries)


H_response_vals_G = []

for G in G_vals
    p.G = G  # update the parameter G for this iteration
    push!(H_response_vals_G, p.sub_pref((10.0, 10.0, 5.0, 5.0, 3.0), p, 0.0))
end
plot(G_vals, H_response_vals_G, xlabel = "G Density", ylabel = "H (foraging on groceries)", title = "Grocery Preference Based on G")

H_vals_C1 = [p.sub_pref((10.0, 10.0, C1, 5.0, 5.0, 3.0), p, 0.0) for C1 in C1_vals]
plot(C1_vals, H_vals_C1, xlabel = "C1 Density", ylabel = "H (foraging on groceries)", title = "Grocery Preference Based on G")

