##Code to look at functional responses and preference functions 
##Date Initiated: October 8, 2025
##Contributor(s): Marie K. Gutgesell

##source model - choose which based on which you want to investigate FR


include("wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P held constant, unique parameters per trophic level, active and passive omnivory parameter structures

##Plotting out preference functions
p = ModelPar_active()
u = [1.0, 1.0, 2.0, 1.0, 0.25]
t = 0.0
##omnivory preference function - active
o_vals = 0.0:0.1:1.0
pref_vals = [(p.o = o; om_i_pref_active(u, p, t)) for o in o_vals]
print(pref_vals)
plot(o_vals, pref_vals, xlabel = "o", ylabel = "omnivory preference (active)", title = "Omnivory Preference over o")

##omnivory preference function - active
o_vals = 0.0:0.1:1.0
pref_vals = [(p.o = o; om_i_pref_fixed(u, p, t)) for o in o_vals]
print(pref_vals)
plot(o_vals, pref_vals, xlabel = "o", ylabel = "omnivory preference (passive)", title = "Omnivory Preference over o")

##what is the difference here between passive and active? 


##habitat preference function
w_vals = 0.0:0.1:1.0
pref_vals = [(p.w = w; hab_pref(u, p, t)) for w in w_vals]
print(pref_vals)
plot(w_vals, pref_vals, xlabel = "w", ylabel = "habitat preference", title = "Habitat Preference over w")

##grocery preference 
H_vals = 0.0:0.1:1.0
pref_vals = [(p.H = H; sub_pref_func(u, p, t)) for H in H_vals]
print(pref_vals)
plot(H_vals, pref_vals, xlabel = "H", ylabel = "Grocery preference", title = "Grocery Preference over H")

##get nonlinearities in preference functions when densities are not equal 

###Plotting functional responses across ranges of o, w, and h
p = ModelPar_active()
u = [1.0, 1.0, 1.0, 1.0, 0.25]
t = 0.0

##functional responses across o 
o_vals = 0.0:0.1:1.0

##R1_P
FR1_P_vals = [(p.o = o; f_R1P(u, p, t)) for o in o_vals]
print(FR1_P_vals)
plot(o_vals, FR1_P_vals, xlabel = "o", ylabel = "FR1_P", title = "FR1_P functional response over o")

##R2_P
FR2_P_vals = [(p.o = o; f_R2P(u, p, t)) for o in o_vals]
print(FR2_P_vals)
plot(o_vals, FR2_P_vals, xlabel = "o", ylabel = "FR2_P", title = "FR2_P functional response over o")

##C1_P
FC1_P_vals = [(p.o = o; f_C1P(u, p, t)) for o in o_vals]
print(FC1_P_vals)
plot(o_vals, FC1_P_vals, xlabel = "o", ylabel = "FC1_P", title = "FC1_P functional response over o")

##C2_P
FC2_P_vals = [(p.o = o; f_C2P(u, p, t)) for o in o_vals]
print(FC2_P_vals)
plot(o_vals, FC2_P_vals, xlabel = "o", ylabel = "FC2_P", title = "FC2_P functional response over o")

##plot all together 
plot(o_vals, FR1_P_vals, xlabel = "o", ylabel = "functional response", title = "functional response over o", label = "FR1_P")
plot!(o_vals, FR2_P_vals, xlabel = "o", label = "FR2_P")
plot!(o_vals, FC1_P_vals, xlabel = "o", label = "FC1_P")
plot!(o_vals, FC2_P_vals, xlabel = "o", label = "FC2_P")


##Range of w 
w_vals = 0.0:0.1:1.0

##R1_P
FR1_P_vals = [(p.w = w; f_R1P(u, p, t)) for w in w_vals]
print(FR1_P_vals)
plot(w_vals, FR1_P_vals, xlabel = "w", ylabel = "FR1_P", title = "FR1_P functional response over w")

##R2_P
FR2_P_vals = [(p.w = w; f_R2P(u, p, t)) for w in w_vals]
print(FR2_P_vals)
plot(w_vals, FR2_P_vals, xlabel = "w", ylabel = "FR2_P", title = "FR2_P functional response over w")

##C1_P
FC1_P_vals = [f_C1P(u, ModelPar_active(; w=w, W=hab_pref), t) for w in w_vals]
print(FC1_P_vals)
plot(w_vals, FC1_P_vals, xlabel = "w", ylabel = "FC1_P", title = "FC1_P functional response over w")

##C2_P
FC2_P_vals = [f_C2P(u, ModelPar_active(; w=w, W=hab_pref), t) for w in w_vals]
print(FC2_P_vals)

plot(w_vals, FC2_P_vals, xlabel = "w", ylabel = "FC2_P", title = "FC2_P functional response over w")

##plot of all together:
plot(w_vals, FR1_P_vals, xlabel = "w", ylabel = "functional response", title = "functional response over w", label = "FR1_P")
plot!(w_vals, FR2_P_vals, xlabel = "w", label = "FR2_P")
plot!(w_vals, FC1_P_vals, xlabel = "w", label = "FC1_P")
plot!(w_vals, FC2_P_vals, xlabel = "w", label = "FC2_P")


##Range of H
H_vals = 0.0:0.1:1.0

##R1_P
FR1_P_vals = [(p.H = H; f_R1P(u, p, t)) for H in H_vals]
print(FR1_P_vals)
plot(H_vals, FR1_P_vals, xlabel = "H", ylabel = "FR1_P", title = "FR1_P functional response over H")

##R2_P
FR2_P_vals = [(p.H = H; f_R2P(u, p, t)) for H in H_vals]
print(FR2_P_vals)
plot(H_vals, FR2_P_vals, xlabel = "H", ylabel = "FR2_P", title = "FR2_P functional response over H")

##C1_P
FC1_P_vals = [f_C1P(u, ModelPar_active(; H=H), t) for H in H_vals]
print(FC1_P_vals)
plot(H_vals, FC1_P_vals, xlabel = "H", ylabel = "FC1_P", title = "FC1_P functional response over H")

##C2_P
FC2_P_vals = [f_C2P(u, ModelPar_active(; H = H), t) for H in H_vals]
print(FC2_P_vals)
plot(H_vals, FC2_P_vals, xlabel = "H", ylabel = "FC2_P", title = "FC2_P functional response over H")

##plot of all together:
plot(H_vals, FR1_P_vals, xlabel = "H", ylabel = "functional response", title = "functional response over H", label = "FR1_P")
plot!(H_vals, FR2_P_vals, xlabel = "H", label = "FR2_P")
plot!(H_vals, FC1_P_vals, xlabel = "H", label = "FC1_P")
plot!(H_vals, FC2_P_vals, xlabel = "H", label = "FC2_P")





##Plotting functional responses across range of R/C
p = ModelPar_active()

##Plot P functional response across range of R1, holding all other variables constant
R1_vals = range(0, 15, length = 10)
R2_vals = range(0, 15, length = 10)
C1_vals = range(0, 15, length = 10)
C2_vals = range(0, 15, length = 10)

response_vals_R1 = [f_R1P((R1, R1, 1.0, 1.0, 0.25), p, 0.0) for R1 in R1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
plot(R1_vals, response_vals_R1,
xlabel = "R1 Density", ylabel = "Predator Consumption Rate of R1", title = "Functional Response to R1")

response_vals_R1_R2 = [f_R1P((R1, R2, 1.0, 1.0, 0.25), p, 0.0) for R1 in R1_vals, R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, R2_vals, response_vals_R1_R2;
    xlabel = "R1 density",
    ylabel = "R2 density",
    title  = "f_R1P across R1–R2",
    colorbar_title = "Predator consumption of R1"
)

response_vals_R1_C1 = [f_R1P((R1, 1.0, C1, 1.0, 0.25), p, 0.0) for R1 in R1_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, C1_vals, response_vals_R1_C1;
    xlabel = "R1 density",
    ylabel = "C1 density",
    title  = "f_R1P across R1–C1",
    colorbar_title = "Predator consumption of R1"
)

response_vals_R2_C1 = [f_R1P((1.0, R2, C1, 1.0, 0.25), p, 0.0) for R2 in R2_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C1_vals, response_vals_R2_C1;
    xlabel = "R2 density",
    ylabel = "C1 density",
    title  = "f_R1P across R2-C1",
    colorbar_title = "Predator consumption of R1"
)

response_vals_R2_C2 = [f_R1P((1.0, R2, 1.0, C2, 0.25), p, 0.0) for R2 in R2_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C2_vals, response_vals_R2_C2;
    xlabel = "R2 density",
    ylabel = "C2 density",
    title  = "f_R1P across R2-C2",
    colorbar_title = "Predator consumption of R1"
)

response_vals_C1_C2 = [f_R1P((1.0, 1.0, C1, C2, 0.25), p, 0.0) for C1 in C1_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    C1_vals, C2_vals, response_vals_C1_C2;
    xlabel = "C1 density",
    ylabel = "C2 density",
    title  = "f_R1P across C1-C2",
    colorbar_title = "Predator consumption of R1"
)




##Plot P functional response across range of R2, holding all other variables constant
R2_vals = range(0, 10, length = 10)
response_vals_R2 = [f_R2P((1.0, R2, 1.0, 1.0, 0.5), p, 0.0) for R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(R2_vals, response_vals_R2,
xlabel = "R2 Density", ylabel = "Predator Consumption Rate of R2", title = "Functional Response to R2")

response_vals_R1_R2 = [f_R2P((R1, R2, 1.0, 1.0, 0.25), p, 0.0) for R1 in R1_vals, R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, R2_vals, response_vals_R1_R2;
    xlabel = "R1 density",
    ylabel = "R2 density",
    title  = "f_R2P across R1–R2",
    colorbar_title = "Predator consumption of R2"
)

response_vals_R1_C1 = [f_R2P((R1, 1.0, C1, 1.0, 0.25), p, 0.0) for R1 in R1_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, C1_vals, response_vals_R1_C1;
    xlabel = "R1 density",
    ylabel = "C1 density",
    title  = "f_R2P across R1–C1",
    colorbar_title = "Predator consumption of R2"
)

response_vals_R2_C1 = [f_R2P((1.0, R2, C1, 1.0, 0.25), p, 0.0) for R2 in R2_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C1_vals, response_vals_R2_C1;
    xlabel = "R2 density",
    ylabel = "C1 density",
    title  = "f_R2P across R2-C1",
    colorbar_title = "Predator consumption of R2"
)

response_vals_R2_C2 = [f_R2P((1.0, R2, 1.0, C2, 0.25), p, 0.0) for R2 in R2_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C2_vals, response_vals_R2_C2;
    xlabel = "R2 density",
    ylabel = "C2 density",
    title  = "f_R2P across R2-C2",
    colorbar_title = "Predator consumption of R2"
)

response_vals_C1_C2 = [f_R2P((1.0, 1.0, C1, C2, 0.25), p, 0.0) for C1 in C1_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    C1_vals, C2_vals, response_vals_C1_C2;
    xlabel = "C1 density",
    ylabel = "C2 density",
    title  = "f_R2P across C1-C2",
    colorbar_title = "Predator consumption of R2"
)



##Plot P functional response across range of R1, holding all other variables constant
#C1_vals = range(0, 10, length = 10)
response_vals_C1 = [f_C1P((1.0, 1.0, C1, 1.0, 0.25), p, 0.0) for C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C1_vals, response_vals_C1,
xlabel = "C1 Density", ylabel = "Predator Consumption Rate of C1", title = "Functional Response to C1")

response_vals_R1_R2 = [f_C1P((R1, R2, 1.0, 1.0, 0.25), p, 0.0) for R1 in R1_vals, R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, R2_vals, response_vals_R1_R2;
    xlabel = "R1 density",
    ylabel = "R2 density",
    title  = "f_C1P across R1–R2",
    colorbar_title = "Predator consumption of C1"
)

response_vals_R1_C1 = [f_C1P((R1, 1.0, C1, 1.0, 0.25), p, 0.0) for R1 in R1_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, C1_vals, response_vals_R1_C1;
    xlabel = "R1 density",
    ylabel = "C1 density",
    title  = "f_C1P across R1–C1",
    colorbar_title = "Predator consumption of C1"
)

response_vals_R2_C1 = [f_C1P((1.0, R2, C1, 1.0, 0.25), p, 0.0) for R2 in R2_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C1_vals, response_vals_R2_C1;
    xlabel = "R2 density",
    ylabel = "C1 density",
    title  = "f_C1P across R2-C1",
    colorbar_title = "Predator consumption of C1"
)

response_vals_R2_C2 = [f_C1P((1.0, R2, 1.0, C2, 0.25), p, 0.0) for R2 in R2_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C2_vals, response_vals_R2_C2;
    xlabel = "R2 density",
    ylabel = "C2 density",
    title  = "f_C1P across R2-C2",
    colorbar_title = "Predator consumption of C1"
)

response_vals_C1_C2 = [f_C1P((1.0, 1.0, C1, C2, 0.25), p, 0.0) for C1 in C1_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    C1_vals, C2_vals, response_vals_C1_C2;
    xlabel = "C1 density",
    ylabel = "C2 density",
    title  = "f_C1P across C1-C2",
    colorbar_title = "Predator consumption of C1"
)




##Plot P functional response across range of C2, holding all other variables constant
#C2_vals = range(0, 10, length = 10)
response_vals_C2 = [f_C2P((1.0, 1.0, 1.0, C2, 0.25), p, 0.0) for C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2

plot(C2_vals, response_vals_C2,
xlabel = "C2 Density", ylabel = "Predator Consumption Rate of C2", title = "Functional Response to C2")

response_vals_R1_R2 = [f_C2P((R1, R2, 1.0, 1.0, 0.25), p, 0.0) for R1 in R1_vals, R2 in R2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, R2_vals, response_vals_R1_R2;
    xlabel = "R1 density",
    ylabel = "R2 density",
    title  = "f_C2P across R1–R2",
    colorbar_title = "Predator consumption of C2"
)

response_vals_R1_C1 = [f_C2P((R1, 1.0, C1, 1.0, 0.25), p, 0.0) for R1 in R1_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R1_vals, C1_vals, response_vals_R1_C1;
    xlabel = "R1 density",
    ylabel = "C1 density",
    title  = "f_C2P across R1–C1",
    colorbar_title = "Predator consumption of C2"
)

response_vals_R2_C1 = [f_C2P((1.0, R2, C1, 1.0, 0.25), p, 0.0) for R2 in R2_vals, C1 in C1_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C1_vals, response_vals_R2_C1;
    xlabel = "R2 density",
    ylabel = "C1 density",
    title  = "f_C2P across R2-C1",
    colorbar_title = "Predator consumption of C2"
)

response_vals_R2_C2 = [f_C2P((1.0, R2, 1.0, C2, 0.25), p, 0.0) for R2 in R2_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    R2_vals, C2_vals, response_vals_R2_C2;
    xlabel = "R2 density",
    ylabel = "C2 density",
    title  = "f_C2P across R2-C2",
    colorbar_title = "Predator consumption of C2"
)

response_vals_C1_C2 = [f_C2P((1.0, 1.0, C1, C2, 0.25), p, 0.0) for C1 in C1_vals, C2 in C2_vals] ##the numbers in u part of function set densities for R2, C1 and C2
surface(
    C1_vals, C2_vals, response_vals_C1_C2;
    xlabel = "C1 density",
    ylabel = "C2 density",
    title  = "f_C2P across C1-C2",
    colorbar_title = "Predator consumption of C2"
)



##Plot P functional response across range of G, holding all other variables constant
G_vals = range(0, 10, length = 10)
response_vals_G = []

for G in G_vals
    p.G = G  # update the parameter G for this iteration
    push!(response_vals_G, f_GP((1.0, 1.0, 1.0, 1.0, 0.25), p, 0.0))
end

plot(G_vals,response_vals_G,
xlabel = "G Density", ylabel = "Predator Consumption Rate of G", title = "Functional Response to G")
