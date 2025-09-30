##Looking at stability across gradient of coupling and omnivory 

##source model 
#include("wf_model_eqs_subsidy.jl")

include("wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P constant/or not (depending on which equation on model structure is silenced), unique parameters per trophic level, active and passive omnivory parameter structures
##WHERE LEFT OFF (JULY 3): trying to understand if dynamics from simpler to more complex model match what i would expect based on theory - working through this


##1) Looking at local stability for unforced model to see changes w/ increasing K
p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0)
#P_fixed = 0.25
K_results = equilibrium_forced(p, P0)

results_K_all_uf = []
#P_fixed = 0.25
for K in 0.1:0.1:6.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, K=K)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_K_all_uf, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all_uf)
##in unforced model, getting inf/NAs when K = 7.3, 9.2, 8.9, 7.4 - when evaluating at t = 100
##in unforced model, getting inf/NAs when K = 0.1-0.7, 7.3, 8.0, 9.0 
##if getting different K values when evaluating at different time points, i think not at equilibrium at t = 100? 

##at k = 1.6, finally get persistence of C1 
##at k = 2.0, get persistence of C2 as well -- and monotonic approaches to equilibrium
##at k = 3.0 and 3.1 start to get tiny wiggles -- potentially still part of transient? could be really long .. 
##when go over 100,000 time steps, still getting same pattern -- not stable limit cycles but does look potentially like it is repeating itself.. 
##then at 3.2 wiggles seem to disappear again 
##then 3.3 get tiny starts of potential oscillations, and slight oscillatory decay
##that does make sense, because that is at bottom of checkmark when start to get imaginary part (i think)
##but it isn't going to a straight stable equilibrium .. well R2 and C2 do, but R1 and C1 don't
##looks like at about 3.9/4 thats when start to get some oscillation in R2/C2
##at k = 5, getting stable limit cycles (i think) -- plot the max/mins after this to see if getting bifurcation 
##at k-5.2 start to get much larger oscillations -- 
##at 6.2 start to get different dynamics at later time periods, so definitely not in a stable oscillation -- 
##woah yea crazy shit going on -- longer wild cycles 
##at k = 7 start to not get repeating patterns 

plot(df_eq.K, df_eq.λ1, xlabel = "K", ylabel = "max eigenvalue")
##why do the eigenvalues go crazy like that? - outside of numerical realm or something.. 


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_cv, label = "R1", xlabel = "K", ylabel = "CV")
plot!(df_eq.K, R2_cv, col = "red", label = "R2")
plot!(df_eq.K, C1_cv, col = "blue", label = "C1")
plot!(df_eq.K, C2_cv, col = "green", label = "C2")
plot(df_eq.K, P_cv)

##spike in C1/C2 CV odd, is before they cross 0, so what is driving that? 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_mean, xlabel = "K", ylabel = "Mean")
plot!(df_eq.K, R2_mean)
plot!(df_eq.K, C1_mean)
plot!(df_eq.K, C2_mean)
plot!(df_eq.K, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_sd, xlabel = "K", ylabel = "SD")
plot!(df_eq.K, R2_sd)
plot!(df_eq.K, C1_sd)
plot!(df_eq.K, C2_sd)
plot(df_eq.K, P_sd)

##so means go up, and SD goes up, so increase in mean i guess is driving decline in CV? 
##Look at min/max Plots - bifurcations
plot(df_eq.K, R1_min, label = "R1 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, R1_max, label = "R1 max")

plot(df_eq.K, R2_min, col = "red", label = "R2 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, R2_max, col = "red", label = "R2 max")

plot(df_eq.K, C1_min, col = "blue", label = "C1 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, C1_max, col = "blue", label = "C1 max")

plot(df_eq.K, C2_min, col = "green", label = "C2 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, C2_max, col = "green", label = "C2 max")

plot(df_eq.K, P_cv)



##forced model

results_K_all_f = []
for K in 1.5:0.1:6.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, K=K, l1 = 1.0, l2 = 1.0)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_K_all_f, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all_f)
#print(df_eq)
plot(df_eq.K, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
print(R1_min)
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.K, R1_cv, label = "R1", xlabel = "K", ylabel = "CV")
plot!(df_eq.K, R2_cv, col = "red", label = "R2")
plot!(df_eq.K, C1_cv, col = "blue", label = "C1")
plot!(df_eq.K, C2_cv, col = "green", label = "C2")
plot(df_eq.K, P_cv)

##spike in C1/C2 CV odd, is before they cross 0, so what is driving that? 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_mean, xlabel = "K", ylabel = "Mean")
plot!(df_eq.K, R2_mean)
plot!(df_eq.K, C1_mean)
plot!(df_eq.K, C2_mean)
plot!(df_eq.K, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_sd, xlabel = "K", ylabel = "SD")
plot!(df_eq.K, R2_sd)
plot!(df_eq.K, C1_sd)
plot!(df_eq.K, C2_sd)
plot(df_eq.K, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.K, R1_min, label = "R1 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, R1_max, label = "R1 max")

plot(df_eq.K, R2_min, col = "red", label = "R2 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, R2_max, col = "red", label = "R2 max")

plot(df_eq.K, C1_min, col = "blue", label = "C1 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, C1_max, col = "blue", label = "C1 max")

plot(df_eq.K, C2_min, col = "green", label = "C2 min", xlabel = "K", ylabel ="min/max")
plot!(df_eq.K, C2_max, col = "green", label = "C2 max")

plot(df_eq.K, P_cv)




##how does a influence? - unforced model
results_a_all_uf = []
for aC_P in 0.0:0.1:20.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, aC_P=aC_P, K = 2.75)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_a_all_uf, (; aC_P=aC_P, eq_data...))
end
df_eq = DataFrame(results_a_all_uf)

plot(df_eq.aC_P, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.aC_P, R1_cv, label = "R1", xlabel = "aC_P", ylabel = "CV")
plot!(df_eq.aC_P, R2_cv, col = "red", label = "R2")
plot!(df_eq.aC_P, C1_cv, col = "blue", label = "C1")
plot!(df_eq.aC_P, C2_cv, col = "green", label = "C2")
plot(df_eq.aC_P, P_cv)
##CV dynamics look very odd, all very low, 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.aC_P, R1_mean, xlabel = "aC_P", ylabel = "Mean")
plot!(df_eq.aC_P, R2_mean)
plot!(df_eq.aC_P, C1_mean)
plot!(df_eq.aC_P, C2_mean)
plot!(df_eq.aC_P, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.aC_P, R1_sd, xlabel = "aC_P", ylabel = "SD")
plot!(df_eq.aC_P, R2_sd)
plot!(df_eq.aC_P, C1_sd)
plot!(df_eq.aC_P, C2_sd)
plot(df_eq.aC_P, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.aC_P, R1_min, label = "R1 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, R1_max, label = "R1 max")

plot(df_eq.aC_P, R2_min, col = "red", label = "R2 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, R2_max, col = "red", label = "R2 max")

plot(df_eq.aC_P, C1_min, col = "blue", label = "C1 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, C1_max, col = "blue", label = "C1 max")

plot(df_eq.aC_P, C2_min, col = "green", label = "C2 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, C2_max, col = "green", label = "C2 max")

plot(df_eq.aC_P, P_cv)

##influence of a in forced model
results_a_all_f = []
for aC_P in 0.0:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, aC_P=aC_P)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_a_all_f, (; aC_P=aC_P, eq_data...))
end
df_eq = DataFrame(results_a_all_f)

plot(df_eq.aC_P, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.aC_P, R1_cv, label = "R1", xlabel = "aC_P", ylabel = "CV")
plot!(df_eq.aC_P, R2_cv, col = "red", label = "R2")
plot!(df_eq.aC_P, C1_cv, col = "blue", label = "C1")
plot!(df_eq.aC_P, C2_cv, col = "green", label = "C2")
plot(df_eq.aC_P, P_cv)

##spike in C1/C2 CV odd, is before they cross 0, so what is driving that? 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.aC_P, R1_mean, xlabel = "aC_P", ylabel = "Mean")
plot!(df_eq.aC_P, R2_mean)
plot!(df_eq.aC_P, C1_mean)
plot!(df_eq.aC_P, C2_mean)
plot!(df_eq.aC_P, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.aC_P, R1_sd, xlabel = "aC_P", ylabel = "SD")
plot!(df_eq.aC_P, R2_sd)
plot!(df_eq.aC_P, C1_sd)
plot!(df_eq.aC_P, C2_sd)
plot(df_eq.aC_P, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.aC_P, R1_min, label = "R1 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, R1_max, label = "R1 max")

plot(df_eq.aC_P, R2_min, col = "red", label = "R2 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, R2_max, col = "red", label = "R2 max")

plot(df_eq.aC_P, C1_min, col = "blue", label = "C1 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, C1_max, col = "blue", label = "C1 max")

plot(df_eq.aC_P, C2_min, col = "green", label = "C2 min", xlabel = "aC_P", ylabel ="min/max")
plot!(df_eq.aC_P, C2_max, col = "green", label = "C2 max")

plot(df_eq.aC_P, P_cv)


##how does aGP influence? - unforced model
results_aG_all_uf = []
for aG_P in 0.0:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1,  aG_P=aG_P)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_aG_all_uf, (; aG_P=aG_P, eq_data...))
end
df_eq = DataFrame(results_aG_all_uf)

plot(df_eq.aG_P, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.aG_P, R1_cv, label = "R1", xlabel = "aG_P", ylabel = "CV")
plot!(df_eq.aG_P, R2_cv, col = "red", label = "R2")
plot!(df_eq.aG_P, C1_cv, col = "blue", label = "C1")
plot!(df_eq.aG_P, C2_cv, col = "green", label = "C2")
plot(df_eq.aG_P, P_cv)
##CV dynamics look very odd, all very low, 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.aG_P, R1_mean, xlabel = "aG_P", ylabel = "Mean")
plot!(df_eq.aG_P, R2_mean)
plot!(df_eq.aG_P, C1_mean)
plot!(df_eq.aG_P, C2_mean)
plot!(df_eq.aG_P, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.aG_P, R1_sd, xlabel = "aG_P", ylabel = "SD")
plot!(df_eq.aG_P, R2_sd)
plot!(df_eq.aG_P, C1_sd)
plot!(df_eq.aG_P, C2_sd)
plot(df_eq.aG_P, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.aG_P, R1_min, label = "R1 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, R1_max, label = "R1 max")

plot(df_eq.aG_P, R2_min, col = "red", label = "R2 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, R2_max, col = "red", label = "R2 max")

plot(df_eq.aG_P, C1_min, col = "blue", label = "C1 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, C1_max, col = "blue", label = "C1 max")

plot(df_eq.aG_P, C2_min, col = "green", label = "C2 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, C2_max, col = "green", label = "C2 max")

plot(df_eq.aG_P, P_cv)

##how does aGP influence - forced model
results_aG_all_f = []
for aG_P in 0.0:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1,  aG_P=aG_P)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_aG_all_f, (; aG_P=aG_P, eq_data...))
end
df_eq = DataFrame(results_aG_all_f)

plot(df_eq.aG_P, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.aG_P, R1_cv, label = "R1", xlabel = "aG_P", ylabel = "CV")
plot!(df_eq.aG_P, R2_cv, col = "red", label = "R2")
plot!(df_eq.aG_P, C1_cv, col = "blue", label = "C1")
plot!(df_eq.aG_P, C2_cv, col = "green", label = "C2")
plot(df_eq.aG_P, P_cv)
##CV dynamics look very odd, all very low, 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.aG_P, R1_mean, xlabel = "aG_P", ylabel = "Mean")
plot!(df_eq.aG_P, R2_mean)
plot!(df_eq.aG_P, C1_mean)
plot!(df_eq.aG_P, C2_mean)
plot!(df_eq.aG_P, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.aG_P, R1_sd, xlabel = "aG_P", ylabel = "SD")
plot!(df_eq.aG_P, R2_sd)
plot!(df_eq.aG_P, C1_sd)
plot!(df_eq.aG_P, C2_sd)
plot(df_eq.aG_P, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.aG_P, R1_min, label = "R1 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, R1_max, label = "R1 max")

plot(df_eq.aG_P, R2_min, col = "red", label = "R2 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, R2_max, col = "red", label = "R2 max")

plot(df_eq.aG_P, C1_min, col = "blue", label = "C1 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, C1_max, col = "blue", label = "C1 max")

plot(df_eq.aG_P, C2_min, col = "green", label = "C2 min", xlabel = "aG_P", ylabel ="min/max")
plot!(df_eq.aG_P, C2_max, col = "green", label = "C2 max")

plot(df_eq.aG_P, P_cv)


##how does e influence? - unforced model
results_e_all_uf = []
for e in 0.0:0.1:1.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, e=e)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_e_all_uf, (; e=e, eq_data...))
end
df_eq = DataFrame(results_e_all_uf)

plot(df_eq.e, df_eq.λ1)

##look at cv, min/max bifurcations
R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.e, R1_cv, label = "R1", xlabel = "e", ylabel = "CV")
plot!(df_eq.e, R2_cv, col = "red", label = "R2")
plot!(df_eq.e, C1_cv, col = "blue", label = "C1")
plot!(df_eq.e, C2_cv, col = "green", label = "C2")
plot(df_eq.e, P_cv)
##CV dynamics look very odd, all very low, 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.e, R1_mean, xlabel = "e", ylabel = "Mean")
plot!(df_eq.e, R2_mean)
plot!(df_eq.e, C1_mean)
plot!(df_eq.e, C2_mean)
plot!(df_eq.e, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.e, R1_sd, xlabel = "e", ylabel = "SD")
plot!(df_eq.e, R2_sd)
plot!(df_eq.e, C1_sd)
plot!(df_eq.e, C2_sd)
plot(df_eq.e, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.e, R1_min, label = "R1 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, R1_max, label = "R1 max")

plot(df_eq.e, R2_min, col = "red", label = "R2 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, R2_max, col = "red", label = "R2 max")

plot(df_eq.e, C1_min, col = "blue", label = "C1 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, C1_max, col = "blue", label = "C1 max")

plot(df_eq.e, C2_min, col = "green", label = "C2 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, C2_max, col = "green", label = "C2 max")

plot(df_eq.e, P_cv)



##how does e influence? - forced model
results_e_all_f = []
for e in 0.0:0.1:1.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, e=e)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_e_all_f, (; e=e, eq_data...))
end
df_eq = DataFrame(results_e_all_f)

plot(df_eq.e, df_eq.λ1)
##eigenvalue flipping likely result of period of time that evaluating across the phase, also this integrated eigenvalue doesnt really make sense anyway so not really worried about it 

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.e, R1_cv, label = "R1", xlabel = "e", ylabel = "CV")
plot!(df_eq.e, R2_cv, col = "red", label = "R2")
plot!(df_eq.e, C1_cv, col = "blue", label = "C1")
plot!(df_eq.e, C2_cv, col = "green", label = "C2")
plot(df_eq.e, P_cv)
##CV dynamics look very odd, all very low, 

R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.e, R1_mean, xlabel = "e", ylabel = "Mean")
plot!(df_eq.e, R2_mean)
plot!(df_eq.e, C1_mean)
plot!(df_eq.e, C2_mean)
plot!(df_eq.e, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.e, R1_sd, xlabel = "e", ylabel = "SD")
plot!(df_eq.e, R2_sd)
plot!(df_eq.e, C1_sd)
plot!(df_eq.e, C2_sd)
plot(df_eq.e, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.e, R1_min, label = "R1 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, R1_max, label = "R1 max")

plot(df_eq.e, R2_min, col = "red", label = "R2 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, R2_max, col = "red", label = "R2 max")

plot(df_eq.e, C1_min, col = "blue", label = "C1 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, C1_max, col = "blue", label = "C1 max")

plot(df_eq.e, C2_min, col = "green", label = "C2 min", xlabel = "e", ylabel ="min/max")
plot!(df_eq.e, C2_max, col = "green", label = "C2 max")

plot(df_eq.e, P_cv)


##how does m influence? - unforced model
results_m_all_uf = []
for mC in 0.0:0.1:1.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, mC=mC)
    eq_data = equilibrium_unforced(p)
    push!(results_m_all_uf, (; mC=mC, eq_data...))
end
df_eq = DataFrame(results_m_all_uf)

plot(df_eq.mC, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.mC, R1_cv, label = "R1", xlabel = "mC", ylabel = "CV")
plot!(df_eq.mC, R2_cv, col = "red", label = "R2")
plot!(df_eq.mC, C1_cv, col = "blue", label = "C1")
plot!(df_eq.mC, C2_cv, col = "green", label = "C2")
plot(df_eq.mC, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.mC, R1_mean, xlabel = "mC", ylabel = "Mean")
plot!(df_eq.mC, R2_mean)
plot!(df_eq.mC, C1_mean)
plot!(df_eq.mC, C2_mean)
plot!(df_eq.mC, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.mC, R1_sd, xlabel = "mC", ylabel = "SD")
plot!(df_eq.mC, R2_sd)
plot!(df_eq.mC, C1_sd)
plot!(df_eq.mC, C2_sd)
plot(df_eq.mC, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.mC, R1_min, label = "R1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, R1_max, label = "R1 max")

plot(df_eq.mC, R2_min, col = "red", label = "R2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, R2_max, col = "red", label = "R2 max")

plot(df_eq.mC, C1_min, col = "blue", label = "C1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, C1_max, col = "blue", label = "C1 max")

plot(df_eq.mC, C2_min, col = "green", label = "C2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, C2_max, col = "green", label = "C2 max")

plot(df_eq.mC, P_cv)


##how does m influence? - forced model
results_m_all_f = []
for mC in 0.0:0.1:1.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, mC=mC)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_m_all_f, (; mC=mC, eq_data...))
end
df_eq = DataFrame(results_m_all_f)

plot(df_eq.mC, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.mC, R1_cv, label = "R1", xlabel = "mC", ylabel = "CV")
plot!(df_eq.mC, R2_cv, col = "red", label = "R2")
plot!(df_eq.mC, C1_cv, col = "blue", label = "C1")
plot!(df_eq.mC, C2_cv, col = "green", label = "C2")
plot(df_eq.mC, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.mC, R1_mean, xlabel = "mC", ylabel = "Mean")
plot!(df_eq.mC, R2_mean)
plot!(df_eq.mC, C1_mean)
plot!(df_eq.mC, C2_mean)
plot!(df_eq.mC, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.mC, R1_sd, xlabel = "mC", ylabel = "SD")
plot!(df_eq.mC, R2_sd)
plot!(df_eq.mC, C1_sd)
plot!(df_eq.mC, C2_sd)
plot(df_eq.mC, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.mC, R1_min, label = "R1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, R1_max, label = "R1 max")

plot(df_eq.mC, R2_min, col = "red", label = "R2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, R2_max, col = "red", label = "R2 max")

plot(df_eq.mC, C1_min, col = "blue", label = "C1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, C1_max, col = "blue", label = "C1 max")

plot(df_eq.mC, C2_min, col = "green", label = "C2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.mC, C2_max, col = "green", label = "C2 max")

plot(df_eq.mC, P_cv)



##how does o influence? - unforced model
results_o_all_uf = []
for o in 0.0:0.1:1.0
    p = ModelPar_active(w = 0.2, o = o, H = 0.1)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_o_all_uf, (; o=o, eq_data...))
end
df_eq = DataFrame(results_o_all_uf)

plot(df_eq.o, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.o, R1_cv, label = "R1", xlabel = "o", ylabel = "CV")
plot!(df_eq.o, R2_cv, col = "red", label = "R2")
plot!(df_eq.o, C1_cv, col = "blue", label = "C1")
plot!(df_eq.o, C2_cv, col = "green", label = "C2")
plot(df_eq.o, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_mean, xlabel = "o", ylabel = "Mean")
plot!(df_eq.o, R2_mean)
plot!(df_eq.o, C1_mean)
plot!(df_eq.o, C2_mean)
plot!(df_eq.o, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_sd, xlabel = "o", ylabel = "SD")
plot!(df_eq.o, R2_sd)
plot!(df_eq.o, C1_sd)
plot!(df_eq.o, C2_sd)
plot(df_eq.o, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.o, R1_min, label = "R1 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, R1_max, label = "R1 max")

plot(df_eq.o, R2_min, col = "red", label = "R2 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, R2_max, col = "red", label = "R2 max")

plot(df_eq.o, C1_min, col = "blue", label = "C1 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, C1_max, col = "blue", label = "C1 max")

plot(df_eq.o, C2_min, col = "green", label = "C2 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, C2_max, col = "green", label = "C2 max")

plot(df_eq.o, P_cv)


##how does o influence? - forced model
results_o_all_f = []
for o in 0.0:0.1:1.0
    p = ModelPar_active(w = 0.2, o = o, H = 0.1)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_o_all_f, (; o=o, eq_data...))
end
df_eq = DataFrame(results_o_all_f)

plot(df_eq.o, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.o, R1_cv, label = "R1", xlabel = "o", ylabel = "CV")
plot!(df_eq.o, R2_cv, col = "red", label = "R2")
plot!(df_eq.o, C1_cv, col = "blue", label = "C1")
plot!(df_eq.o, C2_cv, col = "green", label = "C2")
plot(df_eq.o, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_mean, xlabel = "o", ylabel = "Mean")
plot!(df_eq.o, R2_mean)
plot!(df_eq.o, C1_mean)
plot!(df_eq.o, C2_mean)
plot!(df_eq.o, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_sd, xlabel = "o", ylabel = "SD")
plot!(df_eq.o, R2_sd)
plot!(df_eq.o, C1_sd)
plot!(df_eq.o, C2_sd)
plot(df_eq.o, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.o, R1_min, label = "R1 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, R1_max, label = "R1 max")

plot(df_eq.o, R2_min, col = "red", label = "R2 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, R2_max, col = "red", label = "R2 max")

plot(df_eq.o, C1_min, col = "blue", label = "C1 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, C1_max, col = "blue", label = "C1 max")

plot(df_eq.o, C2_min, col = "green", label = "C2 min", xlabel = "o", ylabel ="min/max")
plot!(df_eq.o, C2_max, col = "green", label = "C2 max")

plot(df_eq.o, P_cv)


##how does w influence? - unforced model
results_w_all_uf = []
for w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = 0.1, H = 0.1)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_w_all_uf, (; w=w, eq_data...))
end
df_eq = DataFrame(results_w_all_uf)

plot(df_eq.w, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.w, R1_cv, label = "R1", xlabel = "w", ylabel = "CV")
plot!(df_eq.w, R2_cv, col = "red", label = "R2")
plot!(df_eq.w, C1_cv, col = "blue", label = "C1")
plot!(df_eq.w, C2_cv, col = "green", label = "C2")
plot(df_eq.w, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_mean, xlabel = "o", ylabel = "Mean")
plot!(df_eq.w, R2_mean)
plot!(df_eq.w, C1_mean)
plot!(df_eq.w, C2_mean)
plot!(df_eq.w, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_sd, xlabel = "o", ylabel = "SD")
plot!(df_eq.w, R2_sd)
plot!(df_eq.w, C1_sd)
plot!(df_eq.w, C2_sd)
plot(df_eq.w, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.w, R1_min, label = "R1 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, R1_max, label = "R1 max")

plot(df_eq.w, R2_min, col = "red", label = "R2 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, R2_max, col = "red", label = "R2 max")

plot(df_eq.w, C1_min, col = "blue", label = "C1 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, C1_max, col = "blue", label = "C1 max")

plot(df_eq.w, C2_min, col = "green", label = "C2 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, C2_max, col = "green", label = "C2 max")

plot(df_eq.w, P_cv)


##how does w influence? - forced model
results_w_all_f = []
for w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = 0.1, H = 0.1)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_w_all_f, (; w=w, eq_data...))
end
df_eq = DataFrame(results_w_all_f)

plot(df_eq.w, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.w, R1_cv, label = "R1", xlabel = "w", ylabel = "CV")
plot!(df_eq.w, R2_cv, col = "red", label = "R2")
plot!(df_eq.w, C1_cv, col = "blue", label = "C1")
plot!(df_eq.w, C2_cv, col = "green", label = "C2")
plot(df_eq.w, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_mean, xlabel = "w", ylabel = "Mean")
plot!(df_eq.w, R2_mean)
plot!(df_eq.w, C1_mean)
plot!(df_eq.w, C2_mean)
plot!(df_eq.w, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_sd, xlabel = "w", ylabel = "SD")
plot!(df_eq.w, R2_sd)
plot!(df_eq.w, C1_sd)
plot!(df_eq.w, C2_sd)
plot(df_eq.w, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.w, R1_min, label = "R1 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, R1_max, label = "R1 max")

plot(df_eq.w, R2_min, col = "red", label = "R2 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, R2_max, col = "red", label = "R2 max")

plot(df_eq.w, C1_min, col = "blue", label = "C1 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, C1_max, col = "blue", label = "C1 max")

plot(df_eq.w, C2_min, col = "green", label = "C2 min", xlabel = "w", ylabel ="min/max")
plot!(df_eq.w, C2_max, col = "green", label = "C2 max")

plot(df_eq.w, P_cv)


##Do over gradient of o and w 
results_o_w_all_uf = []
for o in 0.0:0.1:1.0,  w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = o, H = 0.1)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_o_w_all_uf, (; w=w, o=o, eq_data...))
end
df_eq_o_w = DataFrame(results_o_w_all_uf)

# Get unique values
o_vals = unique(df_eq_o_w.o)
w_vals = unique(df_eq_o_w.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :λ1][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)


##Do over gradient of o and w - forced model
results_o_w_all_f = []
for o in 0.0:0.1:1.0,  w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = o, H = 0.1)
   # t = 100
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_o_w_all_f, (; w=w, o=o, eq_data...))
end
df_eq_o_w = DataFrame(results_o_w_all_f)

# Get unique values
o_vals = unique(df_eq_o_w.o)
w_vals = unique(df_eq_o_w.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :λ1][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

# Create matrix for C1 cv
C1_cv_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :cv][1][3] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, C1_cv_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "CV of C1",
        colorbar_title = "CV",
        c = :viridis)


##how does G influence? - unforced model
results_G_all_uf = []
for G in 0.0:0.1:6.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, G=G)
    P0 = 0.25
    eq_data = equilibrium_unforced(p, P0)
    push!(results_G_all_uf, (; G=G, eq_data...))
end
df_eq = DataFrame(results_G_all_uf)

plot(df_eq.G, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.G, R1_cv, label = "R1", xlabel = "G", ylabel = "CV")
plot!(df_eq.G, R2_cv, col = "red", label = "R2")
plot!(df_eq.G, C1_cv, col = "blue", label = "C1")
plot!(df_eq.G, C2_cv, col = "green", label = "C2")
plot(df_eq.G, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.G, R1_mean, xlabel = "G", ylabel = "Mean")
plot!(df_eq.G, R2_mean)
plot!(df_eq.G, C1_mean)
plot!(df_eq.G, C2_mean)
plot!(df_eq.G, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.G, R1_sd, xlabel = "G", ylabel = "SD")
plot!(df_eq.G, R2_sd)
plot!(df_eq.G, C1_sd)
plot!(df_eq.G, C2_sd)
plot(df_eq.G, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.G, R1_min, label = "R1 min", xlabel = "G", ylabel ="min/max")
plot!(df_eq.G, R1_max, label = "R1 max")

plot(df_eq.G, R2_min, col = "red", label = "R2 min", xlabel = "G", ylabel ="min/max")
plot!(df_eq.G, R2_max, col = "red", label = "R2 max")

plot(df_eq.G, C1_min, col = "blue", label = "C1 min", xlabel = "G", ylabel ="min/max")
plot!(df_eq.G, C1_max, col = "blue", label = "C1 max")

plot(df_eq.G, C2_min, col = "green", label = "C2 min", xlabel = "G", ylabel ="min/max")
plot!(df_eq.G, C2_max, col = "green", label = "C2 max")

plot(df_eq.G, P_cv)

##how does G influence? - forced model
results_G_all_f = []
for G in 0.0:0.1:1.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, G=G)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_G_all_f, (; G=G, eq_data...))
end
df_eq = DataFrame(results_G_all_f)

plot(df_eq.G, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.G, R1_cv, label = "R1", xlabel = "G", ylabel = "CV")
plot!(df_eq.G, R2_cv, col = "red", label = "R2")
plot!(df_eq.G, C1_cv, col = "blue", label = "C1")
plot!(df_eq.G, C2_cv, col = "green", label = "C2")
plot(df_eq.G, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.G, R1_mean, xlabel = "mC", ylabel = "Mean")
plot!(df_eq.G, R2_mean)
plot!(df_eq.G, C1_mean)
plot!(df_eq.G, C2_mean)
plot!(df_eq.G, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.G, R1_sd, xlabel = "mC", ylabel = "SD")
plot!(df_eq.G, R2_sd)
plot!(df_eq.G, C1_sd)
plot!(df_eq.G, C2_sd)
plot(df_eq.G, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.G, R1_min, label = "R1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.G, R1_max, label = "R1 max")

plot(df_eq.G, R2_min, col = "red", label = "R2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.G, R2_max, col = "red", label = "R2 max")

plot(df_eq.G, C1_min, col = "blue", label = "C1 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.G, C1_max, col = "blue", label = "C1 max")

plot(df_eq.G, C2_min, col = "green", label = "C2 min", xlabel = "mC", ylabel ="min/max")
plot!(df_eq.G, C2_max, col = "green", label = "C2 max")

plot(df_eq.G, P_cv)


##how does H influence? - forced model
results_H_all_f = []

for H in 0.0:0.1:1.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = H, K = 4.0)
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_H_all_f, (; H=H, eq_data...))
end
df_eq = DataFrame(results_H_all_f)

plot(df_eq.H, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.H, R1_cv, label = "R1", xlabel = "H", ylabel = "CV")
plot!(df_eq.H, R2_cv, col = "red", label = "R2")
plot!(df_eq.H, C1_cv, col = "blue", label = "C1")
plot!(df_eq.H, C2_cv, col = "green", label = "C2")
plot(df_eq.H, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.H, R1_mean, xlabel = "H", ylabel = "Mean")
plot!(df_eq.H, R2_mean)
plot!(df_eq.H, C1_mean)
plot!(df_eq.H, C2_mean)
plot!(df_eq.H, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.H, R1_sd, xlabel = "H", ylabel = "SD")
plot!(df_eq.H, R2_sd)
plot!(df_eq.H, C1_sd)
plot!(df_eq.H, C2_sd)
plot(df_eq.H, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.H, R1_min, label = "R1 min", xlabel = "H", ylabel ="min/max")
plot!(df_eq.H, R1_max, label = "R1 max")

plot(df_eq.H, R2_min, col = "red", label = "R2 min", xlabel = "H", ylabel ="min/max")
plot!(df_eq.H, R2_max, col = "red", label = "R2 max")

plot(df_eq.H, C1_min, col = "blue", label = "C1 min", xlabel = "H", ylabel ="min/max")
plot!(df_eq.H, C1_max, col = "blue", label = "C1 max")

plot(df_eq.H, C2_min, col = "green", label = "C2 min", xlabel = "H", ylabel ="min/max")
plot!(df_eq.H, C2_max, col = "green", label = "C2 max")

plot(df_eq.H, P_cv)


##how does P influence? - unforced model
results_P_all_uf = []
for P0 in 0.0:0.1:6.0
    p = ModelPar_active(w=0.2, o=0.1, H=0.1)   # no P here; it's a state
    eq_data = equilibrium_unforced(p, P0)
    push!(results_P_all_uf, (; P0, eq_data...))
end
df_eq = DataFrame(results_P_all_uf)

plot(df_eq.P0, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.P0, R1_cv, label = "R1", xlabel = "P", ylabel = "CV")
plot!(df_eq.P0, R2_cv, col = "red", label = "R2")
plot!(df_eq.P0, C1_cv, col = "blue", label = "C1")
plot!(df_eq.P0, C2_cv, col = "green", label = "C2")
plot(df_eq.P0, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.P0, R1_mean, xlabel = "P", ylabel = "Mean")
plot!(df_eq.P0, R2_mean)
plot!(df_eq.P0, C1_mean)
plot!(df_eq.P0, C2_mean)
plot!(df_eq.P0, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.P0, R1_sd, xlabel = "P", ylabel = "SD")
plot!(df_eq.P0, R2_sd)
plot!(df_eq.P0, C1_sd)
plot!(df_eq.P0, C2_sd)
plot(df_eq.P0, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.P0, R1_min, label = "R1 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, R1_max, label = "R1 max")

plot(df_eq.P0, R2_min, col = "red", label = "R2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, R2_max, col = "red", label = "R2 max")

plot(df_eq.P0, C1_min, col = "blue", label = "C1 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, C1_max, col = "blue", label = "C1 max")

plot(df_eq.P0, C2_min, col = "green", label = "C2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, C2_max, col = "green", label = "C2 max")

plot(df_eq.P0, P_cv)

##how does P influence? - forced model
results_P_all_f = []
for P0 in 0.0:0.1:6.0
    p = ModelPar_active(w=0.2, o=0.1, H=0.1)   # no P here; it's a state
    eq_data = equilibrium_forced(p, P0)
    push!(results_P_all_f, (; P0, eq_data...))
end
df_eq = DataFrame(results_P_all_f)

plot(df_eq.P0, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.P0, R1_cv, label = "R1", xlabel = "P", ylabel = "CV")
plot!(df_eq.P0, R2_cv, col = "red", label = "R2")
plot!(df_eq.P0, C1_cv, col = "blue", label = "C1")
plot!(df_eq.P0, C2_cv, col = "green", label = "C2")
plot(df_eq.P0, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.P0, R1_mean, xlabel = "P", ylabel = "Mean")
plot!(df_eq.P0, R2_mean)
plot!(df_eq.P0, C1_mean)
plot!(df_eq.P0, C2_mean)
plot!(df_eq.P0, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.P0, R1_sd, xlabel = "P", ylabel = "SD")
plot!(df_eq.P0, R2_sd)
plot!(df_eq.P0, C1_sd)
plot!(df_eq.P0, C2_sd)
plot(df_eq.P0, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.P0, R1_min, label = "R1 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, R1_max, label = "R1 max")

plot(df_eq.P0, R2_min, col = "red", label = "R2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, R2_max, col = "red", label = "R2 max")

plot(df_eq.P0, C1_min, col = "blue", label = "C1 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, C1_max, col = "blue", label = "C1 max")

plot(df_eq.P0, C2_min, col = "green", label = "C2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.P0, C2_max, col = "green", label = "C2 max")

plot(df_eq.P0, P_cv)


##how does amplitude of K influence? - forced model
results_l1_all_f = []
for l1 in 0.0:0.1:2.0
    p = ModelPar_active(w=0.2, o=0.1, H=0.1, l1 = l1) 
    P0 = 0.25
    eq_data = equilibrium_forced(p, P0)
    push!(results_l1_all_f, (; l1 = l1, eq_data...))
end
df_eq = DataFrame(results_l1_all_f)

plot(df_eq.l1, df_eq.λ1)


R1_cv = [row.cv[1] for row in eachrow(df_eq)]
println(R1_cv)
R1_max = [row.max[1] for row in eachrow(df_eq)]
R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
R2_max = [row.max[2] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C1_max = [row.max[3] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
C2_max = [row.max[4] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]

plot(df_eq.l1, R1_cv, label = "R1", xlabel = "l1", ylabel = "CV")
plot!(df_eq.l1, R2_cv, col = "red", label = "R2")
plot!(df_eq.l1, C1_cv, col = "blue", label = "C1")
plot!(df_eq.l1, C2_cv, col = "green", label = "C2")
plot(df_eq.l1, P_cv)


R1_mean = [row.mean[1] for row in eachrow(df_eq)]
R2_mean = [row.mean[2] for row in eachrow(df_eq)]
C1_mean = [row.mean[3] for row in eachrow(df_eq)]
C2_mean = [row.mean[4] for row in eachrow(df_eq)]
P_mean = [row.mean[5] for row in eachrow(df_eq)]


plot(df_eq.l1, R1_mean, xlabel = "l1", ylabel = "Mean")
plot!(df_eq.l1, R2_mean)
plot!(df_eq.l1, C1_mean)
plot!(df_eq.l1, C2_mean)
plot!(df_eq.l1, P_mean)

##looking at SD
R1_sd = [row.sd[1] for row in eachrow(df_eq)]
R2_sd = [row.sd[2] for row in eachrow(df_eq)]
C1_sd = [row.sd[3] for row in eachrow(df_eq)]
C2_sd = [row.sd[4] for row in eachrow(df_eq)]
P_sd = [row.sd[5] for row in eachrow(df_eq)]


plot(df_eq.l1, R1_sd, xlabel = "l1", ylabel = "SD")
plot!(df_eq.l1, R2_sd)
plot!(df_eq.l1, C1_sd)
plot!(df_eq.l1, C2_sd)
plot(df_eq.l1, P_sd)

##Look at min/max Plots - bifurcations
plot(df_eq.l1, R1_min, label = "R1 min", xlabel = "l1", ylabel ="min/max")
plot!(df_eq.l1, R1_max, label = "R1 max")

plot(df_eq.l1, R2_min, col = "red", label = "R2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.l1, R2_max, col = "red", label = "R2 max")

plot(df_eq.l1, C1_min, col = "blue", label = "C1 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.l1, C1_max, col = "blue", label = "C1 max")

plot(df_eq.l1, C2_min, col = "green", label = "C2 min", xlabel = "P", ylabel ="min/max")
plot!(df_eq.l1, C2_max, col = "green", label = "C2 max")

plot(df_eq.l1, P_cv)


        ##OLD CODE
##how does o influence?
results_o_all = []
for o in 0.1:0.1:1.0
    p = ModelPar_active(w = 0.2, o = o, H = 0.1)
    eq_data = equilibrium_forced(p)
    push!(results_o_all, (; o=o, eq_data...))
end
df_eq = DataFrame(results_o_all)

plot(df_eq.o, df_eq.λ_integrated)

##how does w influence?
results_w_all = []
for w in 0.1:0.1:1.0
    p = ModelPar_active(w = w, o = 0.1, H = 0.1)
    eq_data = equilibrium_forced(p)
    push!(results_w_all, (; w=w, eq_data...))
end
df_eq = DataFrame(results_w_all)

plot(df_eq.w, df_eq.λ_integrated)

##how does H influence?
results_H_all = []
for H in 0.1:0.1:1.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = H)
    eq_data = equilibrium_forced(p)
    push!(results_H_all, (; H=H, eq_data...))
end
df_eq = DataFrame(results_H_all)

plot(df_eq.H, df_eq.λ_integrated)

##how does G influence?
results_G_all = []
for G in 0.1:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, G = G)
    eq_data = equilibrium_forced(p)
    push!(results_G_all, (; G=G, eq_data...))
end
df_eq = DataFrame(results_G_all)

plot(df_eq.G, df_eq.λ_integrated)
 

##how does D (synchrony) influence?
results_D_all = []
for D in 0.0:0.1:0.5
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, D = D)
    eq_data = equilibrium_forced(p)
    push!(results_D_all, (; D=D, eq_data...))
end
df_eq = DataFrame(results_D_all)

plot(df_eq.D, df_eq.λ_integrated)


##Calculate return time 


##Do over gradient of o and w 
results_o_w_all_uf = []
for o in 0.0:0.1:1.0,  w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = o, H = 0.1, K=3, D = 0.0)
    t = 100
    P_fixed = 1.0
    eq_data = equilibrium_unforced(p, t)
    push!(results_o_w_all_uf, (; w=w, o=o, eq_data...))
end
df_eq_o_w = DataFrame(results_o_w_all_uf)

# Get unique values
o_vals = unique(df_eq_o_w.o)
w_vals = unique(df_eq_o_w.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :λ1][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)


##Do over gradient of o and w - with temporal forcomg
results_o_w_all = []
for o in 0.0:0.1:1.0,  w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = o, H = 0.9, K=3, D = 0.5, pf = 1.0)
    t = 100
    P_fixed = 0.25
    eq_data = equilibrium_forced(p)
    push!(results_o_w_all, (; w=w, o=o, eq_data...))
end
df_eq_o_w = DataFrame(results_o_w_all)

# Get unique values
o_vals = unique(df_eq_o_w.o)
w_vals = unique(df_eq_o_w.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :λ_integrated][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

##cv matrix 
C1_cv_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :cv][1][3] for w in w_vals, o in o_vals]

C2_cv_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :cv][1][4] for w in w_vals, o in o_vals]


R1_cv_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :cv][1][1] for w in w_vals, o in o_vals]
R2_cv_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :cv][1][2] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, C2_cv_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "C2-CV",
        colorbar_title = "cv",
        c = :viridis)

##min value matrix 
C1_min_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :min][1][3] for w in w_vals, o in o_vals]

C2_min_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :min][1][4] for w in w_vals, o in o_vals]


R1_min_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :min][1][1] for w in w_vals, o in o_vals]
R2_min_mat = [df_eq_o_w[(df_eq_o_w.o .== o) .& (df_eq_o_w.w .== w), :min][1][2] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, C2_min_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "C2-min",
        colorbar_title = "min",
        c = :viridis)





plot(df_eq_o_w.t_eig,df_eq_o_w.λ_max_vals, label="λₘₐₓ(t)", xlabel="Time", ylabel="Re(λₘₐₓ)")








##how does scale of r to pf influence things?
results_r_pf_all = []
for r in 0.1:0.1:5.0,  pf in 0.1:1.0:10.0
    p = ModelPar_active(r = r, pf = pf, H = 0.3, K=3, D = 0.0, o=0.7, w = 0.5)
    eq_data = equilibrium_forced(p)
    push!(results_r_pf_all, (; r=r, pf=pf, eq_data...))
end
df_eq_r_pf = DataFrame(results_r_pf_all)

# Get unique values
r_vals = unique(df_eq_r_pf.r)
pf_vals = unique(df_eq_r_pf.pf)

# Sort them to be safe
sort!(r_vals)
sort!(pf_vals)

# Create matrix for λmax
λ_mat = [df_eq_r_pf[(df_eq_r_pf.r .== r) .& (df_eq_r_pf.pf .== pf), :λ_integrated][1] for r in r_vals, pf in pf_vals]

# Plot heatmap
heatmap(r_vals, pf_vals, λ_mat;
        xlabel = "growth rate (r)",
        ylabel = "phase length (pf)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)


 ##do over gradient of synchrony/asynchrony 
p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0)
P_fixed = 0.25
D_results = equilibrium_forced(p)

results_D_all = []
for D in 0.0:0.1:1.0
    p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0, D=D)
    eq_data = equilibrium_forced(p)
    push!(results_D_all, (; D=D, eq_data...))
end
df_eq = DataFrame(results_D_all)

plot(df_eq.D, df_eq.λ_integrated)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]


plot(df_eq.D, R1_cv, label = "R1")
plot!(df_eq.D, R2_cv, col = "red", label = "R2")
plot!(df_eq.D, C1_cv, col = "green", label = "C1")
plot!(df_eq.D, C2_cv, col = "purple", label = "C2")
plot(df_eq.D, P_cv)

R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_min = [row.min[5] for row in eachrow(df_eq)]


plot(df_eq.D, R1_min, label = "R1")
plot!(df_eq.D, R2_min, col = "red", label = "R2")
plot!(df_eq.D, C1_min, col = "green", label = "C1")
plot!(df_eq.D, C2_min, col = "purple", label = "C2")
plot(df_eq.D, P_min)


##what about just for w, and what parameter can i use to reduce any temporal forcing? - just make l1 and l2 0 
P_fixed = 0.25

results_w_all = []
for w in 0.0:0.1:1.0
    p = ModelPar_passive(w = w, o = 0.0, H = 0.0, l1 = 0.1, l2 = 0.1)
    eq_data = equilibrium_forced(p)
    push!(results_w_all, (; w=w, eq_data...))
end
df_eq = DataFrame(results_w_all)

plot(df_eq.w, df_eq.λ_integrated)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_cv, label = "R1")
plot!(df_eq.w, R2_cv, col = "red", label = "R2")
plot!(df_eq.w, C1_cv, col = "green", label = "C1")
plot!(df_eq.w, C2_cv, col = "purple", label = "C2")
plot(df_eq.w, P_cv)

R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_min = [row.min[5] for row in eachrow(df_eq)]


plot(df_eq.w, R1_min, label = "R1")
plot!(df_eq.w, R2_min, col = "red", label = "R2")
plot!(df_eq.w, C1_min, col = "green", label = "C1")
plot!(df_eq.w, C2_min, col = "purple", label = "C2")
plot(df_eq.w, P_min)


##what about just for o, and what parameter can i use to reduce any temporal forcing? - just make l1 and l2 0 
P_fixed = 0.25

results_o_all = []
for o in 0.0:0.1:1.0
    p = ModelPar_passive(w = 0.0, o = o, H = 0.0, l1 = 0.1, l2 = 0.1)
    eq_data = equilibrium_forced(p)
    push!(results_o_all, (; o=o, eq_data...))
end
df_eq = DataFrame(results_o_all)

plot(df_eq.o, df_eq.λ_integrated)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_cv, label = "R1")
plot!(df_eq.o, R2_cv, col = "red", label = "R2")
plot!(df_eq.o, C1_cv, col = "green", label = "C1")
plot!(df_eq.o, C2_cv, col = "purple", label = "C2")
plot(df_eq.o, P_cv)

R1_min = [row.min[1] for row in eachrow(df_eq)]
R2_min = [row.min[2] for row in eachrow(df_eq)]
C1_min = [row.min[3] for row in eachrow(df_eq)]
C2_min = [row.min[4] for row in eachrow(df_eq)]
P_min = [row.min[5] for row in eachrow(df_eq)]


plot(df_eq.o, R1_min, label = "R1")
plot!(df_eq.o, R2_min, col = "red", label = "R2")
plot!(df_eq.o, C1_min, col = "green", label = "C1")
plot!(df_eq.o, C2_min, col = "purple", label = "C2")
plot(df_eq.o, P_min)