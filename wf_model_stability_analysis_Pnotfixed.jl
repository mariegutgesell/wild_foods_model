##Looking at stability across gradient of coupling and omnivory 

##source model 
#include("wf_model_eqs_subsidy.jl")

include("wf_model_eqs_subsidy_2.jl") ##model equations with P constant/or not (depending on which equation on model structure is silenced), unique parameters per trophic level, active and passive omnivory parameter structures
##WHERE LEFT OFF (JULY 3): trying to understand if dynamics from simpler to more complex model match what i would expect based on theory - working through this


##1) Looking at local stability for unforced model to see changes w/ increasing K
p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0)

K_results = equilibrium_forced(p)

results_K_all = []
for K in 0.1:0.1:10.0
    p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0, K=K)
    eq_data = equilibrium_forced(p)
    push!(results_K_all, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all)

plot(df_eq.K, df_eq.λ_integrated)

R1_cv = [row.cv[1] for row in eachrow(df_eq)]
R2_cv = [row.cv[2] for row in eachrow(df_eq)]
C1_cv = [row.cv[3] for row in eachrow(df_eq)]
C2_cv = [row.cv[4] for row in eachrow(df_eq)]
P_cv = [row.cv[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_cv)
plot(df_eq.K, R2_cv)
plot(df_eq.K, C1_cv)
plot(df_eq.K, C2_cv)
plot(df_eq.K, P_cv)


R1_mean = [row.mean_state[1] for row in eachrow(df_eq)]
R2_mean = [row.mean_state[2] for row in eachrow(df_eq)]
C1_mean = [row.mean_state[3] for row in eachrow(df_eq)]
C2_mean = [row.mean_state[4] for row in eachrow(df_eq)]
P_mean = [row.mean_state[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_mean)
plot(df_eq.K, R2_mean)
plot(df_eq.K, C1_mean)
plot(df_eq.K, C2_mean)
plot(df_eq.K, P_mean)


##come back to these observations below, and if still relevant, also develop bifurcation diagrams .. 
##feeling slightly more confident about approach of calculating stability? 
##time to dig into the dynamics ... 


##the coupling plot does not look like how i would expect, why max stability at relatively strong coupling? is this related to subsidy? 
##why does omnivory always destabilize at low P values? 

##with the parameters the same across both energy channels, with equal coupling, the P basically just has more access to resources, all acting like one pool 
##essentially with this set up, all one pool ... just increasing productivity? 

##Calculate return time 


