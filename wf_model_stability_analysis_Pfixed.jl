##Looking at stability across gradient of coupling and omnivory 

##source model 
#include("wf_model_eqs_subsidy.jl")

include("wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P constant/or not (depending on which equation on model structure is silenced), unique parameters per trophic level, active and passive omnivory parameter structures
##WHERE LEFT OFF (JULY 3): trying to understand if dynamics from simpler to more complex model match what i would expect based on theory - working through this


##1) Looking at local stability for unforced model to see changes w/ increasing K
p = ModelPar_passive(w = 0.5, o = 0.0, H = 0.0)
P_fixed = 0.25
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


##Do over gradient of o and w 

results_o_w_all = []
for o in 0.0:0.1:1.0,  w in 0.0:0.1:1.0
    p = ModelPar_active(w = w, o = o, H = 0.3, K=3, D = 0.0)
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
##values look very negative, something odd going on in eigenvalue integration, debug
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