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
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, K=K)
    t = 100S
    eq_data = equilibrium_unforced(p, t)
    push!(results_K_all, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all)

plot(df_eq.K, df_eq.λ_integrated)

R1_cv = [row.min[1] for row in eachrow(df_eq)]
R1_max = [row.max[1] for row in eachrow(df_eq)]
R2_cv = [row.min[2] for row in eachrow(df_eq)]
C1_cv = [row.min[3] for row in eachrow(df_eq)]
C2_cv = [row.min[4] for row in eachrow(df_eq)]
P_cv = [row.min[5] for row in eachrow(df_eq)]


plot(df_eq.K, R1_cv, label = "R1")
plot!(df_eq.K, R2_cv, col = "red", label = "R2")
plot!(df_eq.K, C1_cv, col = "blue", label = "C1")
plot!(df_eq.K, C2_cv, col = "green", label = "C2")
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

##unforced model
results_K_all_uf = []
for K in 0.1:0.1:5.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, K=K)
    t = 100
    eq_data = equilibrium_unforced(p, t)
    push!(results_K_all_uf, (; K=K, eq_data...))
end
df_eq = DataFrame(results_K_all_uf)

plot(df_eq.K, df_eq.λ1)

##how does a influence?
results_a_all = []
for aR_P in 0.0:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, aR_P=aR_P)
    eq_data = equilibrium_forced(p)
    push!(results_a_all, (; aR_P=aR_P, eq_data...))
end
df_eq = DataFrame(results_a_all)

plot(df_eq.aR_P, df_eq.λ1)

##how does e influence?
results_e_all = []
for e in 0.1:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, e=e)
    eq_data = equilibrium_forced(p)
    push!(results_e_all, (; e=e, eq_data...))
end
df_eq = DataFrame(results_e_all)

plot(df_eq.e, df_eq.λ_integrated)

##how does r influence?
results_r_all = []
for r in 0.1:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, r=r)
    eq_data = equilibrium_forced(p)
    push!(results_r_all, (; r=r, eq_data...))
end
df_eq = DataFrame(results_r_all)

plot(df_eq.r, df_eq.λ_integrated)

##how does h influence?
results_h_all = []
for hR_P in 0.0:0.1:10.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, hR_P=hR_P)
    eq_data = equilibrium_forced(p)
    push!(results_h_all, (; hR_P=hR_P, eq_data...))
end
df_eq = DataFrame(results_h_all)

plot(df_eq.hR_P, df_eq.λ_integrated)

##how does m influence?
results_m_all = []
for mC in 0.1:0.1:1.0
    p = ModelPar_active(w = 0.2, o = 0.1, H = 0.1, mC=mC)
    eq_data = equilibrium_forced(p)
    push!(results_m_all, (; mC=mC, eq_data...))
end
df_eq = DataFrame(results_m_all)
##high mortality rates cause errors/NAs in jacobian matrix b/c dividing by 0, so no stable solution 
plot(df_eq.mC, df_eq.λ_integrated)

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