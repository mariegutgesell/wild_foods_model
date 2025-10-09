##Looking at CV and min of harvest rates (P consumption)

##source model 
#include("wf_model_eqs_subsidy.jl")

include("wf_model_eqs_subsidy_Pfixed.jl")

# ---- CV helper ----
cv(x) = mean(x) == 0 ? NaN : (std(x) / mean(x))

# ---- run, sample post-warmup, compute FR_into_P series + CVs ----
function fr_cv_unforced(p; u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    tspan = (0.0, t_eval)
    prob  = ODEProblem(model_unforced!, u0, tspan, deepcopy(p))
    sol   = solve(prob; reltol=1e-8, abstol=1e-8)

    t_grid  = range(t_warmup, t_eval; length=ngrid)
    us      = sol.(t_grid)

    # time series of fluxes into P
    fr_total = Float64[]; fr_R1 = Float64[]; fr_R2 = Float64[]
    fr_C1 = Float64[];    fr_C2 = Float64[]; fr_G  = Float64[]

    for (u, t) in zip(us, t_grid)
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t)
        push!(fr_total, tot); push!(fr_R1, r1); push!(fr_R2, r2)
        push!(fr_C1, c1);     push!(fr_C2, c2); push!(fr_G,  g)
    end

    return (; 
        cv_total = cv(fr_total),
        cv_R1    = cv(fr_R1),
        cv_R2    = cv(fr_R2),
        cv_C1    = cv(fr_C1),
        cv_C2    = cv(fr_C2),
        cv_G     = cv(fr_G)
    )
end

function fr_cv_forced(p; u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    tspan = (0.0, t_eval)
    prob  = ODEProblem(model_forced!, u0, tspan, deepcopy(p))
    sol   = solve(prob; reltol=1e-8, abstol=1e-8)

    t_grid  = range(t_warmup, t_eval; length=ngrid)
    us      = sol.(t_grid)

    # time series of fluxes into P
    fr_total = Float64[]; fr_R1 = Float64[]; fr_R2 = Float64[]
    fr_C1 = Float64[];    fr_C2 = Float64[]; fr_G  = Float64[]

    for (u, t) in zip(us, t_grid)
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t)
        push!(fr_total, tot); push!(fr_R1, r1); push!(fr_R2, r2)
        push!(fr_C1, c1);     push!(fr_C2, c2); push!(fr_G,  g)
    end

    return (; 
        cv_total = cv(fr_total),
        cv_R1    = cv(fr_R1),
        cv_R2    = cv(fr_R2),
        cv_C1    = cv(fr_C1),
        cv_C2    = cv(fr_C2),
        cv_G     = cv(fr_G),
        mean_total = mean(fr_total),
        mean_R1 = mean(fr_R1),
        mean_R2 = mean(fr_R2),
        mean_C1 = mean(fr_C1),
        mean_C2 = mean(fr_C2),
        mean_G = mean(fr_G),
        sd_total = std(fr_total),
        sd_R1 = std(fr_R1),
        sd_R2 = std(fr_R2),
        sd_C1 = std(fr_C1),
        sd_C2 = std(fr_C2),
        sd_G = std(fr_G),
    )
end


##Look at across range of parameters
results_K_uf = []
for K in 0.1:0.1:6.5
    pᵢ = ModelPar_active(w=0.2, o=0.1, H=0.1, K=K)  # add other defaults as needed
    cv_nt = fr_cv_unforced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_K_uf, (; K, cv_nt...))  # NamedTuple splat into the row
end
df_cv_K = DataFrame(results_K_uf)

plot(df_cv_K.K, df_cv_K.cv_total)


results_K_f = []
for K in 0.1:0.1:6.5
    pᵢ = ModelPar_active(w=0.2, o=0.1, H=0.1, K=K)  # add other defaults as needed
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_K_f, (; K, cv_nt...))  # NamedTuple splat into the row
end
df_cv_K = DataFrame(results_K_f)

plot(df_cv_K.K, df_cv_K.cv_total, label = "total consumption", xlabel = "K", ylabel = "CV")
plot!(df_cv_K.K, df_cv_K.cv_R1, label = "R1 consumption", xlabel = "K", ylabel = "CV")
plot!(df_cv_K.K, df_cv_K.cv_R2, label = "R2 consumption", xlabel = "K", ylabel = "CV")
plot!(df_cv_K.K, df_cv_K.cv_C1, label = "C1 consumption", xlabel = "K", ylabel = "CV")
plot!(df_cv_K.K, df_cv_K.cv_C2, label = "C2 consumption", xlabel = "K", ylabel = "CV")
plot(df_cv_K.K, df_cv_K.cv_G, label = "G consumption", xlabel = "K", ylabel = "CV")


results_o_f = []
for o in 0.0:0.1:1.0
    u0 = [5.0, 4.0, 3.0, 2.0, 0.25]
    pᵢ = ModelPar_active(w=0.5, o=o, H=0.5, K = 3.0)  # add other defaults as needed
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_o_f, (; o, cv_nt...))  # NamedTuple splat into the row
end
df_cv_o = DataFrame(results_o_f)
print(df_cv_o)
plot(df_cv_o.o, df_cv_o.cv_total, label = "total consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_R1, label = "R1 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_R2, label = "R2 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_C1, label = "C1 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_C2, label = "C2 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_G, label = "G consumption", xlabel = "o", ylabel = "CV")


##trying to plot total harvest cv over values of o 
# use categorical tick labels so bars are discrete by o
o_labels = string.(round.(df_cv_o.o, digits=1))

bar(
    o_labels, df_cv_o.cv_total;
    xlabel = "o (omnivory preference)",
    ylabel = "CV of total harvest",
    title  = "Total harvest CV across o",
    legend = false,
    bar_width = 0.8
)

bar(
    o_labels, df_cv_o.mean_total;
    xlabel = "o (omnivory preference)",
    ylabel = "Mean total harvest",
    title  = "Mean Total harvest across o",
    legend = false,
    bar_width = 0.8
)

bar(
    o_labels, df_cv_o.sd_total;
    xlabel = "o (omnivory preference)",
    ylabel = "SD total harvest",
    title  = "SD Total harvest across o",
    legend = false,
    bar_width = 0.8
)






##do range of w next
results_w_f = []
for w in 0.0:0.1:1.0
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]
    pᵢ = ModelPar_active(w=w, o=0.5, H=0.5, K = 3.0)  # add other defaults as needed
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_w_f, (; w, cv_nt...))  # NamedTuple splat into the row
end
df_cv_w = DataFrame(results_w_f)

plot(df_cv_w.w, df_cv_w.cv_total, label = "total consumption", xlabel = "w", ylabel = "CV")
plot!(df_cv_w.w, df_cv_w.cv_R1, label = "R1 consumption", xlabel = "w", ylabel = "CV")
plot!(df_cv_w.w, df_cv_w.cv_R2, label = "R2 consumption", xlabel = "w", ylabel = "CV")
plot!(df_cv_w.w, df_cv_w.cv_C1, label = "C1 consumption", xlabel = "w", ylabel = "CV")
plot!(df_cv_w.w, df_cv_w.cv_C2, label = "C2 consumption", xlabel = "w", ylabel = "CV")
plot!(df_cv_w.w, df_cv_w.cv_G, label = "G consumption", xlabel = "w", ylabel = "CV")

w_labels = string.(round.(df_cv_w.w, digits=1))

bar(
    w_labels, df_cv_w.cv_total;
    xlabel = "w (habitat preference)",
    ylabel = "CV of total harvest",
    title  = "Total harvest CV across w",
    legend = false,
    bar_width = 0.8
)

bar(
    w_labels, df_cv_w.mean_total;
    xlabel = "w (habitat preference)",
    ylabel = "Mean total harvest",
    title  = "Mean Total harvest across w",
    legend = false,
    bar_width = 0.8
)

bar(
    w_labels, df_cv_w.sd_total;
    xlabel = "w (habitat preference)",
    ylabel = "SD total harvest",
    title  = "SD Total harvest across w",
    legend = false,
    bar_width = 0.8
)


##do range of H next
results_H_f = []
for H in 0.0:0.1:1.0
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]
    pᵢ = ModelPar_active(w=0.5, o=0.5, H=H, K = 3.0)  # add other defaults as needed
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_H_f, (; H, cv_nt...))  # NamedTuple splat into the row
end
df_cv_w = DataFrame(results_H_f)

plot(df_cv_w.H, df_cv_w.cv_total, label = "total consumption", xlabel = "H", ylabel = "CV")
plot!(df_cv_w.H, df_cv_w.cv_R1, label = "R1 consumption", xlabel = "H", ylabel = "CV")
plot!(df_cv_w.H, df_cv_w.cv_R2, label = "R2 consumption", xlabel = "H", ylabel = "CV")
plot!(df_cv_w.H, df_cv_w.cv_C1, label = "C1 consumption", xlabel = "H", ylabel = "CV")
plot!(df_cv_w.H, df_cv_w.cv_C2, label = "C2 consumption", xlabel = "H", ylabel = "CV")
plot!(df_cv_w.H, df_cv_w.cv_G, label = "G consumption", xlabel = "H", ylabel = "CV")

H_labels = string.(round.(df_cv_w.H, digits=1))

bar(
    H_labels, df_cv_w.cv_total;
    xlabel = "H (grocery preference)",
    ylabel = "CV of total harvest",
    title  = "Total harvest CV across H",
    legend = false,
    bar_width = 0.8
)

bar(
    H_labels, df_cv_w.mean_total;
    xlabel = "H (grocery preference)",
    ylabel = "Mean total harvest",
    title  = "Mean Total harvest across H",
    legend = false,
    bar_width = 0.8
)

bar(
    H_labels, df_cv_w.sd_total;
    xlabel = "H (grocery preference)",
    ylabel = "SD total harvest",
    title  = "SD Total harvest across H",
    legend = false,
    bar_width = 0.8
)



##trying heatmap of CV of total harvest

o_vals = 0.0:0.1:1.0
w_vals = 0.0:0.1:1.0

results_grid = []
for o in o_vals, w in w_vals
    pᵢ = ModelPar_active(w=w, o=o, H=0.9, l = 1.0, D = 0.5)  # set/adjust other params as you need
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=600)
    push!(results_grid, (; o, w, cv_nt...))
end
df_cv_o_w = DataFrame(results_grid)


o_vals = unique(df_cv_o_w.o)
w_vals = unique(df_cv_o_w.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for cv total
cv_total_mat = [df_cv_o_w[(df_cv_o_w.o .== o) .& (df_cv_o_w.w .== w), :cv_total][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, cv_total_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "CV of Total Harvest",
        colorbar_title = "CV",
        c = :viridis)