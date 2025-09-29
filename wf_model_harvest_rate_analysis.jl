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
        cv_G     = cv(fr_G)
    )
end


##Look at across range of parameters
results_K_uf = []
for K in 0.0:0.1:6.5
    pᵢ = ModelPar_active(w=0.2, o=0.1, H=0.1, K=K)  # add other defaults as needed
    cv_nt = fr_cv_unforced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_K_uf, (; K, cv_nt...))  # NamedTuple splat into the row
end
df_cv_K = DataFrame(results_K_uf)

plot(df_cv_K.K, df_cv_K.cv_total)


results_K_f = []
for K in 2.0:0.1:6.5
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
    pᵢ = ModelPar_active(w=0.2, o=o, H=0.1, l1 = 0.5, l2 = 0.5)  # add other defaults as needed
    cv_nt = fr_cv_forced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_o_f, (; o, cv_nt...))  # NamedTuple splat into the row
end
df_cv_o = DataFrame(results_o_f)

plot(df_cv_o.o, df_cv_o.cv_total, label = "total consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_R1, label = "R1 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_R2, label = "R2 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_C1, label = "C1 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_C2, label = "C2 consumption", xlabel = "o", ylabel = "CV")
plot!(df_cv_o.o, df_cv_o.cv_G, label = "G consumption", xlabel = "o", ylabel = "CV")



##do range of w next
results_w_f = []
for w in 0.0:0.1:1.0
    pᵢ = ModelPar_active(w=w, o=0.1, H=0.1, l1 = 0.5, l2 = 0.5)  # add other defaults as needed
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






first(df_cv_K, 5)
labels  = ["Total","R1 → P","R2 → P","C1 → P","C2 → P","G → P"]
row     = last(eachrow(df_cv_K))
cv_vals = [row.cv_total, row.cv_R1, row.cv_R2, row.cv_C1, row.cv_C2, row.cv_G]

bar(labels, cv_vals; legend=false, xlabel="Flux into Predator",
    ylabel="Coefficient of Variation (CV)", title="Predator Consumption CV")


    ##trying heatmap of CV of total harvest-- NOT WORKING 

o_vals = 0.0:0.1:1.0
w_vals = 0.0:0.1:1.0

results_grid = []
for o in o_vals, w in w_vals
    pᵢ = ModelPar_active(w=w, o=o, H=0.1, aC_P=1.0, K = 5.0)  # set/adjust other params as you need
    cv_nt = fr_cv_unforced(pᵢ; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=600)
    push!(results_grid, (; o, w, cv_nt...))
end
df_cv_ow = DataFrame(results_grid)


# make a matrix Z with rows = unique(o), cols = unique(w)
sort!(df_cv_ow, [:o, :w])
O = unique(df_cv_ow.o); W = unique(df_cv_ow.w)
Z = Array{Float64}(undef, length(O), length(W))
for (i, oi) in enumerate(O), (j, wj) in enumerate(W)
    Z[i, j] = df_cv_ow[(df_cv_ow.o .== oi) .& (df_cv_ow.w .== wj), :cv_total][1]
end

heatmap(W, O, Z; xlabel="w", ylabel="o", colorbar_title="CV",
        title="CV(total → P) across (o, w)")