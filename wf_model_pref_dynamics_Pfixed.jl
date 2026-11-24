##Looking at dynamics of density dependent preference parameters 
using Pkg
#Pkg.add("CairoMakie")
#Pkg.add("Colors")
using Colors
using CairoMakie
##source model 
#include("wf_model_eqs_subsidy.jl")

include("wf_model_eqs_subsidy_Pfixed.jl")

u0 = [1.5, 1.5, 1.0, 1.0, 0.5]

##Look at how mean o across range of parameters
results_o_uf = []
for o in 0:0.1:1.0
    p = ModelPar_active(w=0.5, o = o, H=0.0)  # add other defaults as needed
    mean_o = pref_dynamics_unforced(p; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_o_uf, (; o, mean_o...))  # NamedTuple splat into the row
end
df_o = DataFrame(results_o_uf)

plot(df_o.time, df_o.W_ts, label="w")
plot(df_o.time, df_o.om_i_ts, label="o_i")
plot(df_o.time, df_o.om_j_ts, label="o_j")
plot(df_o.time, df_o.H_ts, label="H")

xlabel!("Time")
ylabel!("Realized preference")


cols = palette(:viridis, nrow(df_o))

p1 = plot(title="Habitat preference (w)")
p2 = plot(title="Omnivory (o_i)")
p3 = plot(title="Omnivory (o_j)")
p4 = plot(title="Grocery preference (H)")

for (i, row) in enumerate(eachrow(df_o))

    lab = "o = $(round(row.o, digits=1))"

    plot!(p1, row.time, row.W_ts, 
          color=cols[i], label=lab)

    plot!(p2, row.time, row.om_i_ts, 
          color=cols[i], label=lab)

    plot!(p3, row.time, row.om_j_ts, 
          color=cols[i], label=lab)

    plot!(p4, row.time, row.H_ts, 
          color=cols[i], label=lab)
end

plot(p1, p2, p3, p4,
     layout=(2,2),
     size=(1000,700),
     xlabel="Time",
     ylabel="Realized preference",
     ylims=(0,1),
     linewidth=2)


##Look at how mean o across range of parameters
results_o_f = []
for o in 0:0.1:1.0
    p = ModelPar_active(w=0.5, o = o, H=0.0)  # add other defaults as needed
    mean_o = pref_dynamics_forced(p; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)
    push!(results_o_f, (; o, mean_o...))  # NamedTuple splat into the row
end
df_o = DataFrame(results_o_f)

plot(df_o.time, df_o.W_ts, label="w")
plot(df_o.time, df_o.om_i_ts, label="o_i")
plot(df_o.time, df_o.om_j_ts, label="o_j")
plot(df_o.time, df_o.H_ts, label="H")

xlabel!("Time")
ylabel!("Realized preference")


cols = palette(:viridis, nrow(df_o))

p1 = plot(title="Habitat preference (w)")
p2 = plot(title="Omnivory (o_i)")
p3 = plot(title="Omnivory (o_j)")
p4 = plot(title="Grocery preference (H)")

for (i, row) in enumerate(eachrow(df_o))

    lab = "o = $(round(row.o, digits=1))"

    plot!(p1, row.time, row.W_ts, 
          color=cols[i], label=lab)
    hline!(p1, [row.pref_W], 
           color=cols[i], linestyle=:dot, label=false)

    plot!(p2, row.time, row.om_i_ts, 
          color=cols[i], label=lab)
    hline!(p2, [row.pref_om_i], 
           color=cols[i], linestyle=:dot, label=false)

    plot!(p3, row.time, row.om_j_ts, 
          color=cols[i], label=lab)
     hline!(p3, [row.pref_om_j], 
           color=cols[i], linestyle=:dot, label=false)

    plot!(p4, row.time, row.H_ts, 
          color=cols[i], label=lab)
    hline!(p4, [row.pref_H], 
           color=cols[i], linestyle=:dot, label=false)
end

plot(p1, p2, p3, p4,
     layout=(2,2),
     size=(1000,700),
     xlabel="Time",
     ylabel="Realized preference",
     ylims=(0,1),
     linewidth=2)
