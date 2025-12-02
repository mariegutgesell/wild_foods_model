##Function to run pulse experiment -- modifying approach from KC




include("wf_model_eqs_subsidy_Pfixed.jl")


# Entire analysis for one pulse simulation (1 system)
#So: pulse is 200–210;  integrate 0–350;  analyze 175–350
##might need to adjust these to have a longer burn in period 
function pulse_unit_forced(p; p_length = 50.0, p_strength = 0.25, pre_window = 20.0, post_window = 5.0) 
    # Initial conditions
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]            #initial state of R1, R2, C1, C1 and Pfixed -- define P_fixed outside of function
    
    ##Timing
    pulse_start = 200.0               ##time when pulse begins
    pulse_end = pulse_start + p_length      #time when pulse ends
   # pulse_event_times = union(pulse_start, pulse_end) #times where the pulse changes system
    t_end = 350.0         #end of simulation
    t_span = (0.0, t_end)       #time span of simulation
    t_start = 175.0         #start of post-pulse tracking (i.e., when tracking starts to capture pulse and dynamics after)
    len = 1000           ##resolution of the time grid
     t_grid = range(t_start, t_end, length = len)        #time grid for evaluating solution
    t_after = range(pulse_end, t_end, length = len)     #time grid after the pulse
   
   # function for G in parameter set is: G_func(t) = (t < t_pulse || t ≥ t_recover) ? G_pre : G_pulse

    #I think define par outside of function ... 

    # Establish baseline G (try :G_base then :G)
    G_pre = if hasproperty(p, :G_base)
        p.G_base
    elseif hasproperty(p, :G)
        p.G
    else
        error("Param set needs a baseline G (e.g., `G_base` or `G`).")
    end
    p.G_base = G_pre  # keep a record

    G_pulse = p_strength * G_pre

    p.G_func = t -> (t < pulse_start || t ≥ pulse_end) ? G_pre : G_pulse

    ##Pulse Event Logic ##
   
    ##Simulate Pulse for a Given parameter set
   # par_pulse = deepcopy(par)
   # par_pulse.G = p_strength * par.G_base  ##initial pulse value (will be reset via callback)
    prob = ODEProblem(model_forced!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    ##evaluate on grid 
    u_grid = sol(t_grid)
    
   
    # ----- total harvest time series -----
    fr_total = Vector{Float64}(undef, length(t_grid))
    fr_R1 = similar(fr_total); fr_R2 = similar(fr_total)
    fr_C1 = similar(fr_total); fr_C2 = similar(fr_total); fr_G = similar(fr_total)

    for (i, (u, t)) in enumerate(zip(u_grid, t_grid))
        tot, r1, r2, c1, c2, g = total_FR_into_P(u, p, t)
        fr_total[i] = tot
        fr_R1[i] = r1; fr_R2[i] = r2
        fr_C1[i] = c1; fr_C2[i] = c2; fr_G[i] = g
    end
 
    # ----- decline metrics -----
    pre_mask  = (t_grid .>= (pulse_start - pre_window)) .& (t_grid .< pulse_start)
    post_mask = (t_grid .>=  pulse_start)               .& (t_grid .< (pulse_start + post_window))

    @assert any(pre_mask)  "No samples in pre-pulse window; increase `pre_window` or resolution."
    @assert any(post_mask) "No samples in post-pulse window; increase `post_window` or resolution."

    baseline = mean(fr_total[pre_mask])

    post_min, idx_rel = findmin(fr_total[post_mask])
    first_post        = findfirst(post_mask)
    idx_post_min      = first_post - 1 + idx_rel
    t_post_min        = t_grid[idx_post_min]

    decline_abs            = baseline - post_min
    decline_rel_baseline   = decline_abs / baseline
    shock_size             = abs(p_strength - 1.0)                    # dimensionless shock on G
    decline_per_unit_shock = decline_abs / shock_size                 # harvest per unit shock
    decline_dimless        = decline_abs / (baseline * shock_size)    # fully dimensionless sensitivity


    # also compute series on the solver's own (nonuniform) time grid
t_raw = Array(sol.t)
fr_vals_raw = [total_FR_into_P(u, p, t) for (u, t) in zip(sol.u, sol.t)]
fr_total_raw = [x[1] for x in fr_vals_raw]
fr_R1_raw    = [x[2] for x in fr_vals_raw]
fr_R2_raw    = [x[3] for x in fr_vals_raw]
fr_C1_raw    = [x[4] for x in fr_vals_raw]
fr_C2_raw    = [x[5] for x in fr_vals_raw]
fr_G_raw     = [x[6] for x in fr_vals_raw]

  return (;
    baseline,
    post_min,
    decline_abs,
    decline_rel_baseline,
    shock_size,
    decline_per_unit_shock,
    decline_dimless,
    t_grid,
    fr_total,
    fr_components = (; fr_R1, fr_R2, fr_C1, fr_C2, fr_G),
   # pulse = (; start=pulse_start, end=pulse_end, strength=p_strength, G_pre, G_pulse),
    idx_post_min,
    t_post_min,
    raw = (; t=t_raw, fr_total=fr_total_raw, fr_R1=fr_R1_raw, fr_R2=fr_R2_raw,
           fr_C1=fr_C1_raw, fr_C2=fr_C2_raw, fr_G=fr_G_raw)
)
    
end

##function to plot timeseries:
function plot_pulse_timeseries(pulse; 
    tmin=175.0,
    ylims=(0.0, 1.0),
    use_raw=false,
    colors = (:black, :darkgreen, :salmon, :lightgreen, :pink, :blue),
    show_pulse=true
)
    if use_raw
        t   = pulse.raw.t
        ft  = pulse.raw.fr_total
        fR1 = pulse.raw.fr_R1
        fR2 = pulse.raw.fr_R2
        fC1 = pulse.raw.fr_C1
        fC2 = pulse.raw.fr_C2
        fG  = pulse.raw.fr_G
    else
        t   = pulse.t_grid
        ft  = pulse.fr_total
        fR1 = pulse.fr_components.fr_R1
        fR2 = pulse.fr_components.fr_R2
        fC1 = pulse.fr_components.fr_C1
        fC2 = pulse.fr_components.fr_C2
        fG  = pulse.fr_components.fr_G
    end

    m = t .>= tmin
    tt  = t[m]
    ft  = ft[m]
    fR1 = fR1[m]; fR2 = fR2[m]
    fC1 = fC1[m]; fC2 = fC2[m]
    fG  = fG[m]

    plt = plot(tt, ft;  label="total → P", xlabel="Time", ylabel="Community Consumption",
               ylims=ylims, color=colors[1], linewidth=2.5, legend=true, framestyle=:box)
    plot!(tt, fR1; label="R1 → P", color=colors[2], linewidth=2.5)
    plot!(tt, fR2; label="R2 → P", color=colors[3], linewidth=2.5)
    plot!(tt, fC1; label="C1 → P", color=colors[4], linewidth=2.5)
    plot!(tt, fC2; label="C2 → P", color=colors[5], linewidth=2.5)
    plot!(tt, fG;  label="G → P",  color=colors[6], linewidth=2.5)
   # plot!(tt, fG;  label="total → P",  color=colors[1], linewidth=2.5)

    if show_pulse && (:pulse in propertynames(pulse))
        ps = pulse.pulse.start
        pe = pulse.pulse.end
        vspan!([ps, pe]; alpha=0.1, label="", color=:gray)
    end
    return plt
end

##function to plot timeseries so have proper colors for total harvest when only groceries
function plot_pulse_timeseries_2(pulse; 
    tmin=175.0,
    ylims=(0.0, 1.0),
    use_raw=false,
    colors = (:black, :darkgreen, :salmon, :lightgreen, :pink, :blue),
    show_pulse=true
)
    if use_raw
        t   = pulse.raw.t
        ft  = pulse.raw.fr_total
        fR1 = pulse.raw.fr_R1
        fR2 = pulse.raw.fr_R2
        fC1 = pulse.raw.fr_C1
        fC2 = pulse.raw.fr_C2
        fG  = pulse.raw.fr_G
    else
        t   = pulse.t_grid
        ft  = pulse.fr_total
        fR1 = pulse.fr_components.fr_R1
        fR2 = pulse.fr_components.fr_R2
        fC1 = pulse.fr_components.fr_C1
        fC2 = pulse.fr_components.fr_C2
        fG  = pulse.fr_components.fr_G
    end

    m = t .>= tmin
    tt  = t[m]
    ft  = ft[m]
    fR1 = fR1[m]; fR2 = fR2[m]
    fC1 = fC1[m]; fC2 = fC2[m]
    fG  = fG[m]

    plt = plot(tt, ft;  label="total → P", xlabel="Time", ylabel="Community Consumption",
               ylims=ylims, color=colors[1], linewidth=2.5, legend=false, framestyle=:box)
    plot!(tt, fR1; label="R1 → P", color=colors[2], linewidth=2.5)
    plot!(tt, fR2; label="R2 → P", color=colors[3], linewidth=2.5)
    plot!(tt, fC1; label="C1 → P", color=colors[4], linewidth=2.5)
    plot!(tt, fC2; label="C2 → P", color=colors[5], linewidth=2.5)
    plot!(tt, fG;  label="G → P",  color=colors[6], linewidth=2.5)
    plot!(tt, ft;  label="total → P",  color=colors[1], linewidth=2.5)

    if show_pulse && (:pulse in propertynames(pulse))
        ps = pulse.pulse.start
        pe = pulse.pulse.end
        vspan!([ps, pe]; alpha=0.1, label="", color=:gray)
    end
    return plt
end

p = ModelPar_active()
test = pulse_unit_forced(p)

##Calculate decline for 3 different scenarios
p1 = ModelPar_active(o = 0.1, w = 0.5, H = 0.5, K = 3.05, G_base = 1.0)
pulse_1 = pulse_unit_forced(p1)

p2 = ModelPar_active(o = 0.0, w = 1.0, H = 0.5, K = 3.05, G_base = 1.0)
pulse_2 = pulse_unit_forced(p2)


p3 = ModelPar_active(o = 0.0, w = 0.0, H = 1.0, K = 3.05, G_base = 1.0)
pulse_3 = pulse_unit_forced(p3)

plot_pulse_timeseries_2(pulse_1; tmin=150.0)
plot_pulse_timeseries_2(pulse_2; tmin=150.0)
plot_pulse_timeseries_2(pulse_3; tmin=150.0)

plot_pulse_timeseries(pulse_1; tmin=150.0)






###Plot relative decline
# Pull the relative declines (dimensionless)
declines = [
    pulse_1.decline_rel_baseline,
    pulse_2.decline_rel_baseline,
    pulse_3.decline_rel_baseline
]

labels = [
    "o=0.1, w=0.5, H=0.2, K=3.05",
     "o=0.0, w=1.0, H=0.2, K=3.05",
    "o=0.0, w=0.0, H=1.0, K=3.05"
]

# Convert to % for display (optional but nice)
declines_pct = (100 .* declines)
resilience = 1 ./ declines
# Bar plot (no legend)
res_plot = bar(labels, declines_pct;
    legend = false,
    ylabel = "% Food Consumption Decline",
    xlabel = "",
    framestyle = :box,
    xrotation = 15,
    yticks = :auto, ylims = (0.0, 40),  color = :black, size=(600, 800))

savefig(res_plot, "resilience_plot")


###Seeing if i can create plot of decline relative to baseline for all values of H 
results_H_all_f_1 = []
for H in 0.1:0.1:1.0
    p = ModelPar_active(w = 0.0, o = 0.0, H=H, l = 0.5)
    P0 = 0.25
    eq_data = pulse_unit_forced(p)
    push!(results_H_all_f_1, (; H=H, eq_data...))
end
df_pulse_1 = DataFrame(results_H_all_f_1)


rel_decline_1 = [row.decline_rel_baseline[1] for row in eachrow(df_pulse_1)]
plot(df_pulse_1.H, rel_decline_1, col = "red", ylabel = "Relative decline in harvest after perturbation", xlabel = "Grocery Preference (H)", label = "o = 0, w = 0")


results_H_all_f_2 = []
for H in 0.1:0.1:1.0
    p = ModelPar_active(w = 0.5, o = 0.1, H=H, l = 0.5)
    P0 = 0.25
    eq_data = pulse_unit_forced(p)
    push!(results_H_all_f_2, (; H=H, eq_data...))
end
df_pulse_2 = DataFrame(results_H_all_f_2)


rel_decline_2 = [row.decline_rel_baseline[1] for row in eachrow(df_pulse_2)]
plot(df_pulse_2.H, rel_decline_2,  ylabel = "Relative decline in harvest after perturbation", xlabel = "H", label = "o = 0.1, w = 0.5")


results_H_all_f_3 = []
for H in 0.1:0.1:1.0
    p = ModelPar_active(w = 1.0, o = 0.0, H=H, l = 0.5)
    P0 = 0.25
    eq_data = pulse_unit_forced(p)
    push!(results_H_all_f_3, (; H=H, eq_data...))
end
df_pulse_3 = DataFrame(results_H_all_f_3)


rel_decline_3 = [row.decline_rel_baseline[1] for row in eachrow(df_pulse_3)]
plot(df_pulse_3.H, rel_decline_3,  ylabel = "Relative decline in harvest after perturbation", xlabel = "H", label = "o = 0.0, w = 1.0")

##trying to see if i can plot all 3 together
plot(df_pulse_1.H, rel_decline_1, color = :black, linestyle = :solid, ylabel = "Relative decline in harvest \nafter perturbation", xlabel = "Grocery Preference (H)", label = "o = 0, w = 0", linewidth = 2.5, size=(600, 800))
plot!(df_pulse_2.H, rel_decline_2, color = :black, linestyle = :dash,  ylabel = "Relative decline in harvest \nafter perturbation", xlabel = "Grocery Preference (H)", label = "o = 0.1, w = 0.5", linewidth = 2.5, size=(600, 800))
plot!(df_pulse_3.H, rel_decline_3,  color = :black, linestyle =  :dashdotdot, ylabel = "Relative decline in harvest \nafter perturbation", xlabel = "Grocery Preference (H)", label = "o = 0.0, w = 1.0", linewidth = 2.5, size=(600, 800))



##sample code for doing burn in stuff. ...
# 1) Burn-in with constant G_pre to t = pulse_start
p_const = deepcopy(p); p_const.G_func = t -> G_pre
prob1 = ODEProblem(model_forced!, u0, (0.0, pulse_start), p_const, tstops=[pulse_start])
sol1  = solve(prob1; reltol=1e-8, abstol=1e-8)
u_at_pulse = sol1(pulse_start - 1e-9)  # state just before pulse

# 2) Now run the full piecewise G (including pulse window) from t=pulse_start to t_end
p_pulse = deepcopy(p)
p_pulse.G_func = t -> (t < pulse_start || t >= pulse_end) ? G_pre : G_pulse
prob2 = ODEProblem(model_forced!, u_at_pulse, (pulse_start, t_end), p_pulse, tstops=[pulse_start, pulse_end])
sol2  = solve(prob2; reltol=1e-8, abstol=1e-8)

# 3) Stitch burn-in + pulse for plotting if desired
t_all = [sol1.t; sol2.t]
u_all = [sol1.u; sol2.u]



 ###OLD CODE   


tspan = (0.0, 500.0)
G_pre = 1.0 
G_pulse = G_pre
t_pulse = 200.0 ##time when disturbance occurs, want to be once model at equilibirum
t_recover = 250.0

# Perform the pulse for the 3 systems: 
# - Food chain: "chain";C OP (fixed)
# - Passive omnivory: "passive";
# - Active omnivory: "active";
##this structure is from KC, am keeping some of this in so that if i find way to streamline between the different structures i can run this pulse simulataneously 
function pulse(base_par, p_length, p_strength) 
    
   # par_chain = deepcopy(base_par)
   # par_chain.Ω = 0.0
   # par_passive = deepcopy(base_par)
   # par_passive.Ω = Ω
   # par_passive.pref = fixed_pref
   # par_active = deepcopy(par_passive)
   # par_active.pref = adapt_pref
    
   par_active = deepcopy(base_par)

    push!(
        [],
      #  pulse_unit(par_chain, p_length, p_strength),
      #  pulse_unit(par_passive, p_length, p_strength),
        pulse_unit_forced(par_active, P_fixed, p_length, p_strength)
    )
end


##run pulse experiment 
include("wf_model_eqs_subsidy_Pfixed.jl")
par = ModelPar_active()
P_fixed = 0.25
res = pulse(par, 2.0, -2.0)


plot





##OLD - other method based on KC code -- that calculates eigenvalues etc. but not needed for this 
# Entire analysis for one pulse simulation (1 system)
function pulse_unit_forced(par, P_fixed, p_length, p_strength) 
    # Parameters   
    u0 = [1.0, 1.0, 1.5, 1.5, P_fixed]            #initial state of R1, R2, C1, C1 and Pfixed -- define P_fixed outside of function
    pulse_start = 200.0               ##time when pulse begins
    pulse_end = pulse_start + p_length      #time when pulse ends
   # pulse_event_times = union(pulse_start, pulse_end) #times where the pulse changes system
    t_end = 350.0         #end of simulation
    t_span = (0.0, t_end)       #time span of simulation
    t_start = 175.0         #start of post-pulse tracking (i.e., when tracking starts to capture pulse and dynamics after)
    len = 1000           ##resolution of the time grid
    t_grid = range(t_start, t_end, length = len)        #time grid for evaluating solution
    t_after = range(pulse_end, t_end, length = len)     #time grid after the pulse
   
    #I think define par outside of function ... 

    ##Pulse Event Logic ##
    # functions
    #trigger pulse events at pulse_start and pulse_end
    function pulse_event(u, t, integrator) 
        t == pulse_start || t == pulse_end
    end
    
    #at pulse_start: apply pulse by scaling G
    #at pulse_end: restore original G
    # NB: for pulse p_strength is a multiplicator
    function forcing_pulse!(integrator)
        if integrator.t == pulse_start
            # multiplication for pulse
            integrator.p.G = p_strength * integrator.p.G_base  ##this is basically telling it to multiply G by whatever p_strength is, once hit time of pulse_start
        elseif integrator.t == pulse_end
            integrator.p.G = integrator.p.G_base  ##after pulse ends go back to base G
        end
        return
    end
    
    #callback to apply the pulse
    cb_pulse = DiscreteCallback(pulse_event, forcing_pulse!)


    ##Simulate Pulse for a Given parameter set
    par_pulse = deepcopy(par)
    par_pulse.G = p_strength * par.G_base  ##initial pulse value (will be reset via callback)
    prob = ODEProblem(model_forced!, u0, t_span, deepcopy(par), tstops = [pulse_start, pulse_end])
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb_pulse)
    sol_grid = sol(t_grid)
    
    # --- equilibrium before pulse: use the state just before pulse_start ---
    u_pre = sol(pulse_start - 1e-6)            # 5-vector state
    eq_before = find_eq_forced(u_pre, par)     # returns 5-vector [R1,R2,C1,C2,P]

    # --- community matrix for the 4 dynamic variables with P fixed ---
    M = cmat_forced(eq_before[1:4], par, eq_before[5])  # <- (4x4) square

    λ1      = λ1_stability(M)
    λ1_imag = λ1_stability_imag(M)
    react   = ν_stability(M)

    #Find the time after pulse when each species hits equilibrium again --- 
    t_hit_eq = find_times_hit_equil_press(sol(t_after))
    # println(t_hit_eq)
    
    #Overshoot: max deviation from equilibrium after pulse
    R1_OS = overshoot(sol, eq_before, 1, t_hit_eq[1], t_end)
    R2_OS = overshoot(sol, eq_before, 2, t_hit_eq[2], t_end)
    C1_OS = overshoot(sol, eq_before, 3, t_hit_eq[3], t_end)
    C2_OS = overshoot(sol, eq_before, 4, t_hit_eq[4], t_end)
    # Calculate max-min metric - max-min amplitude of post-pulse dynamics
    R1_mm = min_max(sol, 1, t_hit_eq[1], t_end)
    R2_mm = min_max(sol, 2, t_hit_eq[2], t_end)
    C1_mm = min_max(sol, 3, t_hit_eq[3], t_end)
    C2_mm = min_max(sol, 4, t_hit_eq[4], t_end)
    ##return all metrics in a vector
    push!(
        [],
        λ1,
        -1 ./ λ1,
        λ1_imag,
        react,
        R1_OS,
        R2_OS, 
        C1_OS, 
        C2_OS,
        R1_mm, 
        R2_mm, 
        C1_mm, 
        C2_mm,
        eq_before,  # for consistency in output size with `press()`
        eq_before,   
        par,
        sol
    )
    
end