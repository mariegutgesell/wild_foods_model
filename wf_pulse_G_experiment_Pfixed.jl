##Function to run pulse experiment -- modifying approach from KC


using Pkg
#Pkg.add("NLsolve")


#using DifferentialEquations
using DifferentialEquations, NLsolve


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
    len = 100000            ##resolution of the time grid
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