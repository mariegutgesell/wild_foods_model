##Function to run pulse experiment -- modifying approach from KC


using Pkg
#Pkg.add("NLsolve")


#using DifferentialEquations
using DifferentialEquations, NLsolve


# Entire analysis for one pulse simulation (1 system)
function pulse_unit_forced(par, P_fixed, p_length, p_strength) 
    # Parameters   
    u0 = [1.0, 1.0, 1.5, 1.5, P_fixed]            #initial state of R1, R2, C1, C1 and Pfixed -- define P_fixed outside of function
    pulse_start = 200               ##time when pulse begins
    pulse_end = pulse_start + p_length      #time when pulse ends
    pulse_event_times = union(pulse_start, pulse_end) #times where the pulse changes system
    t_end = 350         #end of simulation
    t_span = (0.0, t_end)       #time span of simulation
    t_start = 175.0         #start of post-pulse tracking (i.e., when tracking starts to capture pulse and dynamics after)
    len = 100000            ##resolution of the time grid
    t_grid = range(t_start, t_end, length = len)        #time grid for evaluating solution
    t_after = range(pulse_end, t_end, length = len)     #time grid after the pulse
   
    #I think define par outside of function ... 

    ##Pulse Event Logic ##
    # functions
    #trigger pulse events at pulse_start and pulse_end
    pulse_event(u, t, integrator) = t ∈ pulse_event_times
    
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
    prob = ODEProblem(model_forced!, u0, t_span, deepcopy(par), tstops = pulse_event_times)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb_pulse)
    sol_grid = sol(t_grid)
    
    #Find the time after pulse when each species hits equilibrium again --- where is that function coming from? in press code?
    t_hit_eq = find_times_hit_equil_press(sol(t_after))
    # println(t_hit_eq)
    
    # Same equilibrium before and after
    eq_before = find_eq(sol_grid[1], par)   ##equilibrium before pulse

    # evaluate on 
    λ1 = λ1_stability(cmat(eq_before, par)) #real part of leading eigenvalue
    λ1_imag = λ1_stability_imag(cmat(eq_before, par)) ##imaginary part (oscillation potential)
    react = ν_stability(cmat(eq_before, par)) ##reactivity (initial response to perturbation)

    #Overshoot: max deviation from equilibrium after pulse
    predator_OS = overshoot(sol, eq_before, 3, t_hit_eq[3], t_end)
    consumer_OS = overshoot(sol, eq_before, 2, t_hit_eq[2], t_end)
    resource_OS = overshoot(sol, eq_before, 1, t_hit_eq[1], t_end)

    # Calculate max-min metric - max-min amplitude of post-pulse dynamics
    predator_mm = min_max(sol, 3, t_hit_eq[3], t_end)
    consumer_mm = min_max(sol, 2, t_hit_eq[2], t_end)
    resource_mm = min_max(sol, 1, t_hit_eq[1], t_end)

    ##return all metrics in a vector
    push!(
        [],
        λ1,
        -1 ./ λ1,
        λ1_imag,
        react,
        predator_OS,
        consumer_OS, 
        resource_OS, 
        predator_mm, 
        consumer_mm, 
        resource_mm, 
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
function pulse(base_par, Ω, p_length, p_strength) 
    
    par_chain = deepcopy(base_par)
    par_chain.Ω = 0.0
    par_passive = deepcopy(base_par)
    par_passive.Ω = Ω
    par_passive.pref = fixed_pref
    par_active = deepcopy(par_passive)
    par_active.pref = adapt_pref
    
    push!(
        [],
        pulse_unit(par_chain, p_length, p_strength),
        pulse_unit(par_passive, p_length, p_strength),
        pulse_unit(par_active, p_length, p_strength)
    )

end
