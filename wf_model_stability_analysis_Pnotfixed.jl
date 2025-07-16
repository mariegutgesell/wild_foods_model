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

##calculate local stability (max eigenvalue across values of K)
results_K = []

for K in 0.1:0.1:3.5
    λmax = NaN
    C1_val = NaN
    C2_val = NaN
    try
        p = ModelPar_passive(K = K,
        H = 0, 
        w = 0,
        o = 0)
        x0 = [1.0, 1.0, 0.5, 0.5, 0.5]

        result = nlsolve(x -> rhs_2(x, p), x0)

        if result.f_converged
            x_star = result.zero
            R1, R2, C1_tmp, C2_tmp, P = x_star
            C1_val = C1_tmp
            C2_val = C2_tmp
            J = ForwardDiff.jacobian(x -> rhs_2(x, p), x_star)
            λmax = maximum(real.(eigvals(J)))
        else
            @warn "Solver did not converge for K=$K"
        end
    catch err
        @warn "Error at K=$K: $err"
    end

    push!(results_K, (K = K, λmax = λmax, C1 = C1_val, C2 = C2_val))
end

df_K = DataFrame(results_K)   

plot(df_K.K, df_K.λmax)





# Grid of values for o and w

##Calculating local equilibrium and stability (eigenvalues)
# Initial guess for dynamic variables: R1, R2, C1, C2, P
x0 = [2.0, 2.0, 1.0, 1.0, 1.0]
p = ModelPar_passive(o = 0.75,
                     w = 0.5)

##Local stability analysis
#1) Calculate equilibrium values for system of equation 
result = nlsolve(x -> rhs_2(x, p), x0)

# Extract equilibrium values
x_star = result.zero
R1_star, R2_star, C1_star, C2_star, P_star = x_star
 

#2) Calculate numerical jacobian
J_num = ForwardDiff.jacobian(x -> rhs_2(x, p), x0)
println("Jacobian Matrix:")
display(J_num)

#3) Calculate the eigenvalues 
eigvals(J_num)


##Trying to find multiple equilibriums, as likely more than 1 
#Pkg.add("IterTools")
using IterTools

R_vals = [0.1, 0.5, 1.0, 2.0]
C_vals = [0.1, 0.5, 1.0]
P_vals = [0.1, 0.5, 1.0]
# All combinations of R1, R2, C1, C2
guesses = collect(product(R_vals, R_vals, C_vals, C_vals, P_vals))


tolerance = 1e-4
equilibria_1 = []

for x0 in guesses
    result = nlsolve(x -> rhs_2(x,p), collect(x0))
    if result.f_converged
        x_star = result.zero
        # Check if it's a new equilibrium
        if all(norm(x_star .- eq) > tolerance for eq in equilibria_1)
            push!(equilibria_1, x_star)
        end
    end
end


max_eigenvalues_1 = Float64[]  # to store max real parts

println("Number of distinct equilibria found: ", length(equilibria_1))

for (i, eq) in enumerate(equilibria_1)
    J = ForwardDiff.jacobian(x -> rhs_2(x, p), eq)
    λ = eigvals(J)
    λ_real = real.(λ)

    # For determining stability
    is_stable = all(λ_real .< 0)
    stability = is_stable ? "Stable" : "Unstable or saddle"

    # For assessing strength of stability
    λ_max = maximum(λ_real)  # still needed for stability check
    λ_min = minimum(λ_real)  # use this to define "most negative" eigenvalue

    # If stable, response time based on most negative direction
    RT = is_stable ? 1 / abs(λ_min) : NaN

    println("\nEquilibrium $i:")
    println("  State: ", round.(eq; digits=4))
    println("  Stability: ", stability)
    println("  Max real eigenvalue: ", round(λ_max; digits=6))
    println("  Min real eigenvalue: ", round(λ_min; digits=6))
    println("  RT (based on min eig): ", round(RT; digits=4))
    println("  All eigenvalues: ", round.(λ; digits=4))

    push!(max_eigenvalues_1, λ_max)

end


##tracking oscillatory equilibrium across values of o and w 
using LinearAlgebra

#p = ModelPar_passive()
p = ModelPar_active()
R_vals = [0.1, 0.5, 1.0, 2.0]
C_vals = [0.1, 0.5, 1.0]
P_vals = [0.1, 0.5, 1.0]
# All combinations of R1, R2, C1, C2
guesses = collect(product(R_vals, R_vals, C_vals, C_vals, P_vals))

tolerance = 1e-4
rows = []

for o in 0.0:0.1:1.0
    for w in 0.0:0.1:1.0
        println("\n--- Checking o = $o, w = $w ---")
        p = ModelPar_active(
            w = w,
            o = o,
            # (include all other required fields here or use a helper)- do i need to add the others in or is this calling from right p? 
        )

        equilibria_stable = []

        for x0 in guesses
            result = nlsolve(x -> rhs_2(x, p), collect(x0))
            if result.f_converged
                x_star = result.zero

                if any(isnan.(x_star)) || any(x_star .< 0)
                    continue
                end

                # Check if distinct
                if all(norm(x_star .- eq) > tolerance for eq in equilibria_stable)
                    J = ForwardDiff.jacobian(x -> rhs_2(x, p), x_star)
                    λ = eigvals(J)
                    λ_real = real.(λ)
                    λ_imag = imag.(λ)

                    if all(λ_real .< 0)
                        push!(equilibria_stable, x_star)

                        has_oscillation = any(abs.(λ_imag) .> 1e-6)
                        stability_type = has_oscillation ? "Stable spiral (oscillatory)" : "Stable node"
                        λ_max = maximum(λ_real)
                        λ_min = minimum(λ_real)
                        RT = 1 / abs(λ_min)

                        push!(rows, (
                            o = o,
                            w = w,
                            stability = stability_type,
                            lambda_max = λ_max,
                            lambda_min = λ_min,
                            RT = RT,
                            oscillatory = has_oscillation,
                            C1 = x_star[3], 
                            C2 = x_star[4], # optional: include state vars
                            P = x_star[5]
                        ))
                    end
                end
            end
        end
    end
end

df_stability = DataFrame(rows)


o_vals = unique(df_stability.o)
w_vals = unique(df_stability.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df_stability[(df_stability.o .== o) .& (df_stability.w .== w), :lambda_max][1] for w in w_vals, o in o_vals]
λ_mat = [
    begin
        subset = df_stability[(df_stability.o .== o) .& (df_stability.w .== w), :lambda_max]
        !isempty(subset) ? subset[1] : NaN
    end
    for w in w_vals, o in o_vals
]
# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)


df_stability_o = filter(row -> row.w == 0.0, df_stability)   
plot(df_stability_o.o, df_stability_o.lambda_max, 
    xlabel = "Omnivory Preference (o)", ylabel = "Max Real Eigenvalue")

df_stability_w = filter(row -> row.o == 0.0, df_stability)   

plot(df_stability_w.w, df_stability_w.lambda_max, 
xlabel = "Habitat Preference (o)", ylabel = "Max Real Eigenvalue")

##Plotting P 
P_mat = [
    begin
        subset = df_stability[(df_stability.o .== o) .& (df_stability.w .== w), :P]
        !isempty(subset) ? subset[1] : NaN
    end
    for w in w_vals, o in o_vals
]
# Plot heatmap
heatmap(o_vals, w_vals, P_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Predator Density",
        colorbar_title = "P",
        c = :viridis)

##Plotting C 
C1_mat = [
    begin
        subset = df_stability[(df_stability.o .== o) .& (df_stability.w .== w), :C1]
        !isempty(subset) ? subset[1] : NaN
    end
    for w in w_vals, o in o_vals
]
# Plot heatmap
heatmap(o_vals, w_vals, C1_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "C1 Density",
        colorbar_title = "C",
        c = :viridis)

##Plotting C2 
C2_mat = [
    begin
        subset = df_stability[(df_stability.o .== o) .& (df_stability.w .== w), :C2]
        !isempty(subset) ? subset[1] : NaN
    end
    for w in w_vals, o in o_vals
]
# Plot heatmap
heatmap(o_vals, w_vals, C2_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "C2 Density",
        colorbar_title = "C2",
        c = :viridis)


##QUESTIONS:
##Why getting NAs for some combinations of o and w? 
##Some of the o dynamics look weird...is this indicative of a bifurcation? 
##


##Looking across range of w 
w_vals = 0.0:0.05:1.0
results = []

x0_guesses = [
    [1.0, 1.0, 0.5, 0.5, 0.5],
    [2.0, 2.0, 1.0, 1.0, 1.0],
    [0.5, 0.5, 0.2, 0.2, 0.3]
]

for w in w_vals
    found = false
    λ_min = NaN
    RT = NaN

    p = ModelPar_passive(w=w)

    for x0 in x0_guesses
        try
            result = nlsolve(x -> rhs_2(x, p), x0)

            if result.f_converged
                x_star = result.zero

                if all(x_star .> 0) && all(isfinite, x_star)
                    J = ForwardDiff.jacobian(x -> rhs_2(x, p), x_star)

                    if all(isfinite, J)
                        λ_real = real.(eigvals(J))

                        if all(λ_real .< 0)
                            λ_min = minimum(λ_real)
                            RT = 1 / abs(λ_min)
                            found = true
                            break
                        end
                    end
                end
            end
        catch
            # ignore and continue to next guess
        end
    end

    push!(results, (w=w, λ_min=λ_min, RT=RT, stable=found))
end

df = DataFrame(results)

plot(df.w, df.RT,
    xlabel="Habitat Preference (w)",
    ylabel="Return Time (RT)",
    title="Return Time vs Habitat Coupling",
    legend=false,
    marker=:circle,
    lw=2)

plot(df.w, df.λ_min,
    xlabel="Habitat Preference (w)",
    ylabel="Minimum Eigenvalue",
    title="Minimum eigenvalue vs. Habitat Preference",
    legend=false,
    marker=:circle,
    lw=2)





##looking at bifurcation plot of pw_vals = 0.0:0.01:1.0
##this function just extracts the steady state equilibrium values, not dynamic behaviour from numerical integration -- plots where system settles not oscillation around that point 

x0_guesses = [
    [1.0, 1.0, 0.5, 0.5, 0.5],
    [2.0, 2.0, 1.0, 1.0, 1.0],
    [0.5, 0.5, 0.2, 0.2, 0.3]
]

bifurcation_data = []

for w in w_vals
    p = ModelPar_passive(w=w)
    found = false

    for x0 in x0_guesses
        try
            result = nlsolve(x -> rhs_2(x, p), x0)
            if result.f_converged
                x_star = result.zero

                if all(x_star .> 0) && all(isfinite, x_star)
                    J = ForwardDiff.jacobian(x -> rhs_2(x, p), x_star)
                    if all(isfinite, J)
                        λ_real = real.(eigvals(J))
                        is_stable = all(λ_real .< 0)
                        R1, R2, C1, C2, P = x_star
                        push!(bifurcation_data, (w=w, R1=R1, R2=R2, C1=C1, C2=C2, P=P, stable=is_stable))
                        found = true
                        break
                    end
                end
            end
        catch
            # skip and try next x0
        end
    end

    if !found
        push!(bifurcation_data, (w=w, R1=NaN, R2=NaN, C1=NaN, C2=NaN, P=NaN, stable=false))
    end
end

df_bif = DataFrame(bifurcation_data)


plot(df_bif.w, df_bif.P,
    xlabel = "Habitat Preference (w)",
    ylabel = "Predator Biomass (P*)",
    title = "Bifurcation Diagram",
    label = "",
    markershape = :circle,
    color = :black)

plot(df_bif.w, df_bif.C1,
    xlabel = "Habitat Preference (w)",
    ylabel = "Consumer Biomass (C1*)",
    title = "Bifurcation Diagram",
    label = "",
    markershape = :circle,
    color = :black)


plot(df_bif.w, df_bif.C2,
    xlabel = "Habitat Preference (w)",
    ylabel = "Consumer Biomass (C2*)",
    title = "Bifurcation Diagram",
    label = "",
    markershape = :circle,
    color = :black)



##trying plotting min and max of C1 - integrating across numeric output with different dynamics
function test_one_bif_point(w)
    p = ModelPar_passive(
        w = w,
        o = 0,
        H = 0,
        r = 2.0,
        K = 3.0,
        aR_C = 1.0,
        aR_P = 0.2,
        aC_P = 0.5,
        aG_P = 0.5,
        e = 0.5,
        mC = 0.4,
        mP = 0.2,
        hR_C = 0.4,
        hR_P = 0.6,
        hC_P = 0.3,
        hG_P = 0.3,
        W = hab_pref,
        d_om_i = om_i_pref_fixed,
        d_om_j = om_j_pref_fixed,
        sub_pref = sub_pref_func,
        G = 0,
        f_r1c1 = f_R1C1,
        f_r2c2 = f_R2C2,
        f_r1p = f_R1P,
        f_r2p = f_R2P,
        f_c1p = f_C1P,
        f_c2p = f_C2P,
        f_gp = f_GP,
        f_gp2 = f_GP_2,
        l1 = 1.0,
        l2 = 1.0,
        pf = 10.0,
        D = 0.5,
        e1 = t -> sin(2π / 10.0 * t),
        e2 = t -> sin(2π / 10.0 * (t - 0.5 * 10.0))
    )

    x0 = [0.6, 0.8, 0.45, 0.65, 0.2]
    tspan = (0.0, 1000.0)
    prob = ODEProblem(rhs_2, x0, tspan, p)

    function condition(u, t, integrator)
        any(u .< 0)
    end

    function affect!(integrator)
        integrator.u .= max.(integrator.u, 1e-6)
    end

    cb = DiscreteCallback(condition, affect!)
    sol = solve(prob, callback=cb, reltol=1e-8, abstol=1e-8)

    println("Solver status: ", sol.retcode)
    println("Final time: ", sol.t[end])
    t_eval = sol.t[end] - 100 : 0.1 : sol.t[end]
    C1_vals = sol(t_eval)[3, :]
    println("C1 min = ", minimum(C1_vals))
    println("C1 max = ", maximum(C1_vals))

    C1_min = minimum(C1_vals)
    C1_max = maximum(C1_vals)

    return DataFrame(w = w, C1_min = C1_min, C1_max = C1_max)
end

test_one_bif_point(1.0)

w_vals = 0.0:0.1:1.0

df_bif = for w in w_vals 
    df_bif = test_one_bif_point(w)
return df_bif

end


function run_bifurcation_multi(test_one_bif_point, w_vals)
    dfs = []

    for w in w_vals
        df = test_one_bif_point(w)
        push!(dfs, df)
    end

    return vcat(dfs...)
end

df_bif = run_bifurcation_multi(test_one_bif_point, 0.0:0.1:1.0)


plot(df_bif.w, df_bif.C1_min,
    label = "C1 min", xlabel = "Habitat preference (w)",
    ylabel = "C1 biomass",
    title = "C1 oscillation envelope",
    lw = 2, color = :blue)

plot!(df_bif.w, df_bif.C1_max,
    label = "C1 max", linestyle = :dash, lw = 2, color = :red)

##testing by just running one point:
w = 1.0
x0 = [0.6, 0.8, 0.45, 0.65, 0.2]
p = ModelPar_passive(w=w)

prob = ODEProblem(rhs_2, x0, (0.0, 1000.0), p)
sol = solve(prob, callback = cb, reltol=1e-8, abstol=1e-8)

plot(sol.t, sol[3, :], label="C1", xlabel="Time", ylabel="Biomass")


##Looking across range of o 
##trying plotting min and max of C1 - integrating across numeric output with different dynamics
function test_one_bif_point_o(o)
    p = ModelPar_passive(
        w = 0.5,
        o = o,
        H = 0,
        r = 2.0,
        K = 3.0,
        aR_C = 1.0,
        aR_P = 0.2,
        aC_P = 0.5,
        aG_P = 0.5,
        e = 0.5,
        mC = 0.4,
        mP = 0.2,
        hR_C = 0.4,
        hR_P = 0.6,
        hC_P = 0.3,
        hG_P = 0.3,
        W = hab_pref,
        d_om_i = om_i_pref_fixed,
        d_om_j = om_j_pref_fixed,
        sub_pref = sub_pref_func,
        G = 0,
        f_r1c1 = f_R1C1,
        f_r2c2 = f_R2C2,
        f_r1p = f_R1P,
        f_r2p = f_R2P,
        f_c1p = f_C1P,
        f_c2p = f_C2P,
        f_gp = f_GP,
        f_gp2 = f_GP_2,
        l1 = 1.0,
        l2 = 1.0,
        pf = 10.0,
        D = 0.5,
        e1 = t -> sin(2π / 10.0 * t),
        e2 = t -> sin(2π / 10.0 * (t - 0.5 * 10.0))
    )

    x0 = [0.6, 0.8, 0.45, 0.65, 0.2]
    tspan = (0.0, 1000.0)
    prob = ODEProblem(rhs_2, x0, tspan, p)

    function condition(u, t, integrator)
        any(u .< 0)
    end

    function affect!(integrator)
        integrator.u .= max.(integrator.u, 1e-6)
    end

    cb = DiscreteCallback(condition, affect!)
    sol = solve(prob, callback=cb, reltol=1e-8, abstol=1e-8)

    println("Solver status: ", sol.retcode)
    println("Final time: ", sol.t[end])
    t_eval = sol.t[end] - 100 : 0.1 : sol.t[end]
    C1_vals = sol(t_eval)[3, :]
    println("C1 min = ", minimum(C1_vals))
    println("C1 max = ", maximum(C1_vals))

    C1_min = minimum(C1_vals)
    C1_max = maximum(C1_vals)

    return DataFrame(o = o, C1_min = C1_min, C1_max = C1_max)
end

test_one_bif_point_o(1.0)

o_vals = 0.0:0.1:1.0

df_bif_o = for o in o_vals 
    df_bif_o = test_one_bif_point(o)
return df_bif_o

end


function run_bifurcation_multi_o(test_one_bif_point_o, o_vals)
    dfs = []

    for o in o_vals
        df = test_one_bif_point_o(o)
        push!(dfs, df)
    end

    return vcat(dfs...)
end

df_bif_o = run_bifurcation_multi_o(test_one_bif_point_o, 0.0:0.1:1.0)


plot(df_bif_o.o, df_bif_o.C1_min,
    label = "C1 min", xlabel = "Omnivory preference (o)",
    ylabel = "C1 biomass",
    title = "C1 oscillation envelope",
    lw = 2, color = :blue)

plot!(df_bif_o.o, df_bif_o.C1_max,
    label = "C1 max", linestyle = :dash, lw = 2, color = :red)




    

##OLD CODE BELOW

results = []

# Define multiple biologically plausible initial guesses
x0_guesses = [
    [1.0, 1.0, 0.5, 0.5, 0.5],
    [2.0, 2.0, 1.0, 1.0, 1.0],
    [0.5, 0.5, 0.2, 0.2, 0.3],
    [1.5, 1.5, 0.8, 0.8, 0.6]
]

# Parameter grid
for o in 0.0:0.1:1.0, w in 0.0:0.1:1.0
    λmax = NaN
    RT = NaN
    C1 = NaN
    C2 = NaN

    found = false
    p = ModelPar_passive(o=o, w=w)

    for x0 in x0_guesses
        try
            result = nlsolve(x -> rhs_2(x, p), x0)

            if result.f_converged
                x_star = result.zero

                if all(x_star .> 0) && all(isfinite, x_star)
                    J = ForwardDiff.jacobian(x -> rhs_2(x, p), x_star)

                    if all(isfinite, J)
                        λs = eigvals(J)
                        λ_real = real.(λs)
                        λmax = maximum(λ_real)
                        RT = 1 / abs(λmax)
                        _, _, C1, C2, _ = x_star

                        push!(results, (o=o, w=w, λmax=λmax, RT=RT, C1=C1, C2=C2))
                        found = true
                        break  # Stop looping over x0_guesses
                    end
                end
            end
        catch err
            @warn "Error at o=$o, w=$w, x0=$x0: $err"
        end
    end

    # If no feasible equilibrium found, store NA row
    if !found
        push!(results, (o=o, w=w, λmax=NaN, RT=NaN, C1=NaN, C2=NaN))
    end
end

df = DataFrame(results)

##looking across ranges of o and w 

# Get unique values
o_vals = unique(df.o)
w_vals = unique(df.w)

# Sort them to be safe
sort!(o_vals)
sort!(w_vals)

# Create matrix for λmax
λ_mat = [df[(df.o .== o) .& (df.w .== w), :λmax][1] for w in w_vals, o in o_vals]

# Plot heatmap
heatmap(o_vals, w_vals, λ_mat;
        xlabel = "Omnivory Preference (o)",
        ylabel = "Habitat Preference (w)",
        title = "Max Real Eigenvalue (λmax)",
        colorbar_title = "λmax",
        c = :viridis)

##Plot C1 and C2
C1_mat = [df[(df.o .== o) .& (df.w .== w), :C1][1] for w in w_vals, o in o_vals]
C2_mat = [df[(df.o .== o) .& (df.w .== w), :C2][1] for w in w_vals, o in o_vals]

RT_mat = [df[(df.o .== o) .& (df.w .== w), :RT][1] for w in w_vals, o in o_vals]
heatmap(o_vals, w_vals, C1_mat;
    xlabel = "Omnivory Preference (o)",
    ylabel = "Habitat Preference (w)",
    title = "Equilibrium C1 Biomass",
    colorbar_title = "C1",
    c = :viridis)

# Optional second heatmap for C2
heatmap(o_vals, w_vals, C2_mat;
    xlabel = "Omnivory Preference (o)",
    ylabel = "Habitat Preference (w)",
    title = "Equilibrium C2 Biomass",
    colorbar_title = "C2",
    c = :magma)

##Return time 
heatmap(o_vals, w_vals, RT_mat;
    xlabel = "Omnivory Preference (o)",
    ylabel = "Habitat Preference (w)",
    title = "Return Time",
    colorbar_title = "C1",
    c = :viridis)

##looking at just w 
results_w = []

for w in 0.0:0.1:1.0
    λmax = NaN
    C1 = NaN
    C2 = NaN
    try
        p = ModelPar_passive(w=w)
        x0 = [1.0, 1.0, 0.5, 0.5, 0.5]

        result = nlsolve(x -> wrapped_model!(zeros(4), x, p, 0.5), x0)

        if result.f_converged
            x_star = result.zero
            R1, R2, C1, C2 = x_star
            J = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, 0.5), x_star)
            λmax = maximum(real.(eigvals(J)))
        else
            @warn "Solver did not converge for w=$w"
        end
    catch err
        @warn "Error at w=$w: $err"
    end

    push!(results_w, (w=w, λmax=λmax, C1 = C1, C2=C2))
end

df_w = DataFrame(results_w)        

plot()
plot(df_w.w, df_w.λmax)
plot(df_w.w, df_w.C1)
plot(df_w.w, df_w.C2)

##looking at just o 
results_o = []

for o in 0.0:0.1:1.0
    λmax = NaN
    try
        p = ModelPar_passive(o=o)
        x0 = [1.0, 1.0, 0.5, 0.5, 0.5]

        result = nlsolve(x -> rhs_2(x, p), x0)

        if result.f_converged
            x_star = result.zero
            J = ForwardDiff.jacobian(x -> rhs_2(x,p), x_star)
            λmax = maximum(real.(eigvals(J)))
        else
            @warn "Solver did not converge for o=$o"
        end
    catch err
        @warn "Error at o=$o: $err"
    end

    push!(results_o, (o=o, λmax=λmax))
end

df_o = DataFrame(results_o)        

plot(df_o.o, df_o.λmax)        


##the coupling plot does not look like how i would expect, why max stability at relatively strong coupling? is this related to subsidy? 
##why does omnivory always destabilize at low P values? 

##with the parameters the same across both energy channels, with equal coupling, the P basically just has more access to resources, all acting like one pool 
##essentially with this set up, all one pool ... just increasing productivity? 

##Calculate return time 


