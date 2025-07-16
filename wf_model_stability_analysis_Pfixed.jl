##Looking at stability across gradient of coupling and omnivory 

##source model 
include("wf_model_eqs_subsidy_2.jl") ####model equations with P held constant, unique parameters per trophic level, active and passive omnivory parameter structures, temporal forcing in R
##MAKE SURE PROPER P IS SILENCED IN CODE BEFORE RUNNING HERE 

##Calculating local equilibrium and stability (eigenvalues)
# Initial guess for dynamic variables: R1, R2, C1, C2
x0 = [1.0, 1.0, 0.5, 0.5]
P_fixed = 0.2
p = ModelPar2()

##Local stability analysis
#1) Calculate equilibrium values for system of equation 
result = nlsolve(x -> wrapped_model_vec(x, p, P_fixed), x0)

# Extract equilibrium values
x_star = result.zero
R1_star, R2_star, C1_star, C2_star = x_star
 

#2) Calculate numerical jacobian
J_num = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, P_fixed), x0)
println("Jacobian Matrix:")
display(J_num)

#3) Calculate the eigenvalues 
eigvals(J_num)


##Trying to find multiple equilibriums, as likely more than 1 
#Pkg.add("IterTools")
using IterTools

R_vals = [0.1, 0.5, 1.0, 2.0]
C_vals = [0.1, 0.5, 1.0]

# All combinations of R1, R2, C1, C2
guesses = collect(product(R_vals, R_vals, C_vals, C_vals))

P_fixed = 0.2
tolerance = 1e-4
equilibria_1 = []

for x0 in guesses
    result = nlsolve(x -> wrapped_model_vec(x, p, P_fixed), collect(x0))
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
    J = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, P_fixed), eq)
    λ = eigvals(J)
    λ_real = real.(λ)
    λ_max = maximum(λ_real)

    is_stable = all(λ_real .< 0)
    stability = is_stable ? "Stable" : "Unstable or saddle"
    
    println("\nEquilibrium $i:")
    println("  State: ", round.(eq; digits=4))
    println("  Stability: ", stability)
    println("  Max real eigenvalue: ", round(λ_max; digits=6))
    println("  All eigenvalues: ", round.(λ; digits=4))

    push!(max_eigenvalues_1, λ_max)
end




# Grid of values for o and w

results = []

for o in 0.0:0.1:1.0, w in 0.0:0.1:1.0
    λmax = NaN
    C1 = NaN
    C2 = NaN
    RT = NaN
    try
        p = ModelPar2(o=o, w=w)
        x0 = [1.0, 1.0, 0.5, 0.5]

        result = nlsolve(x -> wrapped_model!(zeros(4), x, p, 0.2), x0)

        if result.f_converged
            x_star = result.zero
            R1, R2, C1, C2 = x_star
            J = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, 0.2), x_star)
            λmax = maximum(real.(eigvals(J)))
            RT = 1/abs(λmax)
        else
            @warn "Solver did not converge for o=$o, w=$w"
        end
    catch err
        @warn "Error at o=$o, w=$w: $err"
    end

    push!(results, (o=o, w=w, λmax=λmax, RT=RT,C1=C1, C2=C2))
end

df = DataFrame(results)

##some issues with convergence, need to dig deeper to figure out what is going on 

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
        x0 = [1.0, 1.0, 0.5, 0.5]

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

for o in 0.0:0.01:1.0
    λmax = NaN
    try
        p = ModelPar2(o=o)
        x0 = [1.0, 1.0, 0.5, 0.5]

        result = nlsolve(x -> wrapped_model!(zeros(4), x, p, 0.5), x0)

        if result.f_converged
            x_star = result.zero
            J = ForwardDiff.jacobian(x -> wrapped_model_vec(x, p, 0.5), x_star)
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
