#Setting up test equations to play around in Julia and understand language/plotting

##Install Packages
using Pkg
#Pkg.add("Plots")
#Pkg.add("Parameters")
#Pkg.add("DifferentialEquations")
#Pkg.add("Symbolics")
#Pkg.add("ForwardDiff")
#Pkg.add("LinearAlgebra")
#Pkg.add("DataFrames")
##Load Libraries needed (ODE for later)
using Parameters: @with_kw, @unpack ##imports Parameters package that provides convenient macros for working with keyword arugemnts, parameter structs and unpacking variables 
using DifferentialEquations
using Plots
using Symbolics
using ForwardDiff
using LinearAlgebra
using NLsolve
using DataFrames
#u = state variables where u[1] = R, u[2] = C, and u[3] = P
#p = Parameters
#t = 

##Population Models -----------------
##Continuous Time Logistic Population Model
#dR/dt = rR(1 - R/K) 

#Set up paramters for a simple logistic growth model for R
@with_kw mutable struct ModelPar
    #logistic Parameters
    r = 2.0
    K = 100.0
end

##Set up model function
function model!(du, u, p, t)
    @unpack r, K = p
    R = u[1]

    #ODE
    du[1]= r * R * (1 - R / K)
    return du
end

#initial condition 
u0 = [10.0]
tspan = (0.0, 100.0)

#parameter
params = ModelPar()

#Solve the equation
prob = ODEProblem(model!, u0 , tspan, params)
sol = solve(prob)

population = [u[1] for u in sol.u]
plot(sol.t, population, label = "R(t)", xlabel = "Time", ylabel = "Population", title = "Logistic Growth Model", lw = 2)


##Consumer-Resource Models ------------
##Type 1 functional response 
#dR/dt = rR(1-r/k) - aRC
#dC/dt = eaRC - mC 

#Set up paramters 
@with_kw mutable struct ModelPar_2
    r = 0.5
    K = 1.0
    a = 0.910
    e = 1.0
    m =0.4
    a2 = 2.75
    m2 = 0.3
    b2 = 0.5
    b = 0.5
end

##
function model!(du, u, p, t)
    @unpack r, K, a, e, m = p 
    R, C = u 

    #ODE
    du[1] = r * R * (1 - R / K) - a * R * C 
    du[2] = e * a * R * C - m * C 
end

# Utilities for doing eigenvalue analysis
##this is a right-hand-side function that is used to solve ODEs 
function rhs(u, p, t=0.0)
    du = similar(u) ##creates a new array du with the same type and size as u but uninitialized
    model!(du, u, p, t) 
    return du ##this returns the computed derivative (or result of model) stored in du 
end

##set initial condition
u0 = [0.66, 0.27]
tspan = (0.0, 1000.0)

##set Parameters
p = ModelPar_2()


##Define the ODE problem
prob2 = ODEProblem(rhs, u0, tspan, p)
sol2 = solve(prob2)

##plot timeseries
plot(sol2, xlabel="Time", ylabel="Population", title="ODE Solution")

##calculate numerical jacobian
J_num = ForwardDiff.jacobian(x -> rhs(x, p), u0)
println("Jacobian Matrix:")
display(J_num)
##find the equilibrium for this system of equations 
eq = nlsolve((du, u) -> model!(du, u, p, 0.0), u0).zero

# Display equilibrium
println("Equilibrium point (R*, C*): ", eq)
##this now makes sense w/ plot 

# Create phase plane grid
R_range = 0:0.1:1.0
C_range = 0:0.1:1.0
#R_grid, C_grid = [R for R in R_range], [C for C in C_range]

# Prepare vector field for quiver
R_grid = repeat(R_range', length(C_range), 1)
C_grid = repeat(C_range, 1, length(R_range))

dR = zeros(size(R_grid))
dC = zeros(size(C_grid))

for i in 1:size(R_grid, 1)
    for j in 1:size(R_grid, 2)
        u = [R_grid[i, j], C_grid[i, j]]
        du = similar(u)
        model!(du, u, p, 0.0)
        norm_du = norm(du) + 1e-8  # Avoid zero division
        dR[i, j] = du[1] / norm_du
        dC[i, j] = du[2] / norm_du
    end
end

# Plot phase portrait
plot(title="Phase Portrait with Isoclines", xlabel="Prey (R)", ylabel="Predator (C)", legend=true)
#quiver!(R_grid, C_grid, quiver=(dR, dC), c=:blue, aspect_ratio=1, label="Vector Field")

##extract parameters from parameter structure to be able to calculate isoclines
r, K, a, e, m = p.r, p.K, p.a, p.e, p.m
# Plot R-isocline (dR/dt = 0)
C_nullcline_R = (r / a) .* (1 .- R_range ./ K)
plot!(R_range, C_nullcline_R, lw=2, lc=:red, label="R-isocline (dR/dt=0)")

# Plot C-isocline (dC/dt = 0)
R_nullcline_C = fill(m / (e * a), length(C_range))
plot!(R_nullcline_C, C_range, lw=2, lc=:green, label="C-isocline (dC/dt=0)")

# Add trajectories for different initial conditions
initial_conditions = [[0.5, 0.5], [1.5, 0.5], [1.0, 1.5], [0.8, 0.8]]

for u0 in initial_conditions
    prob = ODEProblem(model!, u0, tspan, p)
    sol = solve(prob)
    plot!(sol[1, :], sol[2, :], lw=2, label="")
end

# Show plot
plot!()



##Determine maximum real eigenvalues - from numerical jacobian
max_eig = maximum(real.(eigvals(J_num)))

##calculate C:R ratio 
# Function to compute C:R ratio analytically
function CR_ratio(p)
    @unpack r, K, a, e, m = p 

    ratio = (r * e) / m * (1 - m / (e * a * K))
    return ratio
end

# Compute the C:R ratio
cr = CR_ratio(p)
println("C:R Ratio (Analytical) = ", cr)

##Look at C:R ratio over range of parameters a,e,m 
##range of attack rate
a_range = range(0.5, stop = 1.0, length = 50)

##initialize a df to store results
df_a = DataFrame(a=Float64[], CR_ratio=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for a in a_range
    p2 = (r =r, K=K, a=a, e=e, m=m)
    cr = CR_ratio(p2)
    #store in df
    push!(df_a, (a,cr))
end
println(df_a)
plot(df_a.a, df_a.CR_ratio, xlabel = "Attack Rate (a)", ylabel = "C:R Ratio", title = "C:R Ratio vs Attack Rate", legend = false, linewidth = 2)

##sweet this looks exactly like mathematica tutorial 

##range of m 
m_range = range(0.1, stop = 1.0, length = 50)

##initialize a df to store results
df_m = DataFrame(m=Float64[], CR_ratio=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for m in m_range
    p3 = (r =r, K=K, a=a, e=e, m=m)
    cr = CR_ratio(p3)
    #store in df
    push!(df_m, (m,cr))
end
println(df_m)
plot(df_m.m, df_m.CR_ratio, xlabel = "Mortality (m)", ylabel = "C:R Ratio", title = "C:R Ratio vs Mortality", legend = false, linewidth = 2)


##range of conversion efficiency e
e_range = range(0.4, stop = 1.0, length = 50)

##initialize a df to store results
df_e = DataFrame(e=Float64[], CR_ratio=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for e in e_range
    p4 = (r =r, K=K, a=a, e=e, m=m)
    cr = CR_ratio(p4)
    #store in df
    push!(df_e, (e,cr))
end
println(df_e)
plot(df_e.e, df_e.CR_ratio, xlabel = "Conversion Efficienct (e)", ylabel = "C:R Ratio", title = "C:R Ratio vs Conversion Efficiency", legend = false, linewidth = 2)

##range of R's growth rate (r)
r_range = range(0.1, stop = 1.0, length = 50)

##initialize a df to store results
df_r = DataFrame(r=Float64[], CR_ratio=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for r in r_range
    p5 = (r =r, K=K, a=a, e=e, m=m)
    cr = CR_ratio(p5)
    #store in df
    push!(df_r, (r,cr))
end
println(df_r)
plot(df_r.r, df_r.CR_ratio, xlabel = "Growth Rate (r)", ylabel = "C:R Ratio", title = "C:R Ratio vs Growth Rate", legend = false, linewidth = 2)

##Range of K 
K_range = range(0.5, stop = 2.0, length = 50)

##initialize a df to store results
df_K = DataFrame(K=Float64[], CR_ratio=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for K in K_range
    p6 = (r =r, K=K, a=a, e=e, m=m)
    cr = CR_ratio(p6)
    #store in df
    push!(df_K, (K,cr))
end
println(df_K)
plot(df_K.K, df_K.CR_ratio, xlabel = "Growth Rate (r)", ylabel = "C:R Ratio", title = "C:R Ratio vs Growth Rate", legend = false, linewidth = 2)

##To-do: need to work through this and make sure the math actually makes sense and that it 
##looks the same as mathematica code, then also figure out how to plot C:R ratio, and eigenvalues over parameter space

###Looking at eigenvalues across parameter ranges
a_range_eig = range(0.3, stop = 1.0, length = 100)

##initialize a df to store results
df_a_eig = DataFrame(a=Float64[], max_eig=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for a in a_range_eig
    p_7 = (r =r, K=K, a=a, e=e, m=m)
    #compute equilibrium point for each parameter set analytically
    R_star = m / (e * a)
    C_star = (r /a) * (1-R_star / K)

    ##check if equilibrium is biologically feasible (i.e., >0)
    if R_star > 0 && C_star > 0
        u_eq = [R_star, C_star] ##generate the equilibrium point
        J_num1 = ForwardDiff.jacobian(x -> rhs(x, p_7), u_eq)
        max_eig = maximum(real.(eigvals(J_num1)))
         #store in df
         push!(df_a_eig, (a,max_eig))
    else 
        push!(df_a_eig, (a, NaN)) ##add NA if equilibrium point is not biologically feasible   
    end  
end
println(df_a_eig)
plot(df_a_eig.a, df_a_eig.max_eig, xlabel = "Attack Rate (a)", ylabel = "Maximum eigenvalue", title = "Max eigenvalue vs Attack Rate", legend = false, linewidth = 2)
hline!([0], linestyle=:dash, color=:red, label="Stability Threshold")

##Range of e 
e_range_eig = range(0.3, stop = 1.0, length = 100)

##initialize a df to store results
df_e_eig = DataFrame(e=Float64[], max_eig=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for e in e_range_eig
    p_8 = (r =r, K=K, a=a, e=e, m=m)
    #compute equilibrium point for each parameter set analytically
    R_star = m / (e * a)
    C_star = (r /a) * (1-R_star / K)

    ##check if equilibrium is biologically feasible (i.e., >0)
    if R_star > 0 && C_star > 0
        u_eq = [R_star, C_star] ##generate the equilibrium point
        J_num1 = ForwardDiff.jacobian(x -> rhs(x, p_8), u_eq)
        max_eig = maximum(real.(eigvals(J_num1)))
         #store in df
         push!(df_e_eig, (e,max_eig))
    else 
        push!(df_e_eig, (e, NaN)) ##add NA if equilibrium point is not biologically feasible   
    end  
end
println(df_e_eig)
plot(df_e_eig.e, df_e_eig.max_eig, xlabel = "e", ylabel = "Maximum eigenvalue", title = "Max eigenvalue vs e", legend = false, linewidth = 2)
hline!([0], linestyle=:dash, color=:red, label="Stability Threshold")


##range of K 
K_range_eig = range(0.5, stop = 2.0, length = 100)

##initialize a df to store results
df_K_eig = DataFrame(K=Float64[], max_eig=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for K in K_range_eig
    p_9 = (r =r, K=K, a=a, e=e, m=m)
    #compute equilibrium point for each parameter set analytically
    R_star = m / (e * a)
    C_star = (r /a) * (1-R_star / K)

    ##check if equilibrium is biologically feasible (i.e., >0)
    if R_star > 0 && C_star > 0
        u_eq = [R_star, C_star] ##generate the equilibrium point
        J_num1 = ForwardDiff.jacobian(x -> rhs(x, p_9), u_eq)
        max_eig = maximum(real.(eigvals(J_num1)))
         #store in df
         push!(df_K_eig, (K,max_eig))
    else 
        push!(df_K_eig, (K, NaN)) ##add NA if equilibrium point is not biologically feasible   
    end  
end
println(df_K_eig)
plot(df_K_eig.K, df_K_eig.max_eig, xlabel = "K", ylabel = "Maximum eigenvalue", title = "Max eigenvalue vs K", legend = false, linewidth = 2)
hline!([0], linestyle=:dash, color=:red, label="Stability Threshold")


##range of m
m_range_eig = range(0.5, stop = 2.0, length = 100)

##initialize a df to store results
df_m_eig = DataFrame(m=Float64[], max_eig=Float64[])
##Loop over values of a to compute C:R ratio 
#r, K, a, e, m = p.r, p.K, p.a, p.e, p.m

for m in m_range_eig
    p_10 = (r =r, K=K, a=a, e=e, m=m)
    #compute equilibrium point for each parameter set analytically
    R_star = m / (e * a)
    C_star = (r /a) * (1-R_star / K)

    ##check if equilibrium is biologically feasible (i.e., >0)
    if R_star > 0 && C_star > 0
        u_eq = [R_star, C_star] ##generate the equilibrium point
        J_num1 = ForwardDiff.jacobian(x -> rhs(x, p_10), u_eq)
        max_eig = maximum(real.(eigvals(J_num1)))
         #store in df
         push!(df_m_eig, (m,max_eig))
    else 
        push!(df_m_eig, (m, NaN)) ##add NA if equilibrium point is not biologically feasible   
    end  
end
println(df_m_eig)
plot(df_m_eig.m, df_m_eig.max_eig, xlabel = "m", ylabel = "Maximum eigenvalue", title = "Max eigenvalue vs m", legend = false, linewidth = 2)
hline!([0], linestyle=:dash, color=:red, label="Stability Threshold")




##so far yes looks the same as mathematica  code, want to create a new repo that has separate scripts for these different model types etc. 



##finding symbolic jacobian matrix
#J = Symbolics.jacobian(x -> rhs(x, p, t), u)

##testing out getting symbolic jacobian matrix 
# Define symbolic variables for states and parameters
@variables R C r K a e m a2
state_vars = [R, C]
param_vars = (r=r, K=K, a=a, e=e, m=m, a2=a2)

# Create a wrapper for the model! function to handle symbolic inputs
function symbolic_model(u, p)
    du = similar(u)
    model!(du, u, p, 0.0)  # Time is not needed for autonomous systems
    return du
end


# Evaluate the model symbolically
du = symbolic_model(state_vars, param_vars)

# Compute the symbolic Jacobian
J = Symbolics.jacobian(du, state_vars)

# Display the symbolic Jacobian
println("Symbolic Jacobian Matrix:")
display(J)
##so this jacobian looks right - same as in mathematica code 



##helper functions 
find_eq(u, p) = nlsolve((du, u) -> model!(du, u, p, zero(u)), u).zero
cmat(u, p) = ForwardDiff.jacobian(x -> rhs(x, p), u)


"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
λ1_stability(M) = maximum(real.(eigvals(M)))

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
function λ1_stability_imag(M) 
    ev = eigvals(M)
    imag.(ev)[findall(real.(ev) .== maximum(real.(ev)))[1]]
end 


"""M is the community matrix, we can be calculated with `cmat(u, p)`
Note: `\nu` is the what to input `ν` which looks a bit too much like `v` for my taste
"""
ν_stability(M) = λ1_stability((M + M') / 2)