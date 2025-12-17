##Code to look at try out a type of sensitivty analysis across parameter ranges 
##Date Initiated: October 10, 2025
##Contributor(s): Marie K. Gutgesell

#using Pkg
#Pkg.add("GlobalSensitivity")
#Pkg.add("QuasiMonteCarlo")
#Pkg.add("StatsPlots")
#Pkg.add("DataFrames")
#Pkg.add("Parquet")
#using DifferentialEquations, ForwardDiff, LinearAlgebra
#using NLsolve
using GlobalSensitivity
using QuasiMonteCarlo   # for LHS/Sobol sampling
#using Statistics
using Random
using StatsPlots
using Base.Threads
using DataFrames
using Parquet
using Plots
using Statistics
##source model - choose which based on which you want to investigate
#project_root = raw"C:\Users\mccan\Documents\Github\Gutgesell"
include( "wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P held constant, unique parameters per trophic level, active and passive omnivory parameter structures

## base parameter set to copy/override from
const p0 = ModelPar_active()

##1) Sensitivity analysis just for structure -- so keeping H at 0 (constant)
# --- Bounds: focal vs nuisance ---
focal_syms    = (:o, :w)   # for this first one only focusing on o and w 
nuisance_syms = (:H, :r, :K, :aR_C, :aR_P, :aC_P, :aG_P, :hR_P, :hR_C, :hC_P, :hG_P, :e, :mC, :G)

# Ranges (based on univariate stability analysis)
bounds = Dict(
    :o => (0.0, 1.0),
    :w => (0.0, 1.0),
    :H => (0.0, 0.5), ## keeping H at 0 for unforced structure runs
    :r => (0.2, 3.0),
    :K => (1.0, 3.8),
    :aR_P => (1.0, 2.5),
    :aR_C => (1.0, 3.0),
    :aC_P => (0.8, 1.8),
    :aG_P => (1.0, 3.0),
    :hR_C => (0.4, 0.8),
    :hR_P => (0.5, 2.0),
    :hC_P => (1.0, 2.0),
    :hG_P => (0.8, 2.0),
    :e   => (0.7, 1.0),
    :mC  => (0.6, 1.2),
   # :mP  => (0.01, 1.5),
    :G => (1.0, 5.0),
  #  :G_base => (0.0, 10.0),
 #   :l  => (0.0, 1.0),
 #   :pf => (0.5, 10.0),
 #   :D => (0.0, 1.0)
)

# Ranges (bounds -- based on 20% around values used in MS)
bounds = Dict(
    :o => (0.0, 1.0),
    :w => (0.0, 1.0),
    :H => (0.0, 1.0), ##is this a way to make sure it is always 0 ? 
    :r => (0.8, 1.2), ##1.0
    :K => (2.44, 3.66), ##3.05
    :aR_P => (3.2, 4.8), ##4.0
    :aR_C => (2.0, 3.0), ##2.5
    :aC_P => (2.72, 4.08), #3.4
    :aG_P => (2.72, 4.08), #3.4
    :hR_C => (0.32, 0.48), #0.4
    :hR_P => (1.0, 1.5), #1.25
    :hC_P => (1.0, 1.5), #1.25#
    :hG_P => (1.0, 1.5), #1.25
    :e   => (0.64, 0.96), #0.8
    :mC  => (0.8, 1.2), #1.0
   # :mP  => (0.01, 1.5),
    :G => (1.6, 2.4), #2.0
  #  :G_base => (0.0, 10.0),
 #   :l  => (0.0, 1.0),
 #   :pf => (0.5, 10.0),
 #   :D => (0.0, 1.0)
)
##can also try just a 10-20% variance around the values we have in our model 
##could also try drawing from normal distribution rather than from a uniform 
# Optional: which are log-scaled? - good for ones that span orders of magnitude
#logscale = Set([:aR_P, :aR_C, :aG_P, :aC_P, :mC, :hR_C, :hR_P, :hC_P, :hG_P])

##Helper functions - 
# Map unit cube sample x∈[0,1] to parameter in [lo,hi] (linear or log)
map_to_range(x, lo, hi; logscaled=false) =
    logscaled ? exp(log(lo) + x*(log(hi) - log(lo))) : (lo + x*(hi - lo))


# Draw N samples of nuisance params with LHS (unconstrained)
function sample_nuisance(N::Int; rng=Random.default_rng())
    d = length(nuisance_syms)
    X = QuasiMonteCarlo.sample(d, N, LatinHypercubeSample())   # N × d in [0,1]
    samples = Vector{NamedTuple}(undef, N)
    for i in 1:N
        pairs = ntuple(j -> begin
            s = nuisance_syms[j]
            lo, hi = bounds[s]
        #    val = map_to_range(X[i,j], lo, hi; logscaled = (s in logscale))
         val = map_to_range(X[i,j], lo, hi) ##trying out without logscaling 
            (s => val)
        end, d)
        samples[i] = NamedTuple(pairs)
    end
    samples
end

## ===============================
## NEW: parameter constructors
## ===============================
# Unforced version: H taken from ν.H (which is fixed at 0 in current bounds)
function make_par(o::Real, w::Real, ν::NamedTuple)
   ModelPar_active(p0;
      o = o, w = w, H = ν.H,
       r = ν.r, K = ν.K,
      aR_P = ν.aR_P, aC_P = ν.aC_P, aG_P = ν.aG_P, aR_C = ν.aR_C,
        hR_P = ν.hR_P, hC_P = ν.hC_P, hG_P = ν.hG_P, hR_C = ν.hR_C,
        e = ν.e, mC = ν.mC, 
        G = ν.G
    )
end

# General version: H is a focal parameter, not taken from ν
function make_par_forced(o::Real, w::Real, H::Real, ν::NamedTuple)
    ModelPar_active(p0;
        o = o, w = w, H = H,
        r = ν.r, K = ν.K,
        aR_P = ν.aR_P, aC_P = ν.aC_P, aG_P = ν.aG_P, aR_C = ν.aR_C,
        hR_P = ν.hR_P, hC_P = ν.hC_P, hG_P = ν.hG_P, hR_C = ν.hR_C,
        e = ν.e, mC = ν.mC, 
        G = ν.G
    )
end

# ## ============================================
# ## NEW: constrained sampling PER (o,w) grid cell
# ## ============================================

# """
#     sample_feasible_nuisance_for_grid(o, w, N; max_tries=50_000)

# Return a vector of N nuisance parameter NamedTuples that are feasible
# for this (o,w) point in the **unforced** system.

# Feasibility is checked via `is_feasible_paramset`.
# """
# function sample_feasible_nuisance_for_grid(o::Real, w::Real, N::Int;
#                                            rng=Random.default_rng(),
#                                            max_tries::Int = 50_000)

#     feasible = NamedTuple[]
#     tries = 0
#     d = length(nuisance_syms)

#     while length(feasible) < N && tries < max_tries
#         batch = max(N, 50)
#         X = QuasiMonteCarlo.sample(d, batch, LatinHypercubeSample())

#         for i in 1:batch
#             pairs = ntuple(k -> begin
#                 s = nuisance_syms[k]
#                 lo, hi = bounds[s]
#                 val = map_to_range(X[i,k], lo, hi; logscaled = (s in logscale))
#                 (s => val)
#             end, d)
#             ν = NamedTuple(pairs)

#             # build params for unforced system at this (o,w)
#             p = make_par(o, w, ν)

#             # Each candidate counts as one "attempt"
#             tries += 1

#             if is_feasible_paramset(p)
#                 push!(feasible, ν)
#                 length(feasible) == N && break
#             end

#             # stop if we hit max_tries
#             tries >= max_tries && break
#         end
#     end

#     if length(feasible) < N
#         @warn "Grid point (o=$(round(o,digits=3)), w=$(round(w,digits=3))) only reached $(length(feasible)) feasible samples after $tries draws."
#     end

#     return feasible, tries
# end

# ## old constrained sampler (now just a thin wrapper around the new one
# ## using an anchor point, in case you still want it for other tasks)
# function sample_nuisance_constrained(N::Int; rng=Random.default_rng(), max_tries = 10000)
#     sample_feasible_nuisance_for_grid(0.5, 0.5, N; rng=rng, max_tries=max_tries)
# end

# ## ==============================
# ## Grids for the focal parameters
# ## ==============================

# # Grids for the 3 focal parameters
# o_grid = range(bounds[:o]...; length=11)
# w_grid = range(bounds[:w]...; length=11)
# # for now, keep H=0 for the unforced structure runs
# H_grid = range(0.0, 0.0; length=1)  # adjust if you want a true H range later

# Nrep = 3  # random nuisance samples per grid point (tune)

# No, Nw = length(o_grid), length(w_grid)

# ## ============================================
# ## NEW: Feasible nuisance pools for each (o,w)
# ## ============================================
# @info "Generating feasible nuisance pools for each (o,w) grid point..."

# feasible_pool_grid = [Vector{NamedTuple}() for _ in 1:No, _ in 1:Nw]
# attempts_grid      = Array{Int}(undef, No, Nw)   # NEW: attempts per (o,w)

# @threads for idx in 1:(No * Nw)
#     io = fld(idx-1, Nw) + 1
#     iw = mod(idx-1, Nw) + 1
#     o = o_grid[io]
#     w = w_grid[iw]

#     pool, tries = sample_feasible_nuisance_for_grid(o, w, Nrep)
#     feasible_pool_grid[io, iw] = pool
#     attempts_grid[io, iw]      = tries
# end

# @info "Done generating feasible pools."

# feasible_pool_grid
# attempts_grid

# using Plots

# heatmap(o_grid, w_grid, attempts_grid';
#         xlabel = "o", ylabel = "w", colorbar_title = "Attempts",
#         title = "Number of attempts to find Nrep feasible samples")




#############################################
### DEBUG: Test feasibility at one (o,w) ####
#############################################

# function test_one_grid_point(o_test, w_test;
#                              N_test = 50,
#                              P0 = 0.25,
#                              rng = Random.default_rng())

#     println("\n=== Testing grid point (o=$(o_test), w=$(w_test)) ===")
#     println("Sampling $N_test nuisance sets...\n")

#     # draw nuisance samples (unconstrained)
#     ν_list = sample_nuisance(N_test; rng=rng)

#     results = NamedTuple[]

#     for (i, ν) in enumerate(ν_list)

#         # construct p
#         p = make_par(o_test, w_test, ν)

#         # evaluate equilibrium
#         eq, λ1, reason = nothing, nothing, ""

#         out = try
#             r = equilibrium_unforced(p, P0)
#             eq = r.eq
#             λ1 = r.λ1
#             reason = "OK"
#             r
#         catch e
#             eq = fill(NaN, 4)   # placeholder
#             λ1 = NaN
#             reason = "ERROR: $(e)"
#             nothing
#         end

#         feasible = all(isfinite.(eq)) && all(eq .> 0)

#         println("--------------------------------------------------")
#         println("Sample $i")
#         println("ν = $ν")
#         println("Equilibrium eq = $eq")
#         println("Feasible? $feasible   |  reason: $reason")
#         println("--------------------------------------------------\n")

#         push!(results, (
#             ν = ν,
#             eq = eq,
#             feasible = feasible,
#             λ1 = λ1,
#             reason = reason
#         ))
#     end

#     return results
# end

# debug_results = test_one_grid_point(1.0, 1.0; N_test = 50)


function test_one_grid_point_collect_feasible(o_test, w_test;
                                              Nrep = 3,          # number of feasible samples desired
                                              P0 = 0.25,
                                              tol = 10e-3,
                                              max_tries = 50_000,
                                              rng = Random.default_rng())

    println("\n=== Testing grid point (o=$(o_test), w=$(w_test)) ===")
    println("Seeking $Nrep feasible nuisance samples...\n")

    feasible_results = NamedTuple[]
    tries = 0

    while length(feasible_results) < Nrep && tries < max_tries

        # --- draw 1 nuisance set ---
        ν = sample_nuisance(1; rng=rng)[1]
        tries += 1

        # --- construct parameter set ---
        p = make_par(o_test, w_test, ν)

        # --- evaluate equilibrium ---
        eq = nothing
        λ1 = nothing
        reason = ""

        out = try
            r = equilibrium_unforced(p, P0)
            eq = r.eq
            λ1 = r.λ1
            reason = "OK"
            r
        catch e
            eq = fill(NaN, 4)
            λ1 = NaN
            reason = "ERROR: $(e)"
            nothing
        end

        feasible = all(isfinite.(eq)) && all(eq .> tol)

        # --- print diagnostic info ---
        println("--------------------------------------------------")
        println("Attempt $tries")
        println("ν = $ν")
        println("Equilibrium eq = $eq")
        println("Feasible? $feasible   |  reason: $reason")
        println("--------------------------------------------------\n")

        # --- keep only the feasible ones ---
        if feasible
            push!(feasible_results, (
                ν = ν,
                eq = eq,
                λ1 = λ1,
                reason = reason
            ))
        end
    end

    # --- warn if we failed to find enough ---
    if length(feasible_results) < Nrep
        @warn "Only found $(length(feasible_results)) feasible samples after $tries attempts."
    else
        println("SUCCESS: Found $Nrep feasible samples after $tries attempts.")
    end

    return feasible_results
end


feasibles = test_one_grid_point_collect_feasible(1.0, 1.0; Nrep=3)




function test_one_grid_point_collect_feasible(o_test, w_test;
                                              Nrep = 3,          # number of feasible samples desired
                                              P0 = 0.25,
                                              tol = 10e-3,
                                              max_tries = 50_000,
                                              rng = Random.default_rng())

    println("\n=== Testing grid point (o=$(o_test), w=$(w_test)) ===")
    println("Seeking $Nrep feasible nuisance samples...\n")

    feasible_results = NamedTuple[]
    tries = 0

    while length(feasible_results) < Nrep && tries < max_tries

        # draw 1 nuisance set
        ν = sample_nuisance(1; rng=rng)[1]
        tries += 1

        # construct parameter set
        p = make_par(o_test, w_test, ν)

        # evaluate equilibrium
        eq = nothing
        λ1 = nothing
        reason = ""

        out = try
            r = equilibrium_unforced(p, P0)
            eq = r.eq
            λ1 = r.λ1
            reason = "OK"
            r
        catch e
            eq = fill(NaN, 4)
            λ1 = NaN
            reason = "ERROR: $(e)"
            nothing
        end

        feasible = all(isfinite.(eq)) && all(eq .> tol)

        println("--------------------------------------------------")
        println("Attempt $tries")
        println("ν = $ν")
        println("Equilibrium eq = $eq")
        println("Feasible? $feasible   |  reason: $reason")
        println("--------------------------------------------------\n")

        if feasible
            push!(feasible_results, (
                ν = ν,
                eq = eq,
                λ1 = λ1,
                reason = reason
            ))
        end
    end

    if length(feasible_results) < Nrep
        @warn "Only found $(length(feasible_results)) feasible samples after $tries attempts."
    else
        println("SUCCESS: Found $Nrep feasible samples after $tries attempts.")
    end

    # return both the feasible samples and how many attempts it took
    return (results = feasible_results, tries = tries)
end

function build_feasible_pools_grid(o_grid, w_grid;
                                   Nrep = 3,
                                   P0 = 0.25,
                                   tol = 10e-3,
                                   max_tries = 50_000,
                                   rng = Random.default_rng())

    No, Nw = length(o_grid), length(w_grid)

    # Each cell holds a Vector{NamedTuple} of ν for that (o,w)
    feasible_pool_grid = [Vector{NamedTuple}() for _ in 1:No, _ in 1:Nw]
    # Attempts taken at each (o,w)
    attempts_grid      = Array{Int}(undef, No, Nw)

    for io in 1:No
        for iw in 1:Nw
            o = o_grid[io]
            w = w_grid[iw]

            out = test_one_grid_point_collect_feasible(o, w;
                                                       Nrep      = Nrep,
                                                       P0        = P0,
                                                       tol       = tol,
                                                       max_tries = max_tries,
                                                       rng       = rng)

            results = out.results
            tries   = out.tries

            # store only the ν sets in feasible_pool_grid, to match your earlier design
            feasible_pool_grid[io, iw] = [r.ν for r in results]
            attempts_grid[io, iw]      = tries
        end
    end

    return feasible_pool_grid, attempts_grid
end

o_grid = range(bounds[:o]...; length=11)
w_grid = range(bounds[:w]...; length=11)
Nrep   = 3

feasible_pool_grid, attempts_grid =
    build_feasible_pools_grid(o_grid, w_grid; Nrep=Nrep, P0=0.25, tol = 10e-3, max_tries=50_000)

feasible_pool_grid
attempts_grid

# ## =================================================================================
# ## Build library of all parameter sets to test each o-w combination on 
# ## =================================================================================

param_library = [
    (; o = o_grid[i], w = w_grid[j], params = p)
    for i in axes(feasible_pool_grid, 1)
    for j in axes(feasible_pool_grid, 2)
    for p in feasible_pool_grid[i, j]
]

# ## =================================================================================
# ## Unforced robustness grid: use the feasible_pool_grid for each (o,w) in parallel
# ## =================================================================================
No, Nw = length(o_grid), length(w_grid)
# # Allocate a cell array that stores all runs for each (o,w) -- for unforced model 
 runs_stab_unforced = [NamedTuple[] for _ in 1:No, _ in 1:Nw]

 @info "Running robustness grid for UNFORCED model..."
 @threads for idx in 1:(No * Nw)
     io = fld(idx-1, Nw) + 1
     iw = mod(idx-1, Nw) + 1
     o = o_grid[io]
     w = w_grid[iw]

     P0 = 0.25
     u0 = [1.5, 1.5, 1.0, 1.0, 0.25]
     local_runs = NamedTuple[]

     #local_pool = feasible_pool_grid[io, iw]

     for s in param_library
         try
             ν = s.params
            p = make_par(o, w, ν)
             out_1 = equilibrium_unforced(p, P0)
             out_2 = fr_cv_unforced(p; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)

             rec = (; o, w,
               o_src = s.o, w_src = s.w,
                ν...,
                     C1_min = out_1.min[3],
                     C2_min = out_1.min[4],
                     R1_min = out_1.min[1],
                     R2_min = out_1.min[2],
                     λ1     = out_1.λ1,
                     λ1_imag     = out_1.λ1_imag,
                     R1_eq  = out_1.eq[1],
                     R2_eq  = out_1.eq[2],
                     C1_eq  = out_1.eq[3],
                     C2_eq  = out_1.eq[4],
                     cv_total = out_2.cv_total,
                     cv_R1    = out_2.cv_R1,
                     cv_R2    = out_2.cv_R2,
                     cv_C1    = out_2.cv_C1,
                     cv_C2    = out_2.cv_C2,
                     cv_G     = out_2.cv_G)

             push!(local_runs, rec)
         catch err
             @warn "Unforced fail at (o=$o, w=$w): $err"
         end
     end

     runs_stab_unforced[io, iw] = local_runs
 end
 @info "Done UNFORCED robustness grid."

 runs_stab_unforced
 cell = runs_stab_unforced[1,1]
 cell[2]

##flatten results matrix to a DataFrame
all_runs_unforced = reduce(vcat, vec(runs_stab_unforced))
all_runs_unforced = [(; run_id = i, r...) for (i, r) in enumerate(all_runs_unforced)]
all_runs_unforced_df = DataFrame(all_runs_unforced)
write_parquet("runs_stab_unforced.parquet", all_runs_unforced_df)


all_runs_unforced_df

λ = all_runs_unforced_df[:, :λ1]
λ = filter(isfinite, λ)


bw = 0.05
edges = minimum(λ):bw:10

histogram(λ;
    bins = edges,
    xlabel = "λ₁ (max real eigenvalue)",
    ylabel = "Count"
)


##okay sick, so now, want to calculate proportion of feasible, and then median CV
tol = 1e-3

all_runs_unforced_df.feasible = (
    (all_runs_unforced_df.R1_eq .> tol) .&
    (all_runs_unforced_df.R2_eq .> tol) .&
    (all_runs_unforced_df.C1_eq .> tol) .&
    (all_runs_unforced_df.C2_eq .> tol)
)
summary_ow = combine(
    groupby(all_runs_unforced_df, [:o, :w]),
    :feasible => mean => :prop_feasible,
    :feasible => sum  => :n_feasible,
    :cv_total => length => :n_total,
    [:cv_total, :feasible] =>
        ((cv, f) -> median(cv[f])) => :median_cv_feasible,
  [:λ1, :feasible] =>
        ((λ, f) -> any(f) ? median(λ[f]) : missing) => :median_λ1_feasible
)


##transpose back into a grid for Plotting
prop_grid = unstack(summary_ow, :w, :o, :prop_feasible)
median_cv_grid = unstack(summary_ow, :w, :o, :median_cv_feasible)
eig_grid = unstack(summary_ow, :w, :o, :median_λ1_feasible)
# y-axis (rows)
w_vals = median_cv_grid.w

# x-axis (column names, skipping :w)
o_syms = names(median_cv_grid)[2:end]
o_vals = parse.(Float64, string.(o_syms))

# z matrix
Z = Matrix(median_cv_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z;
    xlabel = "o",
    ylabel = "w",
    title  = "Median CV (feasible runs)",
    colorbar_title = "median CV",
    aspect_ratio = :equal
)


# z matrix
Z_2= Matrix(prop_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z_2;
    xlabel = "o",
    ylabel = "w",
    title  = "Proportion Feasible Runs",
    colorbar_title = "Proportion",
    aspect_ratio = :equal
)

Z_3= Matrix(eig_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z_3;
    xlabel = "o",
    ylabel = "w",
    title  = "Median Eigenvalue (Feasible Runs)",
    colorbar_title = "Median Eigenvalue",
    aspect_ratio = :equal
)

##sanity checks
# should be ≤ 1
maximum(summary_ow.prop_feasible)

# cells with zero feasible runs
filter(:n_feasible => ==(0), summary_ow)

# compare raw vs filtered CV
describe(all_runs_unforced_df.cv_total)
describe(all_runs_unforced_df.cv_total[all_runs_unforced_df.feasible])
##so, here i now have a matrix of vectors for each combination of o and w that gives coexistence 
##another alternative rather than the % of times that get lowest cv, is in this setup, since the parameter sets for each comb of o-w are different, can calculate lets say min or mean cv of total harvest, and then plot that on heatmap

# ## =================================================================================
# ## Forced robustness grid: reuse same feasible nuisances per (o,w)
# ## (H is treated as a focal parameter here)
# ## =================================================================================

####With H in model 

H_grid = range(bounds[:H]...; length=2)
# # Allocate a cell array that stores all runs for each (o,w,H) - for forced model
 runs_stab_forced = [NamedTuple[] for _ in 1:No, _ in 1:Nw, _ in 1:length(H_grid)]
 
 @info "Running robustness grid for FORCED model..."
 @threads for idx in 1:(No * Nw * length(H_grid))
     tmp = idx - 1
     io = fld(tmp, (Nw * length(H_grid))) + 1
     rem1 = tmp % (Nw * length(H_grid))
     iw = fld(rem1, length(H_grid)) + 1
     iH = (rem1 % length(H_grid)) + 1

     o = o_grid[io]
     w = w_grid[iw]
     H = H_grid[iH]

     P0 = 0.25
     local_runs = NamedTuple[]
     #local_pool = feasible_pool_grid[io, iw]  # same feasible nuisances as unforced

     for s in param_library
         try
            ν = s.params
             p = make_par_forced(o, w, H, ν)
             out_1 = equilibrium_forced_2(p, P0; 
                                          t_warmup = 300.0, 
                                          t_eval   = 350.0,
                                          ngrid    = 300,
                                          reltol   = 1e-6,
                                          abstol   = 1e-6)

             rec = (; o, w, H, s...,
                     C1_min = out_1.min[3],
                     C2_min = out_1.min[4],
                     R1_min = out_1.min[1],
                     R2_min = out_1.min[2],
                     cv_total = out_1.cv_total,
                     cv_R1    = out_1.cv_R1,
                     cv_R2    = out_1.cv_R2,
                     cv_C1    = out_1.cv_C1,
                     cv_C2    = out_1.cv_C2,
                     cv_G     = out_1.cv_G)

             push!(local_runs, rec)
         catch err
             @warn "Forced fail at (o=$o, w=$w, H=$H): $err"
         end
     end

     runs_stab_forced[io, iw, iH] = local_runs
 end
 @info "Done FORCED robustness grid."

 cell = runs_stab_forced[6,6,1]
 cell[2]
 cell[3]


 ##flatten results matrix to a DataFrame
all_runs_forced = reduce(vcat, vec(runs_stab_forced))
all_runs_forced = [(; run_id = i, r...) for (i, r) in enumerate(all_runs_forced)]
all_runs_forced_df = DataFrame(all_runs_forced)
names(all_runs_forced_df)

df_flat = hcat(
    select(all_runs_forced_df, Not(:params)),
    DataFrame(all_runs_forced_df.params);
    makeunique = true
)
rename!(df_flat, Dict(:H => :H_forcing, :H_1 => :H_nuisance))

names(df_flat)
write_parquet("runs_stab_forced.parquet", df_flat)

all_runs_forced_df
df_flat

filter(n -> occursin("H", String(n)), names(df_flat))
# If you see H and H_1 (or similar), check whether they match:
all(df_flat.H .== df_flat.H_1)
describe(df_flat.H); describe(df_flat.H_1)

# And confirm your H grid is actually represented:
combine(groupby(df_flat, :H), nrow)
combine(groupby(df_flat, :H_1), nrow)


##okay sick, so now, want to calculate proportion of feasible, and then median CV
tol = 1e-3

df_flat.feasible = (
    (df_flat.R1_min .> tol) .&
    (df_flat.R2_min .> tol) .&
    (df_flat.C1_min .> tol) .&
    (df_flat.C2_min .> tol)
)
summary_ow = combine(
    groupby(df_flat, [:o, :w, :H_forcing]),
    :feasible => mean => :prop_feasible,
    :feasible => sum  => :n_feasible,
    :cv_total => length => :n_total,
    [:cv_total, :feasible] =>
        ((cv, f) -> median(cv[f])) => :median_cv_feasible,
)


##transpose back into a grid for Plotting
summary_ow_H0 = filter(:H_forcing => ==(0), summary_ow)
summary = sort(summary_ow_H0, [:w, :o])   # or [:o,:w] depending on your convention

prop_grid = unstack(summary, :w, :o, :prop_feasible)
median_cv_grid = unstack(summary, :w, :o, :median_cv_feasible)

spread_ow = combine(groupby(df_flat, [:o, :w, :H_forcing]),
     [:cv_total, :feasible] => ((cv,f)-> any(f) ? (quantile(cv[f],0.9)-quantile(cv[f],0.1)) : missing) => :q90_q10_cv_feasible
)
spread_H0 = filter(:H_forcing => ==(0), spread_ow)
spread_grid = unstack(spread_H0, :w, :o, :q90_q10_cv_feasible)

# y-axis (rows)
w_vals = median_cv_grid.w

# x-axis (column names, skipping :w)
o_syms = names(median_cv_grid)[2:end]
o_vals = parse.(Float64, string.(o_syms))

# z matrix
Z = Matrix(median_cv_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z;
    xlabel = "o",
    ylabel = "w",
    title  = "Median CV (feasible runs)",
    colorbar_title = "median CV",
    aspect_ratio = :equal
)


# z matrix
Z_2= Matrix(prop_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z_2;
    xlabel = "o",
    ylabel = "w",
    title  = "Proportion Feasible Runs",
    colorbar_title = "Proportion",
    aspect_ratio = :equal
)

# z matrix
Z_3 = Matrix(spread_grid[:, 2:end])

heatmap(
    o_vals,
    w_vals,
    Z_3;
    xlabel = "o",
    ylabel = "w",
    title  = "90 quantile",
    colorbar_title = "90 quantile",
    aspect_ratio = :equal
)


##So when i run the forced model this way, the heat maps look random - why is this? 



###testing to see if maybe issue is how long evaluation window is 
# -----------------------
# User settings
# -----------------------
rng = MersenneTwister(1)

H_forcing = 0.0         # pick one forced level to test
P0 = 0.25
u0 = [1.5, 1.5, 1.0, 1.0, 0.25]

t_warmup = 300.0
t_eval_short = 350.0
t_eval_long  = 2000.0

ngrid_short = 300
ngrid_long  = 2000

n_cells = 5              # test 3–5 cells
n_nuis  = 5              # test 5 nuisance sets

# -----------------------
# Choose cells + nuisance sets
# -----------------------
cells = [(rand(rng, eachindex(o_grid)), rand(rng, eachindex(w_grid))) for _ in 1:n_cells]
nuis_idx = rand(rng, 1:length(param_library), n_nuis)

# -----------------------
# Run comparison
# -----------------------
rows = NamedTuple[]

for (io, iw) in cells
    o = o_grid[io]
    w = w_grid[iw]

    for sid in nuis_idx
        s = param_library[sid]
        ν = s.params

        # Build full forced parameter struct
        p = make_par_forced(o, w, H_forcing, ν)

        # ---- short window ----
        cvS = try
            outS = fr_cv_forced(p; u0=u0, t_warmup=t_warmup, t_eval=t_eval_short, ngrid=ngrid_short)
            outS.cv_total
        catch err
            missing
        end

        # ---- long window ----
        cvL = try
            outL = fr_cv_forced(p; u0=u0, t_warmup=t_warmup, t_eval=t_eval_long, ngrid=ngrid_long)
            outL.cv_total
        catch err
            missing
        end

        push!(rows, (; io, iw, o, w, H_forcing,
                      sid,
                      cv_short=cvS,
                      cv_long=cvL,
                      cv_ratio = (ismissing(cvS) || ismissing(cvL) || cvL == 0) ? missing : cvS / cvL,
                      cv_diff  = (ismissing(cvS) || ismissing(cvL)) ? missing : cvS - cvL))
    end
end

diag_df = DataFrame(rows)

# Inspect
diag_df




























 ##for each o-w combo, select vector that has lowest cv
runs_stab_forced_2 = dropdims(runs_stab_forced; dims = 3)
min_cv = map(v -> minimum(x -> x.cv_total, v), runs_stab_forced_2)
mean_cv = map(v -> mean(x -> x.cv_total, v), runs_stab_forced_2)
#median_cv = map(v -> me(x -> x.cv_total, v), runs_stab_forced_2)


 heatmap(o_grid, w_grid, mean_cv;
         xlabel = "o",
         ylabel = "w",
         title = "Mean CV of total harvest for random feasible nuisance parameter draws",
         colorbar_title = "cv total harvest")

 heatmap(o_grid, w_grid, min_cv;
         xlabel = "o",
         ylabel = "w",
         title = "Minimum CV of total harvest for random feasible nuisance parameter draws",
         colorbar_title = "cv total harvest")


# ## =================================================================================
# ## Filtering for stable coexistence (unforced) and downstream analysis
# ## =================================================================================

# ##filter out only stable coexistence (i.e., where there is some persistence of all species at equilibrium)
# eps = 1e-3
# No, Nw = length(o_grid), length(w_grid)
# predicate = r -> (r.R1_min >eps && r.R2_min > eps && r.C1_min > eps && r.C2_min > eps)
# final_runs_unforced = [ filter(predicate, runs_stab_unforced[i,j])
#                        for i in 1:No, j in 1:Nw]
# #print(final_runs)


# ##See if I can extract cv/eigenvalue from all possible combinations for each set of nuisance parameters
# N_nuisance_final = length(final_runs_unforced[1,1])
# N_nuisance = length(runs_stab_unforced[1,1])
# same_nuisance_results = [ [runs_stab_unforced[i,j][n] for i in 1:No,
#                                          j in 1:Nw]
#                           for n in 1:N_nuisance] ##come back -- do we want this over all nuisance draws or just the final ones?

# ##Seeing if I can plot cv for a single nuisance Draw
# n = 3 ##this is the nuisance draws

# ##extract the data for that draw 
# res_n = same_nuisance_results[n]

# ##turn it back into a matrix over o×w
# cv_total_mat = [res_n[i,j].cv_total for i in 1:No, j in 1:Nw]

# ##now Plot with o on x and w on y, not the transpose of the matrix
# heatmap(o_grid, w_grid, cv_total_mat';
#         xlabel = "o",
#         ylabel = "w",
#         title = "CV of total harvest for nuisance draw $n at H = 0.0",
#         colorbar_title = "cv total harvest")

# ##seeing if i can calculate % of runs where certain combinations of o and w create the lowest cv of total harvest
# No, Nw = size(runs_stab_unforced)
# N_nuisance_final
# final_runs_unforced[1,1][1]
# runs_stab_unforced[1,1][1]

# counts = zeros(Float64, No, Nw)  # use Float64 to split ties fairly if you want

# ##to figure out -- are you doing this for all runs or just the ones w/ values <0.0001? need to double check and think through that filtering again 
# for n in 1:N_nuisance
#     best_cv = Inf
#     winners = Tuple{Int,Int}[]
    
#     # scan all (o,w) for this nuisance draw n
#     for i in 1:No, j in 1:Nw
#         recs = runs_stab_unforced[i,j]

#         # guard in case some cells have fewer draws
#         if n > length(recs)
#             continue
#         end

#         cv_test = recs[n].cv_total

#         if isfinite(cv_test)
#             if cv_test < best_cv - Base.eps(Float64)            # strictly better
#                 best_cv = cv_test
#                 empty!(winners)
#                 push!(winners, (i,j))
#             elseif abs(cv_test - best_cv) ≤ Base.eps(Float64)   # tie
#                 push!(winners, (i,j))
#             end
#         end
#     end

#     # give 1 "vote" split among tied winners
#     credit = 1.0 / max(1, length(winners))
#     for (i,j) in winners
#         counts[i,j] += credit
#     end
# end

# # NOTE: you will want to define n_ids appropriately if you keep this;
# # for now, we use N_nuisance
# pct = 100 .* counts ./ max(1, N_nuisance)

# using Plots
# heatmap(o_grid, w_grid, pct';  # transpose so o→x, w→y
#     xlabel="o", ylabel="w", colorbar_title="Winner %",
#     title="Lowest CV_total frequency at H=0")



# ## The rest of your CV grids, DataFrame aggregation,
# ## boxplots, and global sensitivity analysis code
# ## can follow below essentially unchanged.

# ##Global sensitivity analysis 
# # Parameter order for GS (example: only focal ones)
# gs_syms = (:o, :w, :H)
# lb = [bounds[s][1] for s in gs_syms]
# ub = [bounds[s][2] for s in gs_syms]

# # Fix nuisance at random draws inside g, or better: average over a small inner sample for each x
# function g_focal(x)
#     o, w, H = x
#     P0 = 0.25

#     # inner averaging over nuisance for smoother response
#     rng  = MersenneTwister(2025)           # fixed per eval for reproducibility
#     reps = 8
#     νs   = sample_nuisance(reps; rng)

#     vals = Float64[]
#     for i in 1:reps
#         p   = make_par(o, w, H, νs[i])
#         out = equilibrium_unforced(p, P0)  # your NamedTuple
#         v   = out.cv[3]                    # pick SCALAR: CV of C1
#         if isfinite(v)
#             push!(vals, v)
#         end
#     end
#     return isempty(vals) ? NaN : median(vals)  # or mean(vals)
# end

# # Sobol total indices
# N = 50
# sob = GlobalSensitivity.sobol_sensitivity(g_focal, lb, ub, N; second_order=false)
# # sob.ST gives total-effect indices for (o,w,H)
