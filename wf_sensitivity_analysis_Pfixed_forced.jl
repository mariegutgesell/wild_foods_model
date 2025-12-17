##Code to look at try out a type of sensitivty analysis across parameter ranges 
##Date Initiated: October 10, 2025
##Contributor(s): Marie K. Gutgesell

#using Pkg
#Pkg.add("GlobalSensitivity")
#Pkg.add("QuasiMonteCarlo")
#Pkg.add("StatsPlots")

#using DifferentialEquations, ForwardDiff, LinearAlgebra
#using NLsolve
using GlobalSensitivity
using QuasiMonteCarlo   # for LHS/Sobol sampling
#using Statistics
using Random
using StatsPlots

##source model - choose which based on which you want to investigate
include("wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P held constant, unique parameters per trophic level, active and passive omnivory parameter structures

##not sure if need to redefine model_par, i think may be okay to call from my wf_model code 
##1) Sensitivity analysis just for structure -- so keeping H at 0 (constant)
# --- Bounds: focal vs nuisance ---
focal_syms    = (:o, :w, :H)   # for this first one only focusing on o and w 
nuisance_syms = (:r, :K, :aR_C, :aR_P, :aC_P, :aG_P, :hR_P, :hR_C, :hC_P, :hG_P, :e, :mC, :G)

# Ranges (bounds -- based on univariate stability analysis)
#bounds = Dict(
#    :o => (0.0, 1.0),
#    :w => (0.0, 1.0),
#    :H => (0.0, 0.0), ##is this a way to make sure it is always 0 ? 
#    :r => (0.2, 3.0),
#    :K => (1.0, 3.8),
#    :aR_P => (1.0, 2.5),
#    :aR_C => (1.0, 3.0),
#    :aC_P => (0.8, 1.8),
#    :aG_P => (1.0, 3.0),
#    :hR_C => (0.4, 0.8),
#    :hR_P => (0.5, 2.0),
#    :hC_P => (1.0, 2.0),
#    :hG_P => (0.8, 2.0),
#    :e   => (0.7, 1.0),
#    :mC  => (0.6, 1.2),
   # :mP  => (0.01, 1.5),
#    :G => (1.0, 5.0),
  #  :G_base => (0.0, 10.0),
 #   :l  => (0.0, 1.0),
 #   :pf => (0.5, 10.0),
 #   :D => (0.0, 1.0)
#)

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
# Optional: which are log-scaled? - good for ones that span orders of magnitude
#logscale = Set([:aR_P, :aR_C, :aG_P, :aC_P, :mC, :hR_C, :hR_P, :hC_P, :hG_P])  # e.g., Set([:aR_P, :aC_P, :aG_P])

##Helper functions - 
# Map unit cube sample x∈[0,1] to parameter in [lo,hi] (linear or log)
map_to_range(x, lo, hi; logscaled=false) =
    logscaled ? exp(log(lo) + x*(log(hi) - log(lo))) : (lo + x*(hi - lo)) ##if logscale is true, use this line, otherwise use top line
##this is a helper function to transform all of the random values from sampling scheme (e.g., latin hypercube) which are always between 0 and 1 to meaningful model parameters
##x is the sampled value (between 0 and 1) and lo is the low end of range, and hi is the high end of range 
##log scaled is to give a log-uniform distribution, good for parameters that span orders of magnitude 

##create function to test if parameters yield biologically feasible equilibrium
function is_feasible_paramset(p::ModelPar_active; P0=0.25, tol=1e-3)
    try
        out = equilibrium_unforced(p, P0)
        eq = out.eq
        return all(x -> isfinite(x) && x > tol, eq)
    catch
        return false
    end
end



# Draw N samples of nuisance params with LHS
function sample_nuisance(N::Int; rng=Random.default_rng())
    d = length(nuisance_syms)
    X = QuasiMonteCarlo.sample(d, N, LatinHypercubeSample())   # N × d in [0,1]
    samples = Vector{NamedTuple}(undef, N)
    for i in 1:N
        pairs = ntuple(j -> begin
            s = nuisance_syms[j]
            lo, hi = bounds[s]
         #   val = map_to_range(X[i,j], lo, hi; logscaled = (s in logscale)) 
            val = map_to_range(X[i,j], lo, hi) ##trying out without logscaling 
            (s => val)
        end, d)
        samples[i] = NamedTuple(pairs)
    end
    samples
end

##trying new constrained sampling - constraining 
##set seed
Random.seed!(42)
##setting a series of anchors 
anchors = [(0.5,0.5), (0.0,0.0), (0.0,1.0), (1.0,0.0), (1.0,1.0)]

function sample_nuisance_constrained(N::Int; rng=Random.default_rng(), max_tries=10_000)
    d = length(nuisance_syms)
    feasible = NamedTuple[]
    tries = 0

    while length(feasible) < N && tries < max_tries
        batch_size = max(N, 50)
        X = QuasiMonteCarlo.sample(d, batch_size, LatinHypercubeSample())

        for i in 1:batch_size
            tries += 1

            pairs = ntuple(j -> begin
                s = nuisance_syms[j]
                lo, hi = bounds[s]
                val = map_to_range(X[i,j], lo, hi; logscaled = (s in logscale))
                (s => val)
            end, d)
            ν = NamedTuple(pairs)

            ok = false
            for (oa, wa) in anchors
                p = ModelPar_active(p0; o=oa, w=wa, H=0.0,
                    r=ν.r, K=ν.K,
                    aR_P=ν.aR_P, aC_P=ν.aC_P, aG_P=ν.aG_P, aR_C=ν.aR_C,
                    hR_P=ν.hR_P, hC_P=ν.hC_P, hG_P=ν.hG_P, hR_C=ν.hR_C,
                    e=ν.e, mC=ν.mC, G=ν.G)

                if is_feasible_paramset(p)
                    ok = true
                    break
                end
            end

            if ok
                push!(feasible, ν)
                if length(feasible) >= N
                    break
                end
            end

            if tries >= max_tries
                break
            end
        end
    end

    if length(feasible) < N
        @warn "Only found $(length(feasible)) feasible samples after $tries draws"
    end

    return feasible[1:min(end, N)]
end

##constraining: may need to try different constraints - -with and without H - when is 0 and 0.5 
##could randomly choose anchors - then do different draws for different constrained anchor points - could try center and 4 corners 
##could set seed so make sure all draws are the same/saved 
##do we want to to constrain so that all 121 o-w combinations give coexistence - this might be better than anchoring to o=w=0.5
##Calculate robustness surface for 3 focal parameters

# Grids for the 3 focal parameters
o_grid = range(bounds[:o]...; length=11)
w_grid = range(bounds[:w]...; length=11)
H_grid = range(bounds[:H]...; length=3)

Nrep = 3  # random nuisance samples per grid point (tune)

# Pre-sample nuisance once to reuse (or sample per cell if you prefer)
rng = MersenneTwister(42)
p0 = ModelPar_active()
nuisance_pool = sample_nuisance(Nrep; rng)

##need to add conditions so that sampled parameter combinations satisfy inequalities that allow for persistence in an unforced model 


# Allocate result arrays: (|o|, |w|, |H|)
#CV_med = Array{Float64}(undef, length(o_grid), length(w_grid), length(H_grid))
#CV_iqr = similar(CV_med)
#stab_med = similar(CV_med)
#stab_iqr = similar(CV_med)

# Build parameter from focal + nuisance sample
#const p0 = ModelPar_active()

##this approach copies the parameters from pO and only overrides the ones indicated after the ; (so keeps the function parameters)
function make_par(o, w, H, ν::NamedTuple)
    return ModelPar_active(p0; o=o, w=w, H=H,
        r=ν.r, K=ν.K,
        aR_P=ν.aR_P, aC_P=ν.aC_P, aG_P=ν.aG_P, aR_C=ν.aR_C,  # <- check names
        hR_P=ν.hR_P, hC_P=ν.hC_P, hG_P=ν.hG_P, hR_C=ν.hR_C,  # <- check names
        e=ν.e, mC=ν.mC, 
        G=ν.G)
end


##trying version of function tat only calculates cv if feasible and keeps track if not feasible
# Allocate a cell array that stores all runs for each (o,w,H)-- for unforced model 

runs_stab_unforced = [NamedTuple[] for _ in 1:length(o_grid), _ in 1:length(w_grid), _ in 1:length(H_grid)]
@info "Running robustness grid..."
for (io, o) in enumerate(o_grid), (iw, w) in enumerate(w_grid), (iH, H) in enumerate(H_grid)
    P0 = 0.25
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]

    for s in nuisance_pool
        p = make_par(o, w, H, s)

        feasible = false
        out_eq = nothing
        out_cv = nothing

        try
            out_eq = equilibrium_unforced(p, P0)
            feasible = all(x -> isfinite(x) && x > 1e-3, out_eq.eq)

            if feasible
                out_cv = fr_cv_unforced(
                    p; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800
                )
            end
        catch
            feasible = false
        end

        rec = if feasible
            (; o, w, H, s...,
               feasible = true,
               λ1 = out_eq.λ1,
               R1_eq = out_eq.eq[1], R2_eq = out_eq.eq[2],
               C1_eq = out_eq.eq[3], C2_eq = out_eq.eq[4],
               R1_min = out_eq.min[1], R2_min = out_eq.min[2],
               C1_min = out_eq.min[3], C2_min = out_eq.min[4],
               cv_total = out_cv.cv_total,
               cv_R1 = out_cv.cv_R1, cv_R2 = out_cv.cv_R2,
               cv_C1 = out_cv.cv_C1, cv_C2 = out_cv.cv_C2, cv_G = out_cv.cv_G)
        else
            (; o, w, H, s...,
               feasible = false,
               λ1 = missing,
               R1_eq = missing, R2_eq = missing,
               C1_eq = missing, C2_eq = missing,
               R1_min = missing, R2_min = missing,
               C1_min = missing, C2_min = missing,
               cv_total = missing,
               cv_R1 = missing, cv_R2 = missing,
               cv_C1 = missing, cv_C2 = missing, cv_G = missing)
        end

        push!(runs_stab_unforced[io, iw, iH], rec)
    end
end
@info "Done."

cell = runs_stab_unforced[6,6,1]

feasibility_prop = map(runs_stab_unforced) do cell
    mean(getproperty.(cell, :feasible))
end
feasibility_H0 = feasibility_prop[:, :, 1]
feasibility_H05 = feasibility_prop[:, :, 2]

 heatmap(o_grid, w_grid, feasibility_H0';
         xlabel = "o",
         ylabel = "w",
         title = "Feasibility of random nuisance parameter draws and H = 0",
         colorbar_title = "Proportion feasibility")

 heatmap(o_grid, w_grid, feasibility_H05';
         xlabel = "o",
         ylabel = "w",
         title = "Feasibility of random nuisance parameter draws and H = 0.5",
         colorbar_title = "Proportion feasibility")

cv_median = map(runs_stab_unforced) do cell
    cvs = [r.cv_total for r in cell if r.feasible]
    isempty(cvs) ? missing : median(cvs)
end

cv_median_H0 = cv_median[:, :, 1]
cv_median_H05 = cv_median[:, :, 2]

 heatmap(o_grid, w_grid, cv_median_H0';
         xlabel = "o",
         ylabel = "w",
         title = "Median CV of total harvest for random nuisance parameter draws and H = 0",
         colorbar_title = "cv total harvest")

 heatmap(o_grid, w_grid, cv_median_H05';
         xlabel = "o",
         ylabel = "w",
         title = "Median CV of total harvest for random nuisance parameter draws and H = 0.5",
         colorbar_title = "cv total harvest")

# NEXT:
# 1) Increase Nrep to ~50
# 2) Parallelize over nuisance draws
# 3) Add forced-model feasibility + CV
# Allocate a cell array that stores all runs for each (o,w,H) - for forced model
runs_stab_forced = [NamedTuple[] for _ in 1:length(o_grid), _ in 1:length(w_grid), _ in 1:length(H_grid)]

@info "Running robustness grid..."
for (io, o) in enumerate(o_grid), (iw, w) in enumerate(w_grid), (iH, H) in enumerate(H_grid)
    P0 = 0.25
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]

    for s in nuisance_pool
        p = make_par(o, w, H, s)

        feasible = false
        out_cv = nothing

        try
            out_cv = equilibrium_forced_2(p, P0)
            feasible = all(x -> isfinite(x) && x > 1e-3, out_cv.min)
            # feasible = all(x -> isfinite(x)) ##trying this out just to see what happens 
            if feasible
                feasible = true
            end
        catch
            feasible = false
        end

        rec = if feasible
            (; o, w, H, s...,
               feasible = true,
               R1_min = out_cv.min[1], R2_min = out_cv.min[2],
               C1_min = out_cv.min[3], C2_min = out_cv.min[4],
               cv_total = out_cv.cv_total,
               cv_R1 = out_cv.cv_R1, cv_R2 = out_cv.cv_R2,
               cv_C1 = out_cv.cv_C1, cv_C2 = out_cv.cv_C2, cv_G = out_cv.cv_G)
        else
            (; o, w, H, s...,
               feasible = false,
               R1_min = missing, R2_min = missing,
               C1_min = missing, C2_min = missing,
               cv_total = missing,
               cv_R1 = missing, cv_R2 = missing,
               cv_C1 = missing, cv_C2 = missing, cv_G = missing)
        end

        push!(runs_stab_forced[io, iw, iH], rec)
    end
end

@info "Done."

cell = runs_stab_forced[6,6,1]
cell[2]
cell[3]

feasibility_prop = map(runs_stab_forced) do cell
    mean(getproperty.(cell, :feasible))
end
feasibility_H0 = feasibility_prop[:, :, 1]
feasibility_H05 = feasibility_prop[:, :, 2]

 heatmap(o_grid, w_grid, feasibility_H0';
         xlabel = "o",
         ylabel = "w",
         title = "Feasibility of random nuisance parameter draws and H = 0",
         colorbar_title = "Proportion feasibility")

 heatmap(o_grid, w_grid, feasibility_H05';
         xlabel = "o",
         ylabel = "w",
         title = "Feasibility of random nuisance parameter draws and H = 0.5",
         colorbar_title = "Proportion feasibility")

cv_median = map(runs_stab_forced) do cell
    cvs = [r.cv_total for r in cell if r.feasible]
    isempty(cvs) ? missing : median(cvs)
end

cv_median_H0 = cv_median[:, :, 1]
cv_median_H05 = cv_median[:, :, 2]

 heatmap(o_grid, w_grid, cv_median_H0';
         xlabel = "o",
         ylabel = "w",
         title = "Median CV of total harvest for random nuisance parameter draws and H = 0",
         colorbar_title = "cv total harvest")

 heatmap(o_grid, w_grid, cv_median_H05';
         xlabel = "o",
         ylabel = "w",
         title = "Median CV of total harvest for \nrandom nuisance parameter draws and H = 0.5",
         colorbar_title = "cv total harvest")



##filter out only stable coexistence (i.e., where there is some persistence of all species at equilibrium)
eps = 1e-3
No, Nw = length(o_grid), length(w_grid)
predicate = r -> (r.R1_min >eps && r.R2_min > eps && r.C1_min > eps && r.C2_min > eps)
final_runs_unforced = [ filter(predicate, runs_stab_unforced[i,j])
                       for i in 1:No, j in 1:Nw]
#print(final_runs)


##See if I can extract cv/eigenvalue from all possible combinations for each set of nuisance parameters
N_nuisance_final = length(final_runs_unforced[1,1])
N_nuisance = length(runs_stab_unforced[1,1])
same_nuisance_results = [ [runs_stab_unforced[i,j][n] for i in 1:No,
                                         j in 1:Nw]
                          for n in 1:N_nuisance] ##come back -- do we want this over all nuisane draws or just the final ones?

##Seeing if I can plot cv for a single nuisance Draw
n = 3 ##this is the nuisance draws


##extract the data for that draw 
res_n = same_nuisance_results[n]

##turn it back into a matrix over oxw
cv_total_mat = [res_n[i,j].cv_total for i in 1:No, j in 1:Nw]


##now Plot with o on x and w on y, not the transpose of the matrix
heatmap(o_grid, w_grid, cv_total_mat'; xlabel = "o", ylabel = "w", title = "CV of total harvest for nuisance draw $n at H = 0.5", colorbar_title = "cv total harvest")
##this should be at H = 0, come back to this


##seeing if i can calculate % of runs where certain combinations of o and w create the lowest cv of total harvest
No, Nw = size(runs_stab_unforced)
N_nuisance_final
final_runs_unforced[1,1][1]
runs_stab_unforced[1,1][1]


counts = zeros(Float64, No, Nw)  # use Float64 to split ties fairly if you want

##to figure out -- are you doing this for all runs or just the ones w/ values <0.0001? need to double check and think through that filtering again 
for n in 1:N_nuisance
    best_cv = Inf
    winners = Tuple{Int,Int}[]
    
    # scan all (o,w) for this nuisance draw n
    for i in 1:No, j in 1:Nw
        recs = runs_stab_unforced[i,j]

        # guard in case some cells have fewer draws
        if n > length(recs)
            continue
        end

        cv_test = recs[n].cv_total

        if isfinite(cv_test)
            if cv_test < best_cv - Base.eps(Float64)            # strictly better
                best_cv = cv_test
                empty!(winners)
                push!(winners, (i,j))
            elseif abs(cv_test - best_cv) ≤ Base.eps(Float64)   # tie
                push!(winners, (i,j))
            end
        end
    end

    # give 1 "vote" split among tied winners
    credit = 1.0 / max(1, length(winners))
    for (i,j) in winners
        counts[i,j] += credit
    end
end

# Convert to percentages out of total nuisance draws considered at this H
pct = 100 .* counts ./ max(1, length(n_ids))

# Now pct[i,j] is “% of nuisance draws where (o_i, w_j) had the lowest total-harvest CV”
# e.g., heatmap with o on x, w on y:
using Plots
heatmap(o_grid, w_grid, pct';  # transpose so o→x, w→y
    xlabel="o", ylabel="w", colorbar_title="Winner %", title="Lowest CV_total frequency at H=0)")


##okay this is still not working, get lots of 0s etc. need to dig into this more 













##Create grid of CV and eigenvalue responses for all combinations 
R = length(final_runs[1,1,1])  # number of nuisance runs per cell (10)

cv_grid = Array{Float64}(undef, length(o_grid), length(w_grid), length(H_grid), R)


rows = NamedTuple[]
for io in eachindex(o_grid), iw in eachindex(w_grid), iH in eachindex(H_grid), ir in 1:R
    push!(rows, (; o=o_grid[io], w=w_grid[iw], H=H_grid[iH],
                  cv_total=cv_grid[io,iw,iH,ir], λ1=lam_grid[io,iw,iH,ir],
                  irun=ir))
end
df_all = DataFrame(rows)


#cv_no_na = filter!(r -> isfinite(r.cv_total), df_all)

@df df_all boxplot(string.(:o), :cv_total;
                     xlabel = "o",
                     ylabel = "CV total harvest",
                     title = "CV total harvest across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df df_all boxplot(string.(:w), :cv_total;
                     xlabel = "w",
                     ylabel = "CV total harvest",
                     title = "CV total harvest across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df df_all boxplot(string.(:H), :cv_total;
                     xlabel = "H",
                     ylabel = "CV total harvest",
                     title = "CV total harvest across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)                     


@df df_all boxplot(string.(:o), :λ1;
                     xlabel = "o",
                     ylabel = "λ1",
                     title = "λ1 across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df df_all boxplot(string.(:w), :λ1;
                     xlabel = "w",
                     ylabel = "λ1",
                     title = "λ1 across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)
      
@df df_all boxplot(string.(:H), :λ1;
                     xlabel = "H",
                     ylabel = "λ1",
                     title = "λ1 across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)















##Calculate median and IQR for cvC1 and eigenvalue - filter out NAs
df = filter(r -> isfinite(r.cvC1) && isfinite(r.λ1), df_all)
print(df)
# 1) Helpers
IQR(x) = quantile(x, 0.75) - quantile(x, 0.25)

# 2) Aggregate per (o,w,H)
g = groupby(df, [:o])
agg = combine(g,
    :cvC1 => median  => :cv_med,
    :cvC1 => IQR     => :cv_iqr,
    :λ1   => median  => :lam_med,
    :λ1   => IQR     => :lam_iqr,
)


scatter(agg.o, agg.cv_med;
        yerror = agg.cv_iqr ./ 2,   # error bars = half IQR (Q1–Q3 span)
        xlabel = "o",
        ylabel = "Median ± IQR/2 of CV(C1)",
        title = "Median and variability of CV(C1) across o",
        legend = false,
        markersize = 4)

scatter(agg.o, agg.lam_med;
        yerror = agg.cv_iqr ./ 2,   # error bars = half IQR (Q1–Q3 span)
        xlabel = "o",
        ylabel = "Median ± IQR/2 of max eigenvalue",
        title = "Median and variability of max eigenvalue across o",
        legend = false,
        markersize = 4)

##trying grouping by all 3 params of interested
g = groupby(df, [:o, :w, :H])
agg = combine(g,
    :cvC1 => median  => :cv_med,
    :cvC1 => IQR     => :cv_iqr,
    :λ1   => median  => :lam_med,
    :λ1   => IQR     => :lam_iqr,
)

o_levels = sort(unique(df.o))
w_levels = sort(unique(df.w))
H_levels = sort(unique(df.H))

cv_med_grid  = fill(NaN, length(o_levels), length(w_levels), length(H_levels))
cv_iqr_grid  = similar(cv_med_grid)
lam_med_grid = similar(cv_med_grid)
lam_iqr_grid = similar(cv_med_grid)

# Fill the grids
for row in eachrow(agg)
    io = searchsortedfirst(o_levels, row.o)
    iw = searchsortedfirst(w_levels, row.w)
    iH = searchsortedfirst(H_levels, row.H)
    cv_med_grid[io,iw,iH]  = row.cv_med
    cv_iqr_grid[io,iw,iH]  = row.cv_iqr
    lam_med_grid[io,iw,iH] = row.lam_med
    lam_iqr_grid[io,iw,iH] = row.lam_iqr
end

# pick an H slice (e.g., middle)
iH = round(Int, 9)
heatmap(o_levels, w_levels, (cv_med_grid[:,:,iH])',
        xlabel="o", ylabel="w", colorbar_title="CV(C1) median",
        title = "Median CV at H=$(round(H_levels[iH]; digits=3))")


iH = round(Int, 4)
heatmap(o_levels, w_levels, (lam_med_grid[:,:,iH])',
        xlabel="o", ylabel="w", colorbar_title="eigenvalue median",
        title = "Median eigenvalue at H=$(round(H_levels[iH]; digits=3))")


###Looking at 























##Global sensitivity analysis 
# Parameter order for GS (example: only focal ones)
gs_syms = (:o, :w, :H)
lb = [bounds[s][1] for s in gs_syms]
ub = [bounds[s][2] for s in gs_syms]

# Fix nuisance at random draws inside g, or better: average over a small inner sample for each x
function g_focal(x)
    o, w, H = x
    P0 = 0.25

    # inner averaging over nuisance for smoother response
    rng  = MersenneTwister(2025)           # fixed per eval for reproducibility
    reps = 8
    νs   = sample_nuisance(reps; rng)

    vals = Float64[]
    for i in 1:reps
        p   = make_par(o, w, H, νs[i])
        out = equilibrium_unforced(p, P0)  # your NamedTuple
        v   = out.cv[3]                    # pick SCALAR: CV of C1
        if isfinite(v)
            push!(vals, v)
        end
    end
    return isempty(vals) ? NaN : median(vals)  # or mean(vals)
end

# Sobol total indices
N = 50
sob = GlobalSensitivity.sobol_sensitivity(g_focal, lb, ub, N; second_order=false)
# sob.ST gives total-effect indices for (o,w,H)