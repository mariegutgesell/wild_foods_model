##Code to look at try out a type of sensitivty analysis of harvest rate CV across parameter ranges 
##Date Initiated: October 15, 2025
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
#include("wf_model_eqs_subsidy_Pfixed.jl") ##model equations with P held constant, unique parameters per trophic level, active and passive omnivory parameter structures
include("wf_model_harvest_rate_analysis.jl") ##function to calculate cv of total harvest
##not sure if need to redefine model_par, i think may be okay to call from my wf_model code 

# --- Bounds: focal vs nuisance ---
focal_syms    = (:o, :w, :H)   # change to your three
nuisance_syms = (:r, :K, :aR_C, :aR_P, :aC_P, :aG_P, :hR_P, :hR_C, :hC_P, :hG_P, :e, :mC, :G)

# Ranges (examples—replace with yours)
bounds = Dict(
    :o => (0.0, 1.0),
    :w => (0.0, 1.0),
    :H => (0.0, 1.0),
    :r => (0.2, 2.0),
    :K => (1.0, 5.0),
    :aR_P => (0.1, 3.0),
    :aR_C => (0.1, 3.0),
    :aC_P => (0.1, 3.0),
    :aG_P => (0.1, 3.0),
    :hR_C => (0.1, 3.0),
    :hR_P => (0.1, 3.0),
    :hC_P => (0.1, 3.0),
    :hG_P => (0.1, 3.0),
    :e   => (0.1, 0.9),
    :mC  => (0.01, 1.5),
   # :mP  => (0.01, 1.5),
    :G => (0.0, 10.0),
  #  :G_base => (0.0, 10.0),
  #  :l  => (0.1, 5.0),
  #  :pf => (0.0, 10.0),
  #  :D => (0.0, 1.0)
)

# Optional: which are log-scaled? - good for ones that span orders of magnitude
logscale = Set([:aR_P, :aR_C, :aG_P, :aC_P, :K])  # e.g., Set([:aR_P, :aC_P, :aG_P])


##So, want to randomly select from these ranges when running my model
##Then, evaluate stability - and somehow filter out ones that are biologically unrealistic/give wild eigenvalue or whatever
##then think, how do i focus on my own parameters im interested in? and outcome? 


##Helper functions - 
# Map unit cube sample x∈[0,1] to parameter in [lo,hi] (linear or log)
map_to_range(x, lo, hi; logscaled=false) =
    logscaled ? exp(log(lo) + x*(log(hi) - log(lo))) : (lo + x*(hi - lo)) ##if logscale is true, use this line, otherwise use top line
##this is a helper function to transform all of the random values from sampling scheme (e.g., latin hypercube) which are always between 0 and 1 to meaningful model parameters
##x is the sampled value (between 0 and 1) and lo is the low end of range, and hi is the high end of range 
##log scaled is to give a log-uniform distribution, good for parameters that span orders of magnitude 


# Draw N samples of nuisance params with LHS
function sample_nuisance(N::Int; rng=Random.default_rng())
    d = length(nuisance_syms)
    X = QuasiMonteCarlo.sample(d, N, LatinHypercubeSample())   # N × d in [0,1]
    samples = Vector{NamedTuple}(undef, N)
    for i in 1:N
        pairs = ntuple(j -> begin
            s = nuisance_syms[j]
            lo, hi = bounds[s]
            val = map_to_range(X[i,j], lo, hi; logscaled = (s in logscale))
            (s => val)
        end, d)
        samples[i] = NamedTuple(pairs)
    end
    samples
end



##Calculate robustness surface for 3 focal parameters

# Grids for the 3 focal parameters
o_grid = range(bounds[:o]...; length=11)
w_grid = range(bounds[:w]...; length=11)
H_grid = range(bounds[:H]...; length=11)

Nrep = 2  # random nuisance samples per grid point (tune)

# Pre-sample nuisance once to reuse (or sample per cell if you prefer)
rng = MersenneTwister(42)
nuisance_pool = sample_nuisance(Nrep; rng)

# Allocate result arrays: (|o|, |w|, |H|)
#CV_med = Array{Float64}(undef, length(o_grid), length(w_grid), length(H_grid))
#CV_iqr = similar(CV_med)
#stab_med = similar(CV_med)
#stab_iqr = similar(CV_med)

# Build parameter from focal + nuisance sample
const p0 = ModelPar_active()

##this approach copies the parameters from pO and only overrides the ones indicated after the ; (so keeps the function parameters)
function make_par(o, w, H, ν::NamedTuple)
    return ModelPar_active(p0; o=o, w=w, H=H,
        r=ν.r, K=ν.K,
        aR_P=ν.aR_P, aC_P=ν.aC_P, aG_P=ν.aG_P, aR_C=ν.aR_C,  # <- check names
        hR_P=ν.hR_P, hC_P=ν.hC_P, hG_P=ν.hG_P, hR_C=ν.hR_C,  # <- check names
        e=ν.e, mC=ν.mC, 
        G=ν.G)
end

# Allocate a cell array that stores all runs for each (o,w,H)
runs_stab = [NamedTuple[] for _ in 1:length(o_grid), _ in 1:length(w_grid), _ in 1:length(H_grid)]

@info "Running robustness grid..."
for (io, o) in enumerate(o_grid), (iw, w) in enumerate(w_grid), (iH, H) in enumerate(H_grid)
    #P0 = 0.25
    u0 = [1.5, 1.5, 1.0, 1.0, 0.25]
    stability_df = NamedTuple[]
    for s in nuisance_pool
        p = make_par(o, w, H, s)
        try
           
            out_1 = fr_cv_unforced(p; u0=u0, t_warmup=300.0, t_eval=500.0, ngrid=800)  # e.g. returns (cv=..., λ1=..., ...)
            # Store the full parameter set + outputs in one record
            rec = (; o, w, H, s..., cv_total = out.cv_total, cv_R1 = out.cv_R1, cv_R2 = out.cv_R2, cv_C1 = out.cv_C1,cv_C2 = out.cv_C2,cv_G = out.cv_G)
            push!(stability_df, rec)
            push!(runs_stab[io,iw,iH], rec)
        catch err
            @warn "Fail at (o=$o, w=$w, H=$H): $err"
        end
    end

    # summarize
    cv_total = getindex.(stability_df, :cv_total)
   cv_R1 = getindex.(stability_df, :cv_R1)
    cv_R2 = getindex.(stability_df, :cv_R2)
     cv_C1 = getindex.(stability_df, :cv_C1)
      cv_C2 = getindex.(stability_df, :cv_C2)
       cv_G = getindex.(stability_df, :cv_G)

   # good = .!(isnan.(cvC1) .| isinf.(cvC1) .| isnan.(lam1) .| isinf.(lam1))
   # cvC1 = cvC1[good]; lam1 = lam1[good]

   # CV_med[io,iw,iH]  = median(cvC1)
   # CV_iqr[io,iw,iH]  = quantile(cvC1, 0.75) - quantile(cvC1, 0.25)
   # stab_med[io,iw,iH] = median(lam1)
   # stab_iqr[io,iw,iH] = quantile(lam1, 0.75) - quantile(lam1, 0.25)
end
@info "Done."

cell = runs_stab[1,1,1]
cell[2]
print(CV_med)

k = round(Int, 11)
heatmap(o_grid, w_grid, (cvC1[:, :, k])',
        xlabel="o", ylabel="w", colorbar_title="CV C1",
        title="cvC1 eigenvalue at H=$(round(H_grid[k],digits=2))")

##See if I can extract cv/eigenvalue from all possible combinations for each set of nuisance parameters
N_nuisance = length(runs_stab[1,1,1])
same_nuisance_results = [ [runs_stab[i,j,k][n] for i in 1:length(o_grid),
                                         j in 1:length(w_grid),
                                         k in 1:length(H_grid)]
                          for n in 1:N_nuisance ]

##Seeing if I can plot cv for a single nuisance Draw
n = 1 ##this is the nuisance draws
kH = 10 ##this is the slice of K 

##extract the data for that draw 
res_n = same_nuisance_results[n]

##turn it back into a matrix over oxw for that fixed H 
No, Nw, NH = length(o_grid), length(w_grid), length(H_grid)
CV = [res_n[(i-1)*Nw*NH + (j-1)*NH + kH].cv_total for j in 1:Nw, i in 1:No]

##now Plot 
heatmap(o_grid, w_grid, CV; xlabel = "o", ylabel = "w", title = "CV for nuisance draw $n at H=$(H_grid[kH])", colorbar_title = "CV of total harvest")


##eigenvalue 
eig = [res_n[(i-1)*Nw*NH + (j-1)*NH + kH].λ1 for j in 1:Nw, i in 1:No]
##now Plot 
heatmap(o_grid, w_grid, eig; xlabel = "o", ylabel = "w", title = "Max eigenvalue for nuisance draw $n at H=$(H_grid[kH])", colorbar_title = "max eigenvalue")



##Create grid of CV and eigenvalue responses for all combinations 
R = length(runs_stab[1,1,1])  # number of nuisance runs per cell (10)

cv_grid = Array{Float64}(undef, length(o_grid), length(w_grid), length(H_grid), R)
lam_grid = similar(cv_grid)

for io in eachindex(o_grid), iw in eachindex(w_grid), iH in eachindex(H_grid)
    cell = runs_stab[io, iw, iH]            # Vector{NamedTuple}
    @assert length(cell) == R "Unequal runs per cell at ($io,$iw,$iH)"
    for ir in 1:R
        cv_grid[io,iw,iH,ir]  = cell[ir].cvC1
        lam_grid[io,iw,iH,ir] = cell[ir].λ1
    end
end

cell = runs_stab[1,2,1]
@assert cv_grid[1,1,1,1] == cell[1].cvC1
@assert lam_grid[1,1,1,1] == cell[1].λ1

@show cv_grid[1,1,1,2]

rows = NamedTuple[]
for io in eachindex(o_grid), iw in eachindex(w_grid), iH in eachindex(H_grid), ir in 1:R
    push!(rows, (; o=o_grid[io], w=w_grid[iw], H=H_grid[iH],
                  cvC1=cv_grid[io,iw,iH,ir], λ1=lam_grid[io,iw,iH,ir],
                  irun=ir))
end
df_all = DataFrame(rows)


cv_no_na = filter!(r -> isfinite(r.cvC1), df_all)

@df cv_no_na boxplot(string.(:o), :cvC1;
                     xlabel = "o",
                     ylabel = "CV(C1)",
                     title = "CV(C1) across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df cv_no_na boxplot(string.(:w), :cvC1;
                     xlabel = "w",
                     ylabel = "CV(C1)",
                     title = "CV(C1) across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df cv_no_na boxplot(string.(:H), :cvC1;
                     xlabel = "H",
                     ylabel = "CV(C1)",
                     title = "CV(C1) across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)                     

eig_no_na = filter!(r -> isfinite(r.λ1), df_all)
@df df_all scatter(string.(:o), :λ1;
                     xlabel = "o",
                     ylabel = "λ1",
                     title = "λ1 across all nuisance runs",
                     fillalpha = 0.5,
                     legend = false,
                     linewidth = 0.8,
                     whisker_width = 0.6)

@df df_all scatter(string.(:w), :λ1;
                     xlabel = "w",
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























