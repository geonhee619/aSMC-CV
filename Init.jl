using DataFrames, CSV, JLD
using Dates, LaTeXStrings
using LinearAlgebra, NamedArrays
using Random, Distributions, StatsBase, StatsFuns
using Plots, StatsPlots, ProgressMeter
using AdvancedHMC, LogDensityProblems
using PSIS, Roots

vecvec2mat(_v) = reduce(hcat, _v)' |> Matrix
mat2vecvec(_m) = [_m[i,:] for i in 1:size(_m,1)]
symmetric(_m) = Symmetric((_m + _m') ./ 2)
_round(_x::Float64)::Float64 = round(_x; digits=3)
_sum(_M; dims::Int64) = dropdims(sum(_M; dims=dims); dims=dims)
_mean(_M; dims::Int64) = dropdims(mean(_M; dims=dims); dims=dims)

default(size=(500,200))

(^)(f::Function, i::Int) = i==1 ? f : x->(f^(i-1))(f(x))

;