import Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

# load packages
using XLSX
using StatsBase
using Distributions
using Turing
using Plots
using StatsPlots
using Optim

# read data
root_dir = dirname(@__DIR__)
discount_xl = XLSX.readxlsx(joinpath(root_dir, "data", "discount.xlsx"))
rates = Float64.(discount_xl["Real rates!C7:C231"])
rate_acf = pacf(rates, 1:5)
scatter(rate_acf, linetype=:stem, markersize=4, marker=:circle, markerstrokewidth=3, markeralpha=0.6, color=:red, xlabel="Lag", ylabel="Partial Autocorrelation")

# treat this as an AR(1) model with a drifting level
@model function moving_trend(y)
    ## priors
    # regression for the mean, eta = at+b
    a ~ Normal(0, 0.05)
    b ~ Normal(5, 1)
    # AR(1) residuals
    ρ ~ Uniform(0.75, 1)
    σ ~ truncated(Normal(0, 0.1); lower=0)

    T = length(y)
    y[1] ~ Normal(b, sqrt(σ^2/(1-ρ^2)))
    for t = 2:T
        y[t] ~ Normal(a * (t-1) + b + (ρ * (y[t-1] - a * (t-2) - b)), σ)
    end        
end

model = moving_trend(rates)

chains = sample(model, NUTS(), 10_000)

function predict_values(chains, T, y₀, T₀)
    a = chains[:a].data
    b = chains[:b].data
    ρ = chains[:ρ].data
    σ = chains[:σ].data

    y = zeros(T, length(a))
    T_offset = T₀ - 1798
    for i in eachindex(a)
        y[1, i] = y₀
        for t = 2:T
            y[t, i] = a[i] * (T_offset + t-1) + b[i] + (ρ[i] * (y[t-1, i] - a[i] * (T_offset + t-2) - b[i])) + rand(Normal(0, σ[i]))
        end
    end
    return y
end

y_pred = predict_values(chains, 100, rates[end], 2023)

qs = zeros(3, size(y_pred, 1))

for t = 1:size(y_pred, 1)
    qs[:, t] = quantile(y_pred[t, :], [0.05, 0.5, 0.95])
end

plot(qs[2, :], ribbon=(qs[2, :] .- qs[1, :], qs[3, :] .- qs[2, :]))
plot!(y_pred[:, 1:10], color=:grey, alpha=0.2, label=:false)

discount_factor = exp.(-mapslices(cumsum, y_pred ./ 100; dims=1))
exp_ce_df = mapslices(mean, discount_factor; dims=2)
exp_ce_rate = [100 * (exp_df[t-1] / exp_df[t] - 1) for t=2:length(exp_df)]
plot(2023:2121, exp_ce_rate)