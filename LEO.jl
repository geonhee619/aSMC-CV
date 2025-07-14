Base.@kwdef mutable struct Params
    β::Matrix{Float64} = zeros(data.T, 3)
    β_0::Vector{Float64} = zeros(3)
    Σ_β::Matrix = I(3)
    Σ_y::Matrix = I(data.τ_axis |> length)
end

function Σy_gibbs(c::Params)::Params
    
    ν_0, S_0 = 2data.K, Matrix(I, data.K,data.K)
    
    E = Vector(undef, data.T)
    t = 1
    for t in eachindex(E)
        ŷ = data.X * c.β[t,:]
        e = data.y[t,:] - ŷ
        E[t] = e * e'
    end
    Σ_y_sample = InverseWishart(ν_0 + length(E), S_0 + sum(E)) |> rand
    
    c.Σ_y = Σ_y_sample
    c
end

function Σβ_gibbs(c::Params)::Params
    
    K = 3
    ν_0, S_0 = 2K, Matrix(I, K, K)
    
    E = Vector(undef, data.T)
    e = c.β[1,:] - c.β_0
    E[1] = e * e'
    for t in 2:data.T
        e = c.β[t,:] - c.β[t-1,:]
        E[t] = e * e'
    end
    Σ_β_sample = InverseWishart(
        ν_0 + length(E),
        S_0 + sum(E)
    ) |> rand

    c.Σ_β = Σ_β_sample
    c
end

function β0_gibbs(c::Params)::Params
    
    m_0, J_0 = zeros(3), I(3)/10
    h_0 = J_0 * m_0
    
    Σ_β_inv = inv(c.Σ_β)
    β_0_sample = MvNormalCanon(
        h_0 + Σ_β_inv * c.β[1,:],
        (J_0 + Σ_β_inv) |> symmetric
    ) |> rand
    
    c.β_0 = β_0_sample
    c
end

function β_ffbs(c::Params)::Params
    
    K, ϵ = 3, 1e-10
    
    """Filtering"""
    # initialize state moments
    θₜˌₜ_mean = zeros(data.T, K) # filtered state variable
    θₜˌₜ_cov  = repeat([zeros(K, K)], data.T) # variance
    θₜˌₜ      = c.β_0
    Σ_θₜˌₜ    = c.Σ_β
    G        = I(K)
    F        = I(K)'
    
    for t in 1:data.T
        # State prior
        θₜˌₜ₋₁ = G * θₜˌₜ # transition
        Σ_θₜˌₜ₋₁ = G * Σ_θₜˌₜ * G' + c.Σ_β
        # 1-step prediction
        F = data.X'
        fₜ = F' * θₜˌₜ₋₁ # measurement
        Qₜ = F' * Σ_θₜˌₜ₋₁ * F + c.Σ_y
        eₜ = data.y[t,:] - fₜ
        # State posterior moments
        θₜˌₜ = θₜˌₜ₋₁ + (Σ_θₜˌₜ₋₁ * F * (Qₜ \ eₜ))
        Σ_θₜˌₜ = Σ_θₜˌₜ₋₁ - (Σ_θₜˌₜ₋₁ * F) * (Σ_θₜˌₜ₋₁ * (F / Qₜ))'
        # Store
        θₜˌₜ_mean[t,:] = θₜˌₜ
        θₜˌₜ_cov[t] = Σ_θₜˌₜ
    end
    
    """Sampling"""
    m = zeros(data.T, K); m[data.T,:] = θₜˌₜ_mean[data.T,:]
    C = repeat([zeros(K, K)], data.T); C[data.T] = θₜˌₜ_cov[data.T]
    θ_sample = zeros(data.T, K)
    θ_sample[data.T,:] = MvNormal(m[data.T,:], (C[data.T] + ϵ*I) |> symmetric) |> rand
    for t in data.T-1:-1:1
        θₜ₊₁ˌₜ = G * θₜˌₜ_mean[t,:]
        Σ_θₜ₊₁ˌₜ = G * θₜˌₜ_cov[t] * G' + c.Σ_β
        m[t,:] = θₜˌₜ_mean[t,:] + (θₜˌₜ_cov[t] * G' * (Σ_θₜ₊₁ˌₜ \ (θ_sample[t+1,:] - θₜ₊₁ˌₜ)))
        C[t] = θₜˌₜ_cov[t] - (θₜˌₜ_cov[t] * G' * (θₜˌₜ_cov[t] * (G' / Σ_θₜ₊₁ˌₜ))')
        θ_sample[t,:] = MvNormal(m[t,:], (C[t] + ϵ*I) |> symmetric) |> rand
    end
    
    c.β = θ_sample
    c
end

function Σy_gibbs_LEO(c::Params; rng::AbstractRNG)::Params
    
    T, ν_0, S_0 = data.T-(1-data.ρ == 0 ? -1 : 0), 2data.K, Matrix(I, data.K,data.K)
    
    E = Vector(undef, T)
    t = 1
    for t in eachindex(E)
        ŷ = data.X * c.β[t,:]
        e = data.y[t,:] - ŷ
        E[t] = (e * e') * (t == data.T ? data.ρ : 1)
    end
    Σ_y_sample = rand(rng,
        InverseWishart(ν_0 + data.T-(1-data.ρ), S_0 + sum(E)))
    
    c.Σ_y = Σ_y_sample
    c
end

function Σβ_gibbs_LEO(c::Params; rng::AbstractRNG)::Params
    
    K = 3
    T, ν_0, S_0 = data.T+(data.ρ == 0 ? -1 : 0), 2K, Matrix(I, K, K)
    
    E = Vector(undef, T)
    e = c.β[1,:] - c.β_0
    E[1] = e * e'
    for t in 2:T
        e = c.β[t,:] - c.β[t-1,:]
        E[t] = e * e'
    end
    Σ_β_sample = rand(rng,
            InverseWishart(
            ν_0 + length(E),
            S_0 + sum(E)
        ))

    c.Σ_β = Σ_β_sample
    c
end

function β0_gibbs_LEO(c::Params; rng::AbstractRNG)::Params
    
    m_0, J_0 = zeros(3), I(3)/10
    h_0 = J_0 * m_0
    
    Σ_β_inv = inv(c.Σ_β)
    β_0_sample = rand(rng, MvNormalCanon(
            h_0 + Σ_β_inv * c.β[1,:],
            (J_0 + Σ_β_inv) |> symmetric
    ))
    
    c.β_0 = β_0_sample
    c
end

function β_ffbs_LEO(c::Params; rng::AbstractRNG)::Params
    
    T, K, ϵ = data.T+(data.ρ == 0 ? -1 : 0), 3, 1e-10
    
    """Filtering"""
    # initialize state moments
    θₜˌₜ_mean = zeros(data.T, K) # filtered state variable
    θₜˌₜ_cov  = repeat([zeros(K, K)], data.T) # variance
    θₜˌₜ      = c.β_0
    Σ_θₜˌₜ    = c.Σ_β
    G        = I(K)
    F        = I(K)'
    
    for t in 1:T
        # State prior
        θₜˌₜ₋₁ = G * θₜˌₜ # transition
        Σ_θₜˌₜ₋₁ = G * Σ_θₜˌₜ * G' + c.Σ_β
        # 1-step prediction
        F = data.X'
        fₜ = F' * θₜˌₜ₋₁ # measurement
        Qₜ = F' * Σ_θₜˌₜ₋₁ * F + c.Σ_y / (t == data.T ? data.ρ : 1)
        eₜ = data.y[t,:] - fₜ
        # State posterior moments
        θₜˌₜ = θₜˌₜ₋₁ + (Σ_θₜˌₜ₋₁ * F * (Qₜ \ eₜ))
        Σ_θₜˌₜ = Σ_θₜˌₜ₋₁ - (Σ_θₜˌₜ₋₁ * F) * (Σ_θₜˌₜ₋₁ * (F / Qₜ))'
        # Store
        θₜˌₜ_mean[t,:] = θₜˌₜ
        θₜˌₜ_cov[t] = Σ_θₜˌₜ
    end
    
    """Sampling"""
    m = zeros(data.T, K); m[T,:] = θₜˌₜ_mean[T,:]
    C = repeat([zeros(K, K)], data.T); C[T] = θₜˌₜ_cov[T]
    θ_sample = zeros(data.T, K)
    θ_sample[T,:] = rand(rng, MvNormal(m[T,:], (C[T] + ϵ*I) |> symmetric))
    for t in T-1:-1:1
        θₜ₊₁ˌₜ = G * θₜˌₜ_mean[t,:]
        Σ_θₜ₊₁ˌₜ = G * θₜˌₜ_cov[t] * G' + c.Σ_β
        m[t,:] = θₜˌₜ_mean[t,:] + (θₜˌₜ_cov[t] * G' * (Σ_θₜ₊₁ˌₜ \ (θ_sample[t+1,:] - θₜ₊₁ˌₜ)))
        C[t] = θₜˌₜ_cov[t] - (θₜˌₜ_cov[t] * G' * (θₜˌₜ_cov[t] * (G' / Σ_θₜ₊₁ˌₜ))')
        θ_sample[t,:] = rand(rng, MvNormal(m[t,:], (C[t] + ϵ*I) |> symmetric))
    end
    
    c.β = θ_sample
    c
end

GibbsScan = β_ffbs ∘ β0_gibbs ∘ Σβ_gibbs ∘ Σy_gibbs
GibbsScan_LEO(c::Params; rng::AbstractRNG)::Params = β_ffbs_LEO(β0_gibbs_LEO(Σβ_gibbs_LEO(Σy_gibbs_LEO(c; rng=rng); rng=rng); rng=rng); rng=rng)

;