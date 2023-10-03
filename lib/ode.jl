function densitydependent(du, u, p, t)
    S, r, δ, h, β, M, P, C, ρ, ν = p

    s₀ = u[1:S]
    i₀ = u[(S+1):end]
    H₀ = s₀ .+ i₀

    mutualists = .!iszero.(M)
    predators = P .>= 0.0
    preys = transpose(predators)

    # Terms for each interaction type
    predation = (P .* predators) * (H₀ ./ (h .+ H₀))
    prey = (P .* preys) * H₀
    mutualism = (M .* mutualists) * (H₀ ./ (h .+ H₀))
    competition = C * H₀

    # ODE model
    du[1:S] .= s₀ .* (r - δ .* H₀ + competition + mutualism + predation) + (s₀ ./ (h .+ H₀)) .* prey - s₀ .* (β * i₀) + ρ .* i₀
    du[(S+1):end] .= i₀ .* (r - δ .* H₀ + competition + mutualism + predation) + (i₀ ./ (h .+ H₀)) .* prey + s₀ .* (β * i₀) - (ρ + ν) .* i₀
end
