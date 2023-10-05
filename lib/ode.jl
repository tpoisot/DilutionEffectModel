function RwMASIR(du, u, p, t)
    
    S = length(p.r)

    # Populations at this time step
    s₀ = u[1:S]
    i₀ = u[(S+1):end]
    # We get the total population here too
    H₀ = s₀ .+ i₀

    # This is used to speed-up the lookup of interactions
    mutualists = .!iszero.(p.M)
    predators = p.P .>= 0.0
    preys = transpose(predators)

    # Terms for each interaction type
    predation = (p.P .* predators) * (H₀ ./ (p.h .+ H₀))
    prey = (p.P .* preys) * H₀
    mutualism = (p.M .* mutualists) * (H₀ ./ (p.h .+ H₀))
    competition = p.C * H₀

    # Transmission term
    ℬ = p.density ? (p.β * i₀) : (p.β * (i₀./sum(H₀)))

    # ODE model
    du[1:S] .= s₀ .* (p.r - p.δ .* H₀ + competition + mutualism + predation) + (s₀ ./ (p.h .+ H₀)) .* prey - s₀ .* ℬ + p.ρ .* i₀
    du[(S+1):end] .= i₀ .* (p.r - p.δ .* H₀ + competition + mutualism + predation) + (i₀ ./ (p.h .+ H₀)) .* prey + s₀ .* ℬ - (p.ρ + p.ν) .* i₀
end