function onesim(S, intprop)
    
    A = zeros(Float64, (4, S, S))
    interactions!(A, Co, K, σₓ, σᵢ, intprop)
    
    β = transmissionmatrix(A, δ)
    
    M = A[2,:,:]
    C = A[3,:,:]
    P = A[4,:,:]

    u₀ = 10.0rand(2S)
    p = (S, r, δ, h, β, M, P, C, ρ, ν)

    prob = ODEProblem(densitydependent, u₀, (0.0, timesteps), p)
    sol = solve(prob, callback=ContinuousCallback(onepopislow, extinguishpop!))

    SI = sol[end]
    Sₜ = SI[1:S]
    Iₜ = SI[(S+1):end]
    Hₜ = Sₜ .+ Iₜ
    
    remain = findall((Sₜ .> 10eps()).*(Iₜ .> 10eps()))
    Sₜ = Sₜ[remain]
    Iₜ = Iₜ[remain]
    Hₜ = Sₜ .+ Iₜ
    
    return (Sₜ, Iₜ, Hₜ)
end