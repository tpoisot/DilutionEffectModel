function onepopislow(u, t, integrator)
    return length(u) - count(u .< 100eps())
end

function extinguishpop!(integrator)
    for i in eachindex(integrator.u)
        if integrator.u[i] < 100eps()
            integrator.u[i] = 0.0
        end
    end
end

function onesim(S, intprop)

    A = zeros(Float64, (4, S, S))
    interactions!(A, Co, K, σₓ, σᵢ, intprop)

    β = transmissionmatrix(A, δ)

    M = A[2, :, :]
    C = A[3, :, :]
    P = A[4, :, :]

    u₀ = 10.0rand(2S)
    p = (S, r, δ, h, β, M, P, C, ρ, ν)

    prob = ODEProblem(densitydependent, u₀, (0.0, timesteps), p)
    sol = solve(prob, callback=ContinuousCallback(onepopislow, extinguishpop!))

    SI = sol[end]
    Sₜ = SI[1:S]
    Iₜ = SI[(S+1):end]
    Hₜ = Sₜ .+ Iₜ

    remain = findall((Sₜ .> 0.1) .* (Iₜ .> 0.1))
    Sₜ = Sₜ[remain]
    Iₜ = Iₜ[remain]
    Hₜ = Sₜ .+ Iₜ

    return (Sₜ, Iₜ, Hₜ)
end