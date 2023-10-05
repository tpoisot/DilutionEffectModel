function onepopislow(u, t, integrator)
    return iszero(count(u .< 100eps()))
end

function extinguishpop!(integrator)
    for i in eachindex(integrator.u)
        if integrator.u[i] < 100eps()
            integrator.u[i] = 0.0
        end
    end
end

function onesim(S, ip, mp)

    A = zeros(Float64, (4, S, S))
    interactions!(A, mp, ip)

    β = transmissionmatrix(A, mp)

    M = A[2, :, :]
    C = A[3, :, :]
    P = A[4, :, :]

    u₀ = 10.0rand(2S)

    # This is the full set of parameters we want
    rp = (mp..., M=M, C=C, P=P, β=β)

    prob = ODEProblem(RwMASIR, u₀, (0.0, mp.timesteps), rp)
    sol = solve(prob, callback=ContinuousCallback(onepopislow, extinguishpop!))

    # Get the solution as an array
    t, uₜ = sol.t, Array(sol)
    
    # Remove every population size under the threshold of 0.1
    uₜ[findall(uₜ.<0.1)] .= 0.0

    # Return
    return (vcat(t', uₜ), M, C, P, β)
end