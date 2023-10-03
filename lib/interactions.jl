inversetransformsampling(𝒟) = quantile(𝒟, rand())

function transmissionmatrix(A, δ)
    β = abs.(dropdims(sum(A[2:end, :, :], dims=1), dims=1))
    β[diagind(β)] .= δ
    β .*= A[1, :, :]
end

function interactions!(A::Array{Float64,3}, Co, K, σₓ, σᵢ, p)

    S = size(A, 2)

    𝒟₁ = Truncated(Normal(0.0, σₓ), 0.0, Inf) # Mutualism and predation
    𝒟₂ = Truncated(Normal(0.0, σₓ / K), 0.0, Inf) # Competition
    𝒟₃ = Truncated(Normal(0.0, σᵢ), 0.0, Inf) # Infection

    interaction_types = [:mutualism, :competition, :predation]

    for i in 1:S
        for j in 1:S
            A[1, i, j] = inversetransformsampling(𝒟₃)
        end
    end

    for i in 1:(S-1)
        for j in (i+1):S
            if rand() <= Co
                interaction = sample(interaction_types, Weights([p...]))
                if interaction == :mutualism
                    A[2, i, j] = inversetransformsampling(𝒟₁)
                    A[2, j, i] = inversetransformsampling(𝒟₁)
                elseif interaction == :competition
                    A[3, i, j] = -inversetransformsampling(𝒟₂)
                    A[3, j, i] = -inversetransformsampling(𝒟₂)
                else
                    A[4, i, j] = rand([-1, 1]) * inversetransformsampling(𝒟₁)
                    A[4, j, i] = -(sign(A[4, i, j])) * inversetransformsampling(𝒟₁)
                end
            end
        end
    end

    return A

end