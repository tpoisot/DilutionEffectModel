
function shannon(x)
    p = x ./ sum(x)
    return - sum(p .* log.(p))
end

range(x...) = x ./ sum(x)