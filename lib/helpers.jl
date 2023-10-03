function shannon(x)
    p = x ./ sum(x)
    return -sum(p .* log.(p))
end

function safe_cor(x, y)
    c = cor(x, y)
    return isnan(c) ? 0.0 : c
end

function safe_std(x)
    return length(x) >= 3 ? std(x) : 0.0
end