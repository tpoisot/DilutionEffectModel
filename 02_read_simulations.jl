using DataFrames
import CSV
using DelimitedFiles
using Statistics

output_path = joinpath(dirname(Base.active_project()), "outputs")
result_path = joinpath(dirname(Base.active_project()), "results")
if ~ispath(result_path)
    mkpath(result_path)
end

outputs = readdir(joinpath(output_path), join=true)
filter!(!contains("gitkeep"), outputs)
all_simulations = vcat(readdir.(outputs; join=true)...)
results = [DataFrame() for thr in 1:Threads.nthreads()]

for simulation in all_simulations
    simid, repid = splitpath(simulation)[(end-1):end]
    timeseries = readdlm(joinpath(simulation, "timeseries.dat"))
    competition = readdlm(joinpath(simulation, "competition.dat"))
    infection = readdlm(joinpath(simulation, "infection.dat"))
    mutualism = readdlm(joinpath(simulation, "mutualism.dat"))
    predation = readdlm(joinpath(simulation, "predation.dat"))

    endpoint = timeseries[2:end,end]
    
    sp = round(Int, length(endpoint)/2)

    Sₜ = endpoint[1:sp]
    Iₜ = endpoint[(sp+1):end]
    Hₜ = Sₜ + Iₜ

    keep = findall(Hₜ .> 0.1)

    Sₜ = Sₜ[keep]
    Iₜ = Iₜ[keep]
    Hₜ = Hₜ[keep]
    pₜ = Hₜ./sum(Hₜ)

    prevalence = Iₜ./Hₜ
    degree = vec(sum(infection[keep,keep].>0.0, dims=1))

    push!(results[Threads.threadid()], (
        parameters = simid,
        replicate = repid,
        richness = length(Hₜ),
        susceptibles = sum(Sₜ),
        infectious = sum(Iₜ),
        population = sum(Hₜ),
        correlation = sum(Iₜ) == 0 ? 0.0 : cor(degree, prevalence),
        prevalence = sum(Iₜ)/sum(Hₜ),
        diversity = sum(-pₜ.*log.(pₜ))
    ))

end

CSV.write(joinpath(result_path, "outputs.csv"), vcat(results...))