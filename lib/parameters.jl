S = 25 # Number of species
timesteps = 500 # Timesteps for the simulation
Co = 0.5 # Connectance of the matrix
K = 100.0 # Carrying capacity
σₓ = 0.2 # Dispersal for biotic interactions
σᵢ = 0.001 # Dispersal for disease transmission
δ = 0.02 # Self-regulation parameter for infection
r = ones(S) # Growth rate
ν = 0.05ones(S) # Virulence
ρ = 0.05ones(S) # Recovery
h = 30.0 # Half-saturation constant
