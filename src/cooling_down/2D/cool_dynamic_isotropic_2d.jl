# # Cooling process of 2D plate 
# In this example we simulate the cooling process of a 2D plate with dynamic (temperature-dependent) and isotropic material properties. 

using Hestia 

# Thermal conductivity: λ(Θ)=10 + 0.1Θ
λ = [10.0, 0.1]  # Thermal conductivity: temperature-dependend

# Mass density
ρ = [7800.0]     # Mass density: constant

# Specific heat capacity: c(Θ)=330+0.4Θ
c = [330.0, 0.4] 

# Length and width of plate
L, W = 0.2, 0.1

# Discretized points
Nx, Ny = 40, 20    
Ntotal = Nx*Ny

# Material property with isotropic and temperature-dependend (dynamic) behavior
property = DynamicIsotropic(λ, ρ, c)
plate    = HeatPlate(L, W, Nx, Ny)

# Emission / Stefan-Boltzmann Boundary Condition
# - heat transfer coefficient: k=10
# - heat radiation with emssivity ϵ=0.6
# - ambient temperature  Θamb=300 Kelvin
emission  = Emission(10.0, 0.6, 300.0)

# The emission is assumed on boundary sides west, east and north. Boundary side south is automatically initialized with k=0, ϵ=0.
boundary = Boundary(plate)
setEmission!(boundary, emission, :west )
setEmission!(boundary, emission, :east)
setEmission!(boundary, emission, :north  )

# The discrete heat conduction model is defined as right-hand side function to be solved with a numerical integrator.
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, plate, property, boundary)
end

# Set initial temperature and time span for integration
θinit = 600.0*ones(Ntotal)
tspan = (0.0, 2000.0)

# Save temperates at time each time step
ts    = 0.2       

# The differential equation is solved using the `OrdinaryDiffEq` library.
# Runge-Kutta type integration method `KenCarp5` to solve stiff problems.
# See: https://en.wikipedia.org/wiki/Stiff_equation
import OrdinaryDiffEq
alg = OrdinaryDiffEq.KenCarp5()
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,alg, saveat=ts)

# Plotting the resulting temperature evolution as heatmap
using Plots
xgrid = L/(2Nx) : L/Nx : L;
ygrid = W/(2Ny) : W/Ny : W;
heatmap(xgrid, ygrid, reshape(sol[end], Nx, Ny)', ylabel="Width in [m]", xlabel="Length in [m]", title="Final temperature distribution", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)
