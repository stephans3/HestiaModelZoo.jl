# # Cooling process of 1D rod 
# In this example we simulate the cooling process of a 1D rod with dynamic (temperature-dependent) and isotropic material properties. 
using Hestia 


# Thermal conductivity: λ(Θ)=8 + 0.1Θ
λ = [8.0, 0.1]

# Mass density
ρ = [7800.0]

# Specific heat capacity: c(Θ)=330+0.5Θ
c = [330, 0.5]  

# Length and discretized points
L = 0.2     # Length
Nx = 40    # Number of elements: x direction

# Material property with isotropic and temperature-dependend (dynamic) behavior
property = DynamicIsotropic(λ, ρ, c)
heatrod  = HeatRod(L, Nx)


# Emission / Stefan-Boltzmann Boundary Condition
# - heat transfer coefficient: k=5
# - heat radiation with emssivity ϵ=0.5
# - ambient temperature  Θamb=300
emission = Emission(5.0, 0.5, 300.0)

# The emission is assumed only the right side of the rod ":east". The left side is automatically initialized with k=0, ϵ=0.
boundary = Boundary(heatrod)
setEmission!(boundary, emission, :east)

# The discrete heat conduction model is defined as right-hand side function to be solved with a numerical integrator.
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, heatrod, property, boundary)
end

# Set initial temperature and time span for integration
θinit = 600*ones(Nx)
tspan = (0.0, 2000.0)

# Save solution at
ts = 1.0

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
tgrid = 0 : ts : tspan[2];
heatmap(tgrid, xgrid, sol[:,:], xlabel="Time in [s]", ylabel="Length in [m]", title="Cooling Process of 1D Rod", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)