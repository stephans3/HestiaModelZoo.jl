# # Cooling process of 1D rod 
# In this example we simulate the cooling process of a 1D rod static (temperature-independent) and isotropic material properties. 

using Hestia 
# Thermal conductivity
λ = 45.0    

# Mass density
ρ = 7800.0  

# Specific heat capacity
c = 480.0   

# Material property with isotropic and temperature-independend (static) behavior
property = StaticIsotropic(λ, ρ, c)

# Length and discretized points
L = 0.2
Nx = 40

# 1D rod model
heatrod = HeatRod(L, Nx)

# Emission / Stefan-Boltzmann Boundary Condition
# - heat transfer coefficient: k=5
# - heat radiation with emssivity ϵ=0.5
# - ambient temperature  Θamb=300
emission = Emission(5.0, 0.5, 300.0)

# The emission is assumed only the right side of the rod ":east". The left side is automatically initialized with k=0, ϵ=0.
boundary = Boundary(heatrod)
setEmission!(boundary, emission,:east)


# The discrete heat conduction model is defined as right-hand side function to be solved with a numerical integrator.
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, heatrod, property, boundary)
end


# Set initial temperature and time span for integration
θinit = 600*ones(Nx)
tspan = (0.0, 2000.0)

# Save solution at
ts = 1.0

# Numerical integrator: (forward) Euler method
# If you use the (forward) Euler method you need to check the numerical stability.
# See: https://en.wikipedia.org/wiki/Euler_method
# See: https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis

# Sampling time
Δt = 0.2

# Sampling space       
Δx= heatrod.sampling[1]

# Diffusivity
α = λ/(ρ *c)

# The chosen sampling time Δt needs to be smaller than the maximum possible sampling time. 
# If the chosen sampling time equals the maximum possible sampling time than the numerical integration works close to its stability limit.  
max_dt = 0.5*Δx^2/α

if Δt > max_dt
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

# The differential equation is solved using the `OrdinaryDiffEq` library.
import OrdinaryDiffEq

# Numerical integration method: forward Euler method
alg1 = OrdinaryDiffEq.Euler()
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol1 = OrdinaryDiffEq.solve(prob,alg1, dt=Δt, saveat=ts)

# Better use Runge-Kutta solvers which are able to handle stiff equations like KenCarp5.
# See: https://en.wikipedia.org/wiki/Stiff_equation
alg2 = OrdinaryDiffEq.KenCarp5()
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol2 = OrdinaryDiffEq.solve(prob, alg2, saveat=ts)

# Plotting the resulting temperature evolution as heatmap
using Plots
xgrid = L/(2Nx) : L/Nx : L;
tgrid = 0 : ts : tspan[2];
heatmap(tgrid, xgrid, sol2[:,:], xlabel="Time in [s]", ylabel="Length in [m]", title="Cooling Process of 1D Rod", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)