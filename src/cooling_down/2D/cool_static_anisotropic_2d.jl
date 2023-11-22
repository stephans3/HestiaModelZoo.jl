# # Cooling process of 1D rod 
# In this example we simulate the cooling process of a 2D plate with static (temperature-independent) and anisotropic material properties. 
using Hestia 

# Anisotropic thermal conductivity: λ=[λx 0; 0 λy]
λx = 10.0    
λy = 100.0   

# Mass density
ρ = 7800.0  

# Specific heat capacity
c = 480.0   

# Length and width of plate
L, W = 0.2, 0.2

# Discretized points
Nx, Ny = 40, 40     
Ntotal = Nx*Ny

# Material property with anisotropic and temperature-dependend (dynamic) behavior
property = StaticAnisotropic(λx, λy, ρ, c)
plate    = HeatPlate(L, W, Nx, Ny)


# Emission / Stefan-Boltzmann Boundary Condition
# - heat transfer coefficient: k=10
# - heat radiation with emssivity ϵ=0.6
# - ambient temperature  Θamb=300 Kelvin
emission = Emission(10.0, 0.6, 300.0) 

# The emission is assumed on each boundary side of the 2D plate.
boundary = Boundary(plate)
setEmission!(boundary, emission, :west)
setEmission!(boundary, emission, :east)
setEmission!(boundary, emission, :south)
setEmission!(boundary, emission, :north)

# The discrete heat conduction model is defined as right-hand side function to be solved with a numerical integrator.
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, plate, property, boundary)
end


# Initialize the 2D plate with  temperature for all cells
hot_square = 300*ones(Nx,Ny)
hot_square[10:Nx-10,10:Ny-10] .= 600.0

# Initial conditions of ODE
θinit = reshape(hot_square,Ntotal) # θ₀*ones(Ntotal)
tspan = (0.0, 2000.0)

# Save temperature at
ts = 1.0;


# Numerical integrator: (forward) Euler method
# If you use the (forward) Euler method you need to check the numerical stability.
# See: https://en.wikipedia.org/wiki/Euler_method
# See: https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis
Δx = plate.sampling[1]
Δy = plate.sampling[2]
αx = λx/(ρ *c)
αy = λy/(ρ *c)

# Sampling time
Δt    = 0.2             

# The chosen sampling time Δt needs to be smaller than the maximum possible sampling time. 
# If the chosen sampling time equals the maximum possible sampling time than the numerical integration works close to its stability limit.
max_dt = 0.5*inv(αx/(Δx^2) + αy/(Δy^2))

if Δt > max_dt
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

# Numerical integration method: forward Euler method
import OrdinaryDiffEq
alg1 = OrdinaryDiffEq.Euler();
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol1 = OrdinaryDiffEq.solve(prob,alg1, dt=Δt, saveat=ts)

using Plots
xgrid = L/(2Nx) : L/Nx : L;
ygrid = W/(2Ny) : W/Ny : W;
heatmap(xgrid, ygrid, reshape(sol1[end], Nx, Ny)', ylabel="Width in [m]", xlabel="Length in [m]", title="Final temperature distribution", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)

# Better use Runge-Kutta solvers which are able to handle stiff equations like KenCarp5.
# See: https://en.wikipedia.org/wiki/Stiff_equation
alg2 = OrdinaryDiffEq.KenCarp5()
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol2 = OrdinaryDiffEq.solve(prob, alg2, saveat=ts)
heatmap(ygrid, xgrid, reshape(sol2[end], Nx, Ny), ylabel="Width in [m]", xlabel="Length in [m]", title="Final temperature distribution", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)

# Error between Euler and KenCarp5 solution
err = sol2  - sol1;
heatmap(xgrid, ygrid, reshape(err[end], Nx, Ny)', ylabel="Width in [m]", xlabel="Length in [m]", title="Final temperature distribution", colorbar_title=" \nTemperature in [K]" , right_margin = 3Plots.mm)

