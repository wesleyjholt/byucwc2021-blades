# This file sets up the blade analysis.

# import packages
using CCBlade
using PyPlot
using DelimitedFiles

# specify which blade to analyze
blade = "pennstate2019"

# set the parameters for the specified blade
include(string("params/", blade, ".jl"))

# ---- CREATE NON-DIMENSIONAL POWER CURVE ----
# set inflowing air parameters
Vinf = 10.0
rho = 1.225
hubHt = .75
shearExp = 0.0

# set turbine operating state
rotorR = Rtip*cos(precone)
pitch = 0.0
yaw = 0.0*pi/180
tilt = 0.0*pi/180

# create vector of tip-speed ratios to analyze
ntsr = 20  # number of tip-speed ratios
tsrvec = range(.5, 5, length=ntsr)

# initialize arrays to hold values for power, thrust, and torque coefficients
cpvec = zeros(ntsr)
ctvec = zeros(ntsr)
cqvec = zeros(ntsr)

# create vector of blade azimuth angles at which to analyze turbine performance
azangles = pi/180*[0.0, 90.0, 180.0, 270.0]

# function solve_

# perform analysis for each tip-speed ratio
for i = 1:ntsr

    local Omega = Vinf*tsrvec[i]/rotorR
    local ops = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho)
    local outs = solve.(Ref(rotor), sections, ops)
    local T, Q = thrusttorque(rotor, sections, outs)

    cpvec[i], ctvec[i], cqvec[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "windturbine")
end

# get QBlade data
qblade_data = readdlm("input-files/pennstate2019-nondimpowercurve-qblade.txt", skipstart=3)
qblade_tsrvec = qblade_data[:,12*2-1]
qblade_cpvec = qblade_data[:,12*2]

# plot nondimensional power curve (power coefficient vs. tip-speed ratio)
figure()
plot(tsrvec, cpvec)
# plot(tsrvec, ctvec)
plot(qblade_tsrvec, qblade_cpvec)
xlabel("tip speed ratio")
legend([L"C_P"*" - CCBlade", L"C_T"*" - QBlade"])
savefig(string("figures/", blade, "-nondimpowercurve.pdf"))

maxcp, maxcp_ind = findmax(cpvec)
println("optimal tip-speed ratio is: ", tsrvec[maxcp_ind])
println("max cp: ", maxcp)


# Vinf = 10.0
# tsr = 6.0
# rotorR = Rtip*cos(precone)
# Omega = Vinf*tsr/rotorR
# pitch = 0.0
# azimuth = 0.0*pi/180
# rho = 1.225

# op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)