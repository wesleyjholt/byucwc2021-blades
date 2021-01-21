# This code is an example of simulating the NREL 5 MW wind turbine.

# import packages
using CCBlade
using PyPlot

# define rotor
Rhub = 1.5
Rtip = 63.0
B = 3
precone = 2.5*pi/180

rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)

# define the radial section properties
r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
    28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
    56.1667, 58.9000, 61.6333]
chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
    3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
    6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

# Define airfoils.  In this case we have 8 different airfoils that we load into an array.
# These airfoils are defined in files.
aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# create airfoil array
airfoils = aftypes[af_idx]

# define sections
sections = Section.(r, chord, theta, airfoils)

# operating point for the turbine
yaw = 0.0*pi/180
tilt = 5.0*pi/180
hubHt = 90.0
shearExp = 0.2

Vinf = 10.0
tsr = 7.55
rotorR = Rtip*cos(precone)
Omega = Vinf*tsr/rotorR
pitch = 0.0
azimuth = 0.0*pi/180
rho = 1.225

op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)

# solve
out = solve.(Ref(rotor), sections, op)

# plot distributed loads
figure()
plot(r/Rtip, out.Np/1e3)
plot(r/Rtip, out.Tp/1e3)
xlabel("r/Rtip")
ylabel("distributed loads (kN/m)")
legend(["flapwise", "lead-lag"])
savefig("example-figures/guided-example-turbine-operation-distributed-loads.pdf")

# integrate the loads to get thrust and torque (at azimuth=0)
T, Q = thrusttorque(rotor, sections, out)

# integrate the loads to get thrust and torque (average of four different azimuth angles)
azangles = pi/180*[0.0, 90.0, 180.0, 270.0]
ops = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho)
outs = solve.(Ref(rotor), sections, ops)

T, Q = thrusttorque(rotor, sections, outs)

# create the nondimensional power curve for the turbine (power coefficient vs tip-speed-ratio)
ntsr = 20  # number of tip-speed ratios
tsrvec = range(2, 15, length=ntsr)
cpvec = zeros(ntsr)  # initialize arrays
ctvec = zeros(ntsr)

azangles = pi/180*[0.0, 90.0, 180.0, 270.0]

for i = 1:ntsr
    local Omega = Vinf*tsrvec[i]/rotorR

    local ops = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho)
    local outs = solve.(Ref(rotor), sections, ops)
    local T, Q = thrusttorque(rotor, sections, outs)

    cpvec[i], ctvec[i], _ = nondim(T, Q, Vinf, Omega, rho, rotor, "windturbine")
end

figure()
plot(tsrvec, cpvec)
plot(tsrvec, ctvec)
xlabel("tip speed ratio")
legend([L"C_P", L"C_T"])
savefig("example-figures/guided-example-turbine-operation-nondimpowercurve.pdf")