using CCBlade
using PyPlot

Rtip = 10/2.0 * 0.0254  # inches to meters
Rhub = 0.10*Rtip
B = 2  # number of blades

# create Rotor object
rotor = Rotor(Rhub, Rtip, B)

# normalized radial position, chord, and twist
propgeom = [
0.15   0.130   32.76
0.20   0.149   37.19
0.25   0.173   33.54
0.30   0.189   29.25
0.35   0.197   25.64
0.40   0.201   22.54
0.45   0.200   20.27
0.50   0.194   18.46
0.55   0.186   17.05
0.60   0.174   15.97
0.65   0.160   14.87
0.70   0.145   14.09
0.75   0.128   13.39
0.80   0.112   12.84
0.85   0.096   12.25
0.90   0.081   11.37
0.95   0.061   10.19
1.00   0.041   8.99
]

# scale and organize radial position, chord, and twist into separate vectors
r = propgeom[:, 1] * Rtip
chord = propgeom[:, 2] * Rtip
theta = propgeom[:, 3] * pi/180

# import airfoil data
af = AlphaAF("example-airfoil-data/naca4412.dat")

# define sections
sections = Section.(r, chord, theta, Ref(af))

# define operating points
Vinf = 5.0
Omega = 5400*pi/30  # convert to rad/s
rho = 1.225
op = simple_op.(Vinf, Omega, r, rho)    # constant inflow across rotor

# BEM analysis
out = solve.(Ref(rotor), sections, op)

# display results of interest
# out.alpha*180/pi  # angles of attack in degrees at each radial station
# out.Np  # normal force per unit length (N/m)
# out.Tp  # tangential force per unit length (N/m)
# out.a   # axial induction factor
# out.ap  # tangential induction factor
# out.u   # axial induced velocity
# out.v   # tangential induced velocity
# out.phi # inflow angle
# out.W   # inflow velocity
# out.cl  # lift coefficient
# out.cd  # drag coefficient
# out.cn  # normal force coefficient
# out.ct  # tangential force coefficient
# out.F   # hub/tip loss correction
# out.G   # effective hub/tip loss correction for induced velocities: u = Vx * a * G, v = Vy * ap * G

figure()
plot(r/Rtip, out.Np)
plot(r/Rtip, out.Tp)
xlabel("r/Rtip")
ylabel("distributed loads (N/m)")
legend(["flapwise", "lead-lag"])
savefig("example-figures/quick-start-example-distributed-loads.pdf")

T, Q = thrusttorque(rotor, sections, out)
cp, ct, _ = nondim(T, Q, Vinf, Omega, rho, rotor, "windturbine")
println(cp)