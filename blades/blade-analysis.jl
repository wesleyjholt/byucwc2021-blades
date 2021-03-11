# This file sets up the blade analysis.

cd("blades")

# import packages
using CCBlade
using PyPlot
using DelimitedFiles
using Statistics

# specify which blade to analyze
# blade = "testblades01"
blade = "pennstate2019"

# get QBlade data
qblade_data = readdlm("input-files/testblades01-nondimpowercurve-8mps-qblade.txt", skipstart=3)
# qblade_data = readdlm("input-files/pennstate2019-nondimpowercurve-qblade.txt", skipstart=3)
qblade_tsrvec = qblade_data[:,1]
qblade_cpvec = qblade_data[:,2]

# set the parameters for the specified blade
include(string("params/", blade, ".jl"))

# ---- CREATE POWER CURVES ----
# set inflowing air parameters
Vinf_vec = range(5.0, stop=10.0, length=5)
rho = 1.225
mu = 1.647e-5
hubHt = .75
shearExp = 0.0 #0.0

# set turbine operating state
rotorR = Rtip*cos(precone)
pitch = 0.0*pi/180
yaw = 0.0*pi/180
tilt = 0.0*pi/180 #0.0*pi/180

# create vector of tip-speed ratios to analyze
ntsr = 20  # number of tip-speed ratios
tsrvec = range(.5, 7, length=ntsr) #range(.5, 5, length=ntsr)

# initialize arrays to hold values for power, thrust, and torque coefficients
cpvec = zeros(ntsr, length(Vinf_vec))
ctvec = zeros(ntsr, length(Vinf_vec))
cqvec = zeros(ntsr, length(Vinf_vec))

# create vector of blade azimuth angles at which to analyze turbine performance
azangles = pi/180*[0.0] #, 90.0, 180.0, 270.0]


figure()





figure()

for j = 1:length(Vinf_vec)

    # perform analysis for each tip-speed ratio
    for i = 1:ntsr

        local Omega = Vinf*tsrvec[i]/rotorR
        local ops = windturbine_op.(Vinf_vec[j], Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho, mu)
        local outs = solve.(Ref(rotor), sections, ops)
        local T, Q = thrusttorque(rotor, sections, outs)

        cpvec[i,j], ctvec[i,j], cqvec[i,j] = nondim(T, Q, Vinf_vec[j], Omega, rho, rotor, "windturbine")
    end

    # plot nondimensional power curve (power coefficient vs. tip-speed ratio)
    plot(tsrvec, cpvec[:,j])
    # plot(tsrvec, ctvec)
    # plot(qblade_tsrvec, qblade_cpvec)
    xlabel("tip speed ratio")
    ylabel(L"C_p")
    savefig(string("figures/", blade, "-nondimpowercurve.pdf"))

    maxcp, maxcp_ind = findmax(cpvec)
    println("optimal tip-speed ratio is: ", tsrvec[maxcp_ind])
    println("max cp: ", maxcp)

end

legend(string.(pitch_vec) .* " m/s wind speed")
savefig(string("figures/", blade, "-nondimpowercurve.pdf"))


cp_vec_mean = mean(cp_vec, axis=1)
plot(tsrvec, cp_vec_mean)
title("Non-dimensional Power Curve")
xlabel("tip speed ratio")
ylabel(L"C_p")
savefig(string("figures/", blade, "-nondimpowercurve-avg.pdf"))


# Vinf = Vinf_vec[1]

# run analysis at wind speed 8 m/s, tsr 3.0
azangles_1 = pi/180*[0.0, 90.0, 180.0, 270.0]
tsr_1 = 3.0
Omega_1 = Vinf*tsr_1/rotorR
# op_1 = windturbine_op(Vinf, Omega_1, pitch, r, precone, yaw, tilt, azangles_1', hubHt, shearExp, rho)
op_1 = simple_op.(Vinf, Omega_1, r, rho, mu=1.81e-5)    # constant inflow across rotor
out_1 = solve.(Ref(rotor), sections, op_1)
T_1, Q_1 = thrusttorque(rotor, sections, out_1)

cp_1, ct_1, cq_1 = nondim(T_1, Q_1, Vinf, Omega_1, rho, rotor, "windturbine")

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

# plot nondimensional power curve (power coefficient vs. tip-speed ratio)
figure()
plot(r/Rtip, out_1.Np)
plot(r/Rtip, out_1.Tp)
xlabel("r/Rtip")
legend([L"N_P"*" - CCBlade", L"T_P"*" - CCBlade"])
savefig(string("figures/", blade, "-normaltangentialforces.pdf"))
