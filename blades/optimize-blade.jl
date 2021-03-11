#=
AUTHOR(S): Dr. Andrew Ning (modified by: Wesley Holt)
DESCRIPTION: Blade optimization script. (for BYU's blade design in the 2021 Collegiate Wind Competition)
NOTES:

(future contributors: feel free to add your name to "AUTHOR(S)")
=#

using CCBlade
using FLOWMath
using SpecialFunctions: gamma
import LinearAlgebra: dot
using Snopt
import ForwardDiff
import ReverseDiff
import SparseArrays, SparseDiffTools
import FiniteDiff
import NLsolve
using DelimitedFiles

# =================================== I/O FUNCTIONS =======================================

function write_blade_to_params_file(blade, chord, theta, cspline=[], tspline=[], tsr=[], pitch=[]; append_to_file_name="-")

    # write blade parameters to a .jl file
    open("params/" * blade * append_to_file_name * ".jl", "w") do io
        write(io, "# ================================================\n")
        write(io, "# ============ MACHINE GENERATED FILE ============\n")
        write(io, "# ================================================\n\n\n")
        write(io, read("params/nochordtwist/" * blade * "-nochordtwist.jl"))
        write(io, "\n\n# chord and twist\nchord = $chord\ntheta = $theta")
        write(io, "\n\n# chord and twist spline parameters\ncspline = $cspline\ntspline = $tspline")
        write(io, "\n\n# optimal tip-speed ratio\ntsr = $tsr")
        write(io, "\n\n# optimal pitch angles\npitch = $pitch")
    end

end


# ============= FUNCTIONS REQUIRED FOR OPTIMIZATION/DERIVATIVE CALCULATION ================

function residual2d(as, aps, rotor, section, op)

    phi = atan((1 + as)*op.Vx/((1 - aps)*op.Vy))
    _, outputs = CCBlade.residual(phi, rotor, section, op)
    
    R = [-outputs.a - as; 10*(-outputs.ap - aps)]

    return R, outputs

end

function solve2d(rotor, section, op)

    # ----- solve residual function ------    

    # wrapper to residual function
    R(x) = residual2d(x[1], x[2], rotor, section, op)[1]

    a0 = -1.0/3
    ap0 = -0.1
    x0 = [a0; ap0]  

  
    results = NLsolve.nlsolve(R, x0, method=:trust_region, ftol=1e-14, autodiff=:forward)
    if !NLsolve.converged(results)
        println("not converged")
        println(rotor)
        println(section)
        println(op)
    end
    x = results.zero
    
    _, outputs = residual2d(x[1], x[2], rotor, section, op)

    return outputs
end

function residualadjoint(R, ymin, ymax, x)
    y, _ = FLOWMath.brent(y -> R(x, y), ymin, ymax)
    return y
end
 
function residualadjoint(R, ymin, ymax, d::Vector{D}) where D <: ForwardDiff.Dual

    # solve root
    x = ForwardDiff.value.(d)
    y = residualadjoint(R, ForwardDiff.value(ymin), ForwardDiff.value(ymax), x)

    # compute derivatives
    wrap(xy) = R(xy[1:end-1], xy[end])
    g = Vector{Float64}(undef, length(x) + 1)
    ForwardDiff.gradient!(g, wrap, [x; y])
    drdx = g[1:end-1]
    drdy = g[end]
    dydx = -drdx/drdy

    # assemble dual number
    b_in = zip(collect.(ForwardDiff.partials.(d))...)
    b_arr = map(x->dot(dydx, x), b_in)
    p = ForwardDiff.Partials((b_arr...,))
    return D(y, p)
end

function solveadjoint(rotor, section, op)

    # error handling
    if typeof(section) <: Vector
        error("You passed in an vector for section, but this funciton does not accept an vector.\nProbably you intended to use broadcasting (notice the dot): solve.(Ref(rotor), sections, ops)")
    end

    # check if we are at hub/tip
    if isapprox(section.r, rotor.Rhub, atol=1e-6) || isapprox(section.r, rotor.Rtip, atol=1e-6)
        return Outputs()  # no loads at hub/tip
    end

    # parameters
    npts = 10  # number of discretization points to find bracket in residual solve

    # unpack
    Vx = op.Vx
    Vy = op.Vy
    theta = section.theta + op.pitch

    # ---- determine quadrants based on case -----
    Vx_is_zero = isapprox(Vx, 0.0, atol=1e-6)
    Vy_is_zero = isapprox(Vy, 0.0, atol=1e-6)

    # quadrants
    epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    if Vx_is_zero && Vy_is_zero
        return Outputs()

    elseif Vx_is_zero

        startfrom90 = false  # start bracket at 0 deg.

        if Vy > 0 && theta > 0
            order = (q1, q2)
        elseif Vy > 0 && theta < 0
            order = (q2, q1)
        elseif Vy < 0 && theta > 0
            order = (q3, q4)
        else  # Vy < 0 && theta < 0
            order = (q4, q3)
        end

    elseif Vy_is_zero

        startfrom90 = true  # start bracket search from 90 deg

        if Vx > 0 && abs(theta) < pi/2
            order = (q1, q3)
        elseif Vx < 0 && abs(theta) < pi/2
            order = (q2, q4)
        elseif Vx > 0 && abs(theta) > pi/2
            order = (q3, q1)
        else  # Vx < 0 && abs(theta) > pi/2
            order = (q4, q2)
        end

    else  # normal case

        startfrom90 = false

        if Vx > 0 && Vy > 0
            order = (q1, q2, q3, q4)
        elseif Vx < 0 && Vy > 0
            order = (q2, q1, q4, q3)
        elseif Vx > 0 && Vy < 0
            order = (q3, q4, q1, q2)
        else  # Vx[i] < 0 && Vy[i] < 0
            order = (q4, q3, q2, q1)
        end

    end

        

    # ----- solve residual function ------

    # # wrapper to residual function to accomodate format required by fzero
    R(phi) = CCBlade.residual(phi, rotor, section, op)[1]

    function R2(x, y)
        section2 = Section(section.r, x[2], x[3], section.af)
        op2 = OperatingPoint(x[4], x[5], op.rho, x[1], op.mu, op.asound)

        return CCBlade.residual(y, rotor, section2, op2)[1]
    end
    x2 = [op.pitch, section.chord, section.theta, op.Vx, op.Vy]

    success = false
    for j = 1:length(order)  # quadrant orders.  In most cases it should find root in first quadrant searched.
        phimin, phimax = order[j]

        # check to see if it would be faster to reverse the bracket search direction
        backwardsearch = false
        if !startfrom90
            if phimin == -pi/2 || phimax == -pi/2  # q2 or q4
                backwardsearch = true
            end
        else
            if phimax == pi/2  # q1
                backwardsearch = true
            end
        end
        
        # force to dual numbers if necessary
        phimin = phimin*one(section.chord)
        phimax = phimax*one(section.chord)

        # find bracket
        success, phiL, phiU = CCBlade.firstbracket(R, phimin, phimax, npts, backwardsearch)

        # once bracket is found, solve root finding problem and compute loads
        if success
            # phistar, _ = FLOWMath.brent(R, phiL, phiU)
            phistar = residualadjoint(R2, phiL, phiU, x2)
            _, outputs = CCBlade.residual(phistar, rotor, section, op)
            return outputs
        end    
    end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return Outputs()
end


# =========================== FUNCTIONS FOR TURBINE ANALYSIS ==============================

function shape(cspline, tspline, r, Rhub, Rtip)
    #= 
    ---- DESCRIPTION ----
    The function interpolates the chord and twist splines.
    ---- INPUTS ---- 
    cspline (array): chord spline parameters
    tspline (array): twist spline parameters
    r (array): radial values at which blade sections are defined
    Rhub (float): hub radius
    Rtip (float): hub tip
    ---- OUTPUTS ----
    chord (array): the chord at each section of the blade
    theta (array): the twist angle at each section of the blade
    =#

    # interpolate chord lengths using akima spline
    rspline = [0.00, 0.25, 0.5, 0.75, 1.0]*(Rtip-Rhub) .+ Rhub 
    chord = akima(rspline, cspline, r) 
    
    # interpolate twist angles using akima spline
    rspline2 = [0.11111, 0.4, 0.7, 1.0]*(Rtip-Rhub) .+ Rhub 
    theta = [[tspline[1], tspline[1]]; akima(rspline2, tspline, r[3:end])]  # fixed root twist b/c it is a cylinder

    return chord, theta
end


function turbineopt(cspline, tspline, tsr, pitch, method, blade)
    #= 
    ---- DESCRIPTION ----
    This is the turbine analysis function.
    ---- INPUTS ---- 
    cspline (array): chord spline parameters
    tspline (array): twist spline parameters
    tsr (float): tip-speed ratio
    pitch (array): pitch angle at each wind speed
    method (string): gradient calculation/implementation method
    blade (string): name of the blade file (in the "params" directory)
    ---- OUTPUTS ----
    AEP, V, T, P, Q, Np
    AEP: annual energy production
    V: wind velocities
    T: thrust values at each wind velocity
    P: power values at each wind velocity
    Q: torque values at each wind velocity
    Np: normal force per unit length (N/m) near the rated speed
    =#

    # println("cspline = ", cspline)
    # println("tspline = ", tspline)
    # println("tsr = ", tsr)
    # println("pitch = ", pitch)

    # ----- SET TURBINE PARAMETERS -----

    # define pitch to start at zero, add root pitch to twist instead.
    p0 = pitch[1]
    pitch .-= p0
    tspline .+= p0

    # set rotor and blade parameters
    include("params/nochordtwist/" * blade * "-nochordtwist.jl")

    # get chord and twist at the defined radial locations
    chord, theta = shape(cspline, tspline, r, Rhub, Rtip)   # get blade geometric parameters for aerodynamic analysis

    # rotational correction
    # du = DuSeligEggers()  # Du-Selig correction for lift, Eggers correction for drag
    
    # define rotor and sections
    # rotor = Rotor(Rhub, Rtip, B, rotation=du, precone=precone, turbine=true)
    rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)
    sections = Section.(r, chord, theta, airfoils)


    # ----- SET OPERATING CONDITIONS -----

    # set operating point for the turbine
    yaw = 0.0*pi/180    # yaw angle
    tilt = 0.0*pi/180   # tilt angle
    hubHt = 90.0        # hub height
    shearExp = 0.0      # wind shear exponent

    rotorR = Rtip*cos(precone)  # rotor radius
    azimuth = 0.0*pi/180    # azimuth angles at which the blade is on
    rho = 1.225  # air density #*one(pitch[1])

    V = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]   # wind velocities tested for the power curve performance task in the 2021 CWC
    nV = length(V)  # number of wind velocities tested

    # Omega_min = 0.0
    Omega_max = 12.0*pi/30.0    # blade rotation speed (radians/second)
    Omega = min.(collect(V).*tsr/rotorR, Omega_max)  # collect needed for ForwardDiff
    P = similar(V, eltype(chord[1]))
    Q = similar(V, eltype(chord[1]))
    T = similar(V, eltype(chord[1]))


    # ----- RUN TURBINE ANALYSIS -----

    # set the solver type based on the derivative calculation method specified
    if method == "traditional"
        solver = solve2d
    elseif method == "adjoint"
        solver = solveadjoint
    else
        solver = solve
    end

    # iterate through wind velocities
    for i = 1:nV
        # create operation state object
        ops = windturbine_op.(V[i], Omega[i], pitch[i], r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
        # run BEM analysis
        outputs = solver.(Ref(rotor), sections, ops)
        # get thrust, torque, and power values
        T[i], Q[i] = thrusttorque(rotor, sections, outputs)
        P[i] = Q[i]*Omega[i]
    end

    # find the max normal force on the blade
    i = nV      # this should be the velocity at which the blade will experience the largest normal force
    ops = windturbine_op.(V[i], Omega[i], zero(rho), r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
    outputs = solver.(Ref(rotor), sections, ops)
    Np = outputs.Np


    # ----- GET ANNUAL ENERGY PRODUCTION (AEP) -----

    # define wind speed weightings for the power curve task (set by the 2021 CWC)
    w = [0.7, 0.8, 0.8, 0.7, 0.4, 0.3, 0.1]

    # calculate AEP (weighted sum of the powers at each wind velocity)
    AEP = sum(P .* w)

    # println("AEP = ", AEP)
    # println("V = ", V)
    # println("T = ", T)
    # println("P = ", P)
    # println("Q = ", Q)
    # println("Np = ", Np)

    return AEP, V, T, P, Q, Np
end


# ========================= FUNCTION THAT RUNS THE OPTIMIZATION ===========================

function run_opt(method, nP, blade, x0, lb, ub)
    #= 
    ---- DESCRIPTION ----
    Sets up and runs the optimization
    ---- INPUTS ---- 
    method: the method used for gradient calculation/implementation ("fd", "traditional", "dense", "adjoint", or "sparse")
    nP: number of pitch angles
    ---- OUTPUTS ----
    xopt: optimal values for the design variables
    fopt: optimal objective value
    =#

    # # bounds for the design variables
    # lb = [0.001*ones(5); -5*pi/180*ones(4); 1.0; zeros(nP)]
    # ub = [0.2*ones(5); 20*pi/180*ones(4); 15.0; 30.0*pi/180*ones(nP)]
    # lb[1] = 0.01 #3.542  # lower bound on root chord, for pitch bearing

    # rated power
    Prated = 5e6

    # set the objective scale factor
    obj_scale = 1e7

    function wrappersparse!(y, x)
        #= 
        ---- DESCRIPTION ----
        Gets the objective and constraint values by calling the turbine analysis function
        ---- INPUTS ---- 
        x: vector of design variable values
        y: vector that gets populated with the constraint values (the function mutates this argument)
        ---- OUTPUTS ----
        scaled_AEP: annual energy production of the turbine, scaled down to be a value appropriate for optimization 
        =#

        # separate design variables vector into separate vectors for chord, twist, tip-speed ratio, and pitch angles
        cspline2 = x[1:5]
        tspline2 = x[6:9]
        tsr2 = x[10]
        idx = 10
        pitch2 = x[idx+1:idx+nP]

        # calculate objective (AEP) and other values needed for contraints
        AEP, V, T, P, Q, Np = turbineopt(cspline2, tspline2, tsr2, pitch2, method, blade)
        scaled_AEP = -AEP/obj_scale   # scale the objective here so that the objective derivatives are between 0 and 1

        # calculate constraints (listed below)
        # -> power is below the rated power at each wind speed
        # -> torque is below the maximum torque
        # -> pitch angles are increasing
        # -> normal forces per unit length at each radial section are below 6500 N/m at the rated speed

        # Qmax = Prated/Omega_max
        Tmax = 600e3
        y .= [P/Prated .- 1; T/Tmax .- 1; pitch2[1:end-1] .- pitch2[2:end]; Np/6500 .- 1]  #; Omega[1:end-1] .- Omega[2:end]]

        return scaled_AEP  # maximize AEP
    end

    function wrapper!(y, x)
        #= 
        ---- DESCRIPTION ----
        Calls another wrapper function to get objective and constraint values
        ---- INPUTS ---- 
        x: vector of design variable values
        y: vector that gets populated with the objective value as the first element and the constraint values as the rest (the function mutates this argument)
        =#

        # separate the constraint values from the objective value
        yv = @view y[2:end]
        # call another wrapper function to get the objective and constraint values
        y[1] = wrappersparse!(yv, x)
    end

    # initialize arrays to store values needed for optimization
    nsections = 17                          # number of radial sections
    y = zeros(1 + 2*nP + nP-1 + nsections)  # objective and contraint values (all in one vector)
    config = ForwardDiff.JacobianConfig(wrapper!, y, x0)
    J = zeros(length(y), length(x0))        # Jacobian array
    dfdx = zeros(length(x0))                # objective gradients

    # initialize additional arrays needed for sparse derivatives method
    if method == "sparse"
        ysp = zeros(2*nP + nP-1 + nsections)   # objective and contraint values (all in one vector)
        Jspstart = zeros(length(ysp), length(x0))   # Jacobian array
        configsp = ForwardDiff.JacobianConfig(wrappersparse!, ysp, x0)
        ForwardDiff.jacobian!(Jspstart, wrappersparse!, ysp, x0, configsp)
        Jsp = SparseArrays.sparse(Jspstart)  
        colors = SparseDiffTools.matrix_colors(Jsp)
        cache = SparseDiffTools.ForwardColorJacCache(wrappersparse!, x0, dx=ysp, colorvec=colors, sparsity=Jsp)
        
        # setup for AEP derivative
        # Prated = 5e6
        # Vin = 3.0
        # Vout = 25.0
        # Vvec = range(Vin, Vout, length=nP)
        # k = 2.0  # weibull shape
        # Vbar = 6.0  # mean speed of weibull
        # A = Vbar / gamma(1.0 + 1.0/k)
        # cdf = 1.0 .- exp.(-(Vvec/A).^k)
    end


    function optimize(x)
        #= 
        ---- DESCRIPTION ----
        This is the function passed into the optimizer (SNOPT)
        ---- INPUTS ---- 
        x: vector of design variable values
        ---- OUTPUTS ----
        f: objective value
        g: constraint values
        dfdx: objective gradients
        dgdx: constraint gradients
        fail: output Boolean variable needed for SNOPT
        =#

        # call the wrapper to calculate the objective and constraint values
        wrapper!(y, x)
        f = y[1]
        g = y[2:end]

        # choose a method and calculate gradients
        if method == "fd" || method == "traditional"
            FiniteDiff.finite_difference_jacobian!(J, wrapper!, x)
            dfdx .= J[1, :]
            dgdx = J[2:end, :] 
        elseif method == "dense" || method == "adjoint"
            ForwardDiff.jacobian!(J, wrapper!, y, x, config)
            dfdx .= J[1, :]
            dgdx = J[2:end, :] 
        elseif method == "sparse"
            SparseDiffTools.forwarddiff_color_jacobian!(Jsp, wrappersparse!, x, cache)
            dgdx = Jsp
            dPdx = Prated*Jsp[1:nP, :]
            # dfdx .= 0.0
            # for i = 1:nP-1
            #     dfdx .+= (cdf[i+1]-cdf[i])*0.5*(dPdx[i, :] + dPdx[i+1, :]) *365*24/-1e10
            # end
            w = [0.7, 0.8, 0.8, 0.7, 0.4, 0.3, 0.1]
            # println("dfdx_zeros = ", dfdx)
            # println("dfdx = ", permutedims(permutedims(w)*dPdx))
            for i = 1:nP
                dfdx .+= w[i] .* 0.5*dPdx[i, :]/-obj_scale
            end
        end

        fail = false

        return f, g, dfdx, dgdx, fail
    end

    # set snopt options
    options = Dict{String, Any}()
    options["Derivative option"] = 1
    options["Verify level"] = 3
    
    # call snopt
    xopt, fopt, info = snopt(optimize, x0, lb, ub, options)
    println(info)
    
    return xopt, fopt
end


# =========================================================================================
# ================================= SET UP PARAMETERS =====================================
# =========================================================================================

# ================= This block contains all the user inputs. ==============================

# set which blade to optimize
blade = "NREL5MW"   # this should match the blade name given by "params/nochordtwist/BLADE_NAME-nochordtwist.jl"
# (ex: if the blade parameter file is params/nochordtwist/NREL5MW-nochordtwist.jl then set blade="NREL5MW")
# (the parameters file "BLADE_NAME-nochordtwist.jl" should not specify any chord or twist value for the blade)

# number of pitch angles
nP = 7      # 7 pitch angles for each of the 7 wind speeds for the rated power task of the 2021 CWC

# # define the starting point for the design variables
# cspline = [0.050974, 0.061282, 0.046799, 0.037665, 0.031924] # chord
# tspline = [42.8, 15.579, 8.861, 6.506] * pi/180   # twist
# tsr = 4.0   # tip-speed ratio
# pitch = 20*pi/180*ones(nP)  # pitch

cspline = [3.542, 4.6035, 3.748, 2.8255, 1.419] # chord
tspline = [13.308, 8.0, 4.1, 0.05]*pi/180   # twist
tsr = 7.55  # tip-speed ratio
pitch = 20*pi/180*ones(nP)  # pitch

# # define bounds for the design variables
# cspline_lb = 0.001*ones(5)
# cspline_ub = 0.2*ones(5)
# tspline_lb = -5*pi/180*ones(4)
# tspline_ub = 40*pi/180*ones(4)
# tsr_lb = 1.0
# tsr_ub = 15.0
# pitch_lb = zeros(nP)
# pitch_ub = 30.0*pi/180*ones(nP)

cspline_lb = 0.5*ones(5)
cspline_ub = 10*ones(5)
tspline_lb = -5*pi/180*ones(4)
tspline_ub = 20*pi/180*ones(4)
tsr_lb = 1.0
tsr_ub = 15.0
pitch_lb = zeros(nP)
pitch_ub = zeros(nP)#30.0*pi/180*ones(nP)
cspline_lb[1] = 3.542  # lower bound on root chord, for pitch bearing


# # bounds for the design variables
# lb = [0.001*ones(5); -5*pi/180*ones(4); 1.0; zeros(nP)]
# ub = [0.2*ones(5); 20*pi/180*ones(4); 15.0; 30.0*pi/180*ones(nP)]
# lb[1] = 0.01 #3.542  # lower bound on root chord, for pitch bearing

# rated power
Prated = 5e7

# =========================================================================================
# ========= (EVERYTHING BELOW THIS POINT NORMALLY DOES NOT NEED TO BE EDITED) =============
# =========================================================================================


# --------------- WRITE THE INITIAL BLADE PARAMETERS TO A PARAMS FILE ---------------------

# create new blade file names
blade_init = blade * "-init"
blade_opt = blade * "-opt"

# get initial blade parameters
include("params/nochordtwist/" * blade * "-nochordtwist.jl")
chord, theta = shape(cspline, tspline, r, Rhub, Rtip)

# write the initial blade parameters to a params file
write_blade_to_params_file(blade, chord, theta, cspline, tspline, tsr, pitch, append_to_file_name="-init")

# ------------------------------- GET THE INITIAL AEP -------------------------------------

# combine starting point values into a single vector
x0 = [cspline; tspline; tsr; pitch]

# get the initial AEP value
AEP_init, V_init, T_init, P_init, Q_init, Np_init = turbineopt(cspline, tspline, tsr, pitch, "dense", blade)
println("Initial AEP value = ", AEP_init)

# ------------------------------- RUN THE OPTIMIZATION ------------------------------------

# combine bounds for all design variables into two arrays: a lower bound array and an upper bound array
lb = [cspline_lb; tspline_lb; tsr_lb; pitch_lb]
ub = [cspline_ub; tspline_ub; tsr_ub; pitch_ub]

# call the optimizer wrapper function
xopt, fopt = run_opt("dense", nP, blade, x0, lb, ub) # options for derivative calculation (first argument), from fastest to slowest: sparse, adjoint, dense, fd, traditional 
                                                    # (note: "sparse" is not yet set up for CWC 2021 blade optimization)

# ------------------------------ GET THE OPTIMIZED AEP ------------------------------------    

# separate the design variables into appropriate vectors
cspline = xopt[1:5]
tspline = xopt[6:9]
tsr = xopt[10]
idx = 10
pitch = xopt[idx+1:idx+nP]

# get the final AEP value
AEP, V, T, P, Q, Np = turbineopt(cspline, tspline, tsr, pitch, "dense", blade)
println("Final AEP value = ", AEP)

# -------------- CREATE PLOTS TO VISUALIZE THE BLADE DESIGN AND PERFORMANCE ---------------

# import plotting package
using PyPlot

# create a directory to store all plots for the optimized blade (if needed)
try mkdir("figures/" * blade_opt); catch; end

# close any existing plots
close("all")

# power curve
figure()
plot(V, P/1e6)
title("Optimized " * blade * " Power Curve")
xlabel(L"V" * " (m/s)")
ylabel(L"P" * " (MW)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-powercurve.pdf")

# thrust curve
figure()
plot(V, T)
title("Optimized " * blade * " Thrust Curve")
xlabel(L"V" * " (m/s)")
ylabel(L"T" * " (N)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-thrustcurve.pdf")

# torque curve
figure()
plot(V, Q)
title("Optimized " * blade * " Torque Curve")
xlabel(L"V" * " (m/s)")
ylabel(L"Q" * " (N-m)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-thrustcurve.pdf")

# pitch
figure()
plot(V, pitch*180/pi)
title("Optimized " * blade * " Pitch")
xlabel(L"V" * " (m/s)")
ylabel("pitch (deg)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-optimal-pitch.pdf")

# get the chord and twist values at each radial section
chord, theta = shape(cspline, tspline, r, Rhub, Rtip)

# normal forces
figure()
plot(r, Np)
title("Optimized " * blade * " Normal Forces")
xlabel(L"r" * " (m)")
ylabel(L"N^\prime" * " (N/m)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-normalforces.pdf")

# chord
figure()
plot(r, chord)
title("Optimized " * blade * " Chord")
xlabel(L"r" * " (m)")
ylabel("chord (m)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-optimal-chord.pdf")

# twist
figure()
plot(r, theta*180/pi)
title("Optimized " * blade * " Twist")
xlabel(L"r" * " (m)")
ylabel("theta (deg)")
savefig("figures/" * blade_opt * "/" * blade_opt * "-optimal-twist.pdf")

# ---------------- WRITE THE FINAL BLADE PARAMETERS TO A PARAMS FILE ---------------------

# write the results to a params file
write_blade_to_params_file(blade, chord, theta, cspline, tspline, tsr, pitch, append_to_file_name="-opt")

# ----------------------- PLOT THE NONDIMENSIONAL POWER CURVE ----------------------------

# import the file containing power curve plotting functions
include("plot-powercurve.jl")

# specify the plotting range for tip-speed ratio and wind velocity
tsr_vec = range(1, stop=15.0, length=50)
Vinf_vec = range(1, stop=11.0, step=1.0)

# plot the nondimensional power curve
plot_nondimpowercurve(blade_opt, tsr_vec, Vinf_vec; file_extension=".pdf", transparent=false)
plot_nondimpowercurve(blade_opt, tsr_vec, Vinf_vec; file_extension=".png", transparent=false)




println("cspline = ", cspline)
println("tspline = ", tspline * 180/pi, " degrees")
println("tsr = ", tsr)
println("pitch = ", pitch)