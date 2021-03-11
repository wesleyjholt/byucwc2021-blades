#=
AUTHOR(S): Wesley Holt
DESCRIPTION: Contains functions for creating the dimensional power curve (Power vs. Wind velocity) and 
            nondimensional power curve (Power coefficient vs. tip-speed ratio) for a given set of blades.
NOTES: Also has a convenience function for interpolating chord and twist values from a default blade shape file.

(future contributors: feel free to add your name to "AUTHOR(S)")
=#

using CCBlade
using Statistics
using FLOWMath
using PyPlot
using DelimitedFiles
using FillArrays

# ================= PLOT DIMENSIONAL POWER CURVE (P VS. V) ===================
function plot_powercurve(blade, tsr, Vinf_vec; opt=true, file_extension=".pdf", transparent=false, chord=[0.0], theta=[0.0])

    # ----- SET TURBINE PARAMETERS -----

    # set rotor and blade parameters
    if opt
        include("params/" * blade * "-opt.jl")
    else
        include("params/" * blade * ".jl")
    end

    # define rotor and sections
    rotor = Rotor(Rhub, Rtip, B; turbine, precone)
    sections = Section.(r, chord, theta, airfoils)


    # ----- SET OPERATING CONDITIONS -----

    # set inflowing air parameters
    rho = 1.225
    mu = 1.647e-5
    hubHt = .75
    shearExp = 0.0

    # set turbine operating state
    rotorR = Rtip*cos(precone)
    pitch = 0.0*pi/180
    yaw = 0.0*pi/180
    tilt = 0.0*pi/180

    # get number of wind velocities to analyze
    nVinf = length(Vinf_vec)

    # initialize arrays to hold values for power
    P_vec = zeros(nVinf)

    # create vector of blade azimuth angles at which to analyze turbine performance
    azangles = pi/180*[0.0]


    # ----- GENERATE THE DIMENSIONAL POWER CURVE(S) -----

    figure()
    # iterate through each wind velocity
    for j = 1:nVinf
        # perform BEM analysis for each condition
        local Omega = Vinf_vec[j]*tsr/rotorR
        local ops = windturbine_op.(Vinf_vec[j], Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho, mu)
        local outs = solve.(Ref(rotor), sections, ops)
        local T, Q = thrusttorque(rotor, sections, outs)
        P_vec[j] = Q*Omega
    end

    # plot power curve (power vs. wind velocity)
    plot(Vinf_vec, P_vec, linewidth=3)
    title("Turbine Power Curve", fontsize=18)
    xlabel("Wind Velocity (m/s)", fontsize=14)
    ylabel("Power (W)", fontsize=14)
    if opt
        savefig(string("figures/", blade, "/" * blade * "-opt-powercurve" * file_extension), transparent=transparent, dpi=600)
    else
        savefig(string("figures/", blade, "/" * blade * "-orig-powercurve" * file_extension), transparent=transparent, dpi=600)
    end

end

# ============== PLOT NON-DIMENSIONAL POWER CURVE (C_P VS. TSR) ==============
function plot_nondimpowercurve(blades, tsr_vec, Vinf_vec=[8.0]; file_extension=".pdf", transparent=false, smooth=false)

    # if only a single string was passed in for the blades, convert the string into a vector with one string element
    if typeof(blades) == String
        blades = [blades]
    end

    # get number of blades, tip-speed ratios, and wind velocities to analyze
    nblades = length(blades)
    ntsr = length(tsr_vec)
    nVinf = length(Vinf_vec)

    # initialize vectors to store the data for all the blades
    tsr_vec_fine = collect(range(tsr_vec[1], stop=tsr_vec[end], length=100))
    cp_vec_all = fill(zeros(2),2)
    cp_vec_fine_all = fill(zeros(1),2)

    # iterate through each blade
    for b = 1:nblades
        blade = blades[b]

        # ----- SET TURBINE PARAMETERS -----
        # set rotor and blade parameters
        include("params/" * blade * ".jl")

        # define rotor and sections
        rotor = Rotor(Rhub, Rtip, B; turbine, precone)
        sections = Section.(r, chord, theta, airfoils)


        # ----- SET OPERATING CONDITIONS -----

        # set inflowing air parameters
        rho = 1.225
        mu = 1.647e-5
        hubHt = .75
        shearExp = 0.0

        # set turbine operating state
        rotorR = Rtip*cos(precone)
        pitch = 0.0*pi/180
        yaw = 0.0*pi/180
        tilt = 0.0*pi/180

        # initialize arrays to hold values for power, thrust, and torque coefficients
        cp_vec = zeros(ntsr, nVinf)
        ct_vec = zeros(ntsr, nVinf)
        cq_vec = zeros(ntsr, nVinf)

        # create vector of blade azimuth angles at which to analyze turbine performance
        azangles = pi/180*[0.0]


        # ----- PLOT NONDIMENSIONAL POWER CURVE(S) FOR EACH WIND SPEED -----

        figure()
        # iterate through each wind velocity
        for j = 1:nVinf
            # iterate through each tip-speed ratio
            for i = 1:ntsr
                # perform BEM analysis for each condition
                local Omega = Vinf_vec[j]*tsr_vec[i]/rotorR
                local ops = windturbine_op.(Vinf_vec[j], Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho, mu)
                local outs = solve.(Ref(rotor), sections, ops)
                local T, Q = thrusttorque(rotor, sections, outs)

                # nondimensionalize the power, thrust, and torque
                cp_vec[i,j], ct_vec[i,j], cq_vec[i,j] = nondim(T, Q, Vinf_vec[j], Omega, rho, rotor, "windturbine")
            end

            # plot nondimensional power curve (power coefficient vs. tip-speed ratio)
            plot(tsr_vec, cp_vec[:,j])
            title("Non-dimensional Power Curves\nat Various Wind Speeds: " * blade, fontsize=14)
            xlabel("Tip-speed Ratio")
            ylabel(L"C_p")
        end

        legend(string.(Vinf_vec) .* " m/s", title="Wind Speeds")
        try mkdir("figures/" * blade); catch; end
        savefig("figures/" * blade * "/" * blade * "-nondimpowercurve-windspeeds" * file_extension, transparent=transparent, dpi=600)

        # ---- PLOT A NONDIMENSIONAL POWER CURVE (AVERAGING ALL WIND SPEEDS) ----

        figure()
        # calculate average power coefficient values
        cp_vec_mean = mean(cp_vec, dims=2)[:,1]
        if smooth
            # plot smoothed power curve (for the poster)
            cp_vec_fine = akima(tsr_vec, cp_vec_mean, tsr_vec_fine)
            plot(tsr_vec_fine, cp_vec_fine, linewidth=3)
        else
            # plot nonsmoothed power curve
            plot(tsr_vec, cp_vec_mean, linewidth=3)
        end
        title("Non-dimensional Power Curve: " * blade, fontsize=14)
        xlabel("Tip-speed Ratio", fontsize=12)
        ylabel(L"C_p", fontsize=12)
        savefig("figures/" * blade * "/" * blade * "-nondimpowercurve" * file_extension, transparent=transparent, dpi=600)

        # save the power coefficient values into vector holding data for all blades
        cp_vec_all[b] = cp_vec_mean
        if smooth
            cp_vec_fine_all[b] = cp_vec_fine
        end

        # save power coefficient values into a text file
        try mkdir("test-data/" * blade); catch; end
        open("test-data/" * blade * "/" * blade * "-nondimpowercurve.txt", "w") do io
            write(io, "# tip-speed ratio, power coefficient\n")
            writedlm(io, [tsr_vec cp_vec_mean])
        end

    end
    
    # ------- PLOT NONDIMENSIONAL POWER CURVE FOR EACH BLADE -------

    if nblades > 1
        figure()
        blades_string = ""
        for b = 1:nblades
            if b == 1
                blades_string *= blades[b]
            else
                blades_string *= "_" * blades[b]
            end
            if smooth
                plot(tsr_vec_fine, cp_vec_fine_all[b], linewidth=3, alpha=0.7)
            else
                plot(tsr_vec, cp_vec_all[b], linewidth=3, alpha=0.7)
            end
        end
        title("Non-dimensional Power Curve Comparison", fontsize=14)
        xlabel("Tip-speed Ratio", fontsize=12)
        ylabel(L"C_p", fontsize=12)
        legend(blades)
        savefig("figures/blade-comparison/nondimpowercurve-" * blades_string * file_extension, transparent=transparent, dpi=600)
    end

end

# ===================== INTERPOLATE CHORD AND TWIST ==========================
function interp_chord_twist(r, default_file; default_twist_in_radians=true)
    #= 
    Use this function to get interpolated chord and twist values if the the params file doesn't specify it.
    This is meant to be a quick-and-dirty way to get some values for chord and twist from some existing file.
    Final values for chord and twist should be determined from an optimization.
    =#

    # get default chord and twist values
    default_data = readdlm(default_file, skipstart=1)
    chord_default = default_data[:,2]
    twist_default = default_data[:,3]
    if !default_twist_in_radians
        twist_default .*= pi/180
    end

    # scale the chord if it is way smaller than the length of the blade
    if Rtip > 1e2*maximum(chord_default)
        chord_default *= Rtip/(10*maximum(chord_default))
    end

    # interpolate chord and twist values for each radial section
    r_default = range(r[1], stop=r[end], length=length(chord_default))
    chord = linear(r_default, chord_default, r)
    twist = linear(r_default, twist_default, r)

    return chord, twist
end