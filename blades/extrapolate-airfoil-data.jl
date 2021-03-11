#=
AUTHOR(S): Wesley Holt
DESCRIPTION: Contains wrappers for CCBlade's Viterna extrapolation function. This file makes it easy to extrapolate 
            for multiple Reynolds numbers at once. Requires the user to supply existing airfoil data.
NOTES: Airfoil extrapolation takes existing 2D airfoil data and extrapolates it to angles of attack ranging from 
      -180 to 180 degrees. This is required for the Blade Element Momentum (BEM) solver in the wind turbine blade analysis.

(future contributors: feel free to add your name to "AUTHOR(S)")
=#

using CCBlade
using PyPlot
using DelimitedFiles
import FLOWMath


# ----- AIRFOIL EXTRAPOLATION FUNCTION WITHOUT ROTATIONAL CORRECTIONS -----

function extrapolate_no_rot(xfoildata, alpha_interpolate, Re; af_name="(name not specified)", af_directory="", ext_af_file="", cr75=0.128)
  # This function is used if rotational corrections will be applied on the fly during the BEM analysis.
  # (slightly more accurate than precomputing rotational corrections)

  # separate angle of attack, coefficient of lift, and coefficient of drag into individual vectors
  alpha = xfoildata[:, 1] * pi/180
  cl = xfoildata[:, 2]
  cd = xfoildata[:, 3]

  # extrapolate airfoil data to all angles of attack
  alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)  # Viterna extrapolation

  # interpolate onto specified alpha values
  cl_interpolate = FLOWMath.linear(alpha_ext, cl_ext, alpha_interpolate)
  cd_interpolate = FLOWMath.linear(alpha_ext, cd_ext, alpha_interpolate)

  # fix any wacky stuff in the airfoil data
  println(alpha_interpolate[abs(cl_interpolate) .> 10.0 * abs(cd_interpolate) .> 10.0])

  # save values to txt file
  open(af_directory * ext_af_file, "w") do io
    write(io, af_name * " airfoil data, extrapolated using CCBlade's Viterna extrapolation function\n$Re\n0\n")
    writedlm(io, [alpha_interpolate cl_interpolate cd_interpolate])
  end

end

# ----- AIRFOIL EXTRAPOLATION FUNCTION WITH ROTATIONAL CORRECTIONS -----

function extrapolate_rot(xfoildata, alpha_interpolate, Re; af_name="(name not specified)", af_directory="", ext_af_file="", cr75=0.128, tsr=5.0)
  # This function is used for factoring in precomputed rotational corrections into the extrapolated airfoil data.
  # (slightly less accurate than calculating rotational corrections on the fly)

  # separate angle of attack, coefficient of lift, and coefficient of drag into individual vectors
  alpha = xfoildata[:, 1] * pi/180
  cl = xfoildata[:, 2]
  cd = xfoildata[:, 3]

  # extrapolate airfoil data to all angles of attack
  alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)  # Viterna extrapolation

  # interpolate onto specified alpha values
  cl_interpolate = FLOWMath.linear(alpha_ext, cl_ext, alpha_interpolate)
  cd_interpolate = FLOWMath.linear(alpha_ext, cd_ext, alpha_interpolate)

  # define rotationally corrected alpha, cl, and cd vectors
  cl_rot = similar(cl_interpolate)
  cd_rot = similar(cd_interpolate)

  # choose representative blade location
  rR = 0.75  # r/R = 75%

  # apply rotational correction for lift and drag
  for i = 1:length(cl_ext)
    cl_rot[i], cd_rot[i] = rotation_correction(DuSeligEggers(), cl_interpolate[i], cd_interpolate[i], cr75, rR, tsr, alpha_interpolate[i])
  end

  # # save final data
  # airfoil = AlphaAF(alpha_interpolate, cl_rot, cd_rot, afdescription, Re, 0.0)
  # write_af(string(path, afname, ".dat"), airfoil)
  
  # save values to txt file
  open(af_directory * ext_af_file, "w") do io
    write(io, af_name * " airfoil data with rotational corrections, extrapolated using CCBlade's Viterna extrapolation function\n$Re\n0\n")
    writedlm(io, [alpha_interpolate cl_rot cd_rot])
  end

end


# =========================================================================================
# ======================= SET UP AIRFOIL DATA EXTRAPOLATION ===============================
# =========================================================================================

# =================== This block contains all the user inputs. ============================

# set airfoil data directory
af_directory = "airfoil-data/"

# specify name of airfoil
af_name = "NACA 4412"    # this should match exactly with the name of the directory that contains all the data for the airfoil of interest

# set the Reynolds numbers
Re_vec = [.005, .01, .06, .08, .1, .13, .15] .* 1e6    # you should already have nonextrapolated airfoil data for these Re numbers, saved in the appropriate location

# specify the form for the nonextrapolated airfoil data file name
get_afdata_nonext_filename(Re) = af_name * "_T1_Re" * rpad(Re/1e6, 5, '0') * "_M0.00_N9.0.dat"    # (we are creating a function here that takes in the Reynolds number as an argument)
                                                                                                  # QBlade airfoil data naming convention. Example: "NACA 4412_T1_Re0.100_M0.00_N9.0.dat" (for Re=0.1e6, Ma=0.0, and N_crit=9)

# specify how many header lines are in the nonextrapolated airfoil data file
n_headerlines = 14

# specify whether or not to apply rotational corrections to the extrapolated airfoil data
rotational_corrections = false    # set to false if rotational corrections will be applied on the fly (which is slightly more accurate)

# set tip-speed ratio (only applicable if applying rotational corrections to the airfoil data)
tsr = 5.0

# =========================================================================================
# =========================================================================================
# =========================================================================================


# ========= Everything below this point normally does not need to be edited. ==============


# get a vector of angles of attack
alpha_interpolate = readdlm("airfoil-data/alpha_360airfoildata.txt", skipstart=1)[:,1]    # all 360-degree airfoil data will be interpolated/extrapolated to these angles

# loop through each Reynolds number
for i = 1:length(Re_vec)
  # specify file name for the nonextrapolated airfoil data
  nonext_af_file = af_name * "/nonextrapolated/" * get_afdata_nonext_filename(Re_vec[i])
  # set file name for the extrapolated airfoil data
  ext_af_file = af_name * "/extrapolated/" * af_name * "_Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
  # read in nonextrapolated airfoil data
  nonext_af = readdlm(af_directory * nonext_af_file, skipstart=n_headerlines)
  # extrapolate airfoil data
  if rotational_corrections
    extrapolate_rot(nonext_af, alpha_interpolate, Re_vec[i]; af_name=af_name, af_directory=af_directory, ext_af_file=ext_af_file, cr75=0.128, tsr=tsr)
  else
    extrapolate_no_rot(nonext_af, alpha_interpolate, Re_vec[i]; af_name=af_name, af_directory=af_directory, ext_af_file=ext_af_file, cr75=0.128)
  end
end