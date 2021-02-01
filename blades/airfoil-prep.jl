using CCBlade
using PyPlot
using DelimitedFiles
import FLOWMath

function extrapolate_no_rot(xfoildata, alpha_interpolate, Re; afname="airfoil", afdescription="", path="", cr75=0.128)

  # separate angle of attack, coefficient of lift, and coefficient of drag into individual vectors
  alpha = xfoildata[:, 1] * pi/180
  cl = xfoildata[:, 2]
  cd = xfoildata[:, 3]

  # extrapolate airfoil data to all angles of attack
  alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)  # Viterna extrapolation

  # interpolate onto specified alpha values
  cl_interpolate = FLOWMath.linear(alpha_ext, cl_ext, alpha_interpolate)
  cd_interpolate = FLOWMath.linear(alpha_ext, cd_ext, alpha_interpolate)

  # save values to txt file
  open(path * afname, "w") do io
    write(io, "Wortmann FX 63-137 airfoil data\n$Re\n0\n")
    writedlm(io, [alpha_interpolate cl_interpolate cd_interpolate])
  end

end

function extrapolate_rot(xfoildata, alpha_interpolate, Re; afname="airfoil", afdescription="", path="", tsr=6.0, cr75=0.128)

  # separate angle of attack, coefficient of lift, and coefficient of drag into individual vectors
  alpha = xfoildata[:, 1] * pi/180
  cl = xfoildata[:, 2]
  cd = xfoildata[:, 3]

  # extrapolate airfoil data to all angles of attack
  alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)  # Viterna extrapolation

  # interpolate onto specified alpha values
  cl_ext = FLOWMath.linear(alpha_ext, cl_ext, alpha_interpolate)
  cd_ext = FLOWMath.linear(alpha_ext, cd_ext, alpha_interpolate)

  # define rotationally corrected alpha, cl, and cd vectors
  alpha_rot = alpha_interpolate
  cl_rot = similar(cl_ext)
  cd_rot = similar(cd_ext)

  # choose representative blade location
  rR = 0.75  # r/R = 75%

  # apply rotational correction for lift and drag
  for i = 1:length(cl_ext)
    cl_rot[i], cd_rot[i] = rotation_correction(DuSeligEggers(), cl_ext[i], cd_ext[i], cr75, rR, tsr, alpha_rot[i])
  end

  # save final data
  airfoil = AlphaAF(alpha_rot, cl_rot, cd_rot, afdescription, Re, 0.0)
  write_af(string(path, afname, ".dat"), airfoil)
  
end


xfoildata_directory = "airfoil-data/wortmann-fx63137/"
Re_vec = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 
    0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3] .* 1e6
alpha_interpolate = readdlm("airfoil-data/alpha_360airfoildata.txt", skipstart=1)[:,1]

for i = 1:length(Re_vec)
  xfoildata_filename = "WORTMANN FX 63-137 AIRFOIL EXTRA REFINED 2_T1_Re" * rpad(Re_vec[i]/1e6, 5, '0') * "_M0.00_N9.0.dat"
  extrapolateddata_filename = "wortmann-fx63137-extrapolated-Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
  afdescription = "Wortmann FX 63-137 airfoil data"
  xfoildata = readdlm(xfoildata_directory * xfoildata_filename, skipstart=14)
  extrapolate_no_rot(xfoildata, alpha_interpolate, Re_vec[i]; afname=extrapolateddata_filename, afdescription=afdescription, path=xfoildata_directory, cr75=0.128)
end