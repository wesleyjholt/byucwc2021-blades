using CCBlade
using PyPlot
import FLOWMath

function extrapolate_no_rot(xfoildata, alpha_interpolate, Re, afname="airfoil", afdescription="", path="", cr75=0.128)

  # separate angle of attack, coefficient of lift, and coefficient of drag into individual vectors
  alpha = xfoildata[:, 1] * pi/180
  cl = xfoildata[:, 2]
  cd = xfoildata[:, 3]

  # extrapolate airfoil data to all angles of attack
  alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)  # Viterna extrapolation

  # interpolate onto specified alpha values
  cl_ext = FLOWMath.linear(alpha_ext, cl_ext, alpha_interpolate)
  cd_ext = FLOWMath.linear(alpha_ext, cd_ext, alpha_interpolate)

  # save final data
  airfoil = AlphaAF(alpha_ext, cl_ext, cd_ext, afdescription, Re, 0.0)
  write_af(string(path, afname, ".dat"), airfoil)
  
end

function extrapolate_rot(xfoildata, alpha_interpolate, Re, afname="airfoil", afdescription="", path="", tsr=6.0, cr75=0.128)

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

