using DelimitedFiles
using PyPlot

# set airfoil data file directory and name
data_file_directory = "airfoil-data/wortmann-fx63137/"
data_file_name = "wortmann-fx63137-extrapolated-Re0.005e6.dat"

# import airfoil data
airfoil_data = readdlm(data_file_directory * data_file_name, skipstart=3)
alpha = airfoil_data[:,1] .* 180/pi
c_l = airfoil_data[:,2]
c_d = airfoil_data[:,3]

