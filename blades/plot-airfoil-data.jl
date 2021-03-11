using DelimitedFiles
using PyPlot
# matplotlib.use("tkagg")

# set airfoil data directory
af_directory = "airfoil-data/"

# specify name of airfoil
af_name = "NACA 4412"    # this should match exactly with the name of the directory that contains all the data for the airfoil of interest

# set the Reynolds numbers
Re_vec = 0#[.005, .01, .06, .08, .1, .13, .15] .* 1e6    # you should already have airfoil data for these Re numbers, saved in the appropriate location

# loop through each Reynolds number
for i = 1:length(Re_vec)
    
    # get airfoil file name
    # af_file = af_name * "/extrapolated/" * af_name * "_Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
    af_file = "example-airfoil-data/DU25_A17.dat"

    # import airfoil data
    af_data = readdlm(af_directory * af_file, skipstart=3)
    alpha = af_data[:,1]
    if maximum(af_data[:,1]) < 10.0
        # radians
        alpha *= 180/pi
    end
    c_l = af_data[:,2]
    c_d = af_data[:,3]

    figure()
    plot(alpha, c_l)
    plot(alpha, c_d)
    xlabel(L"\alpha")
    ylim([-1., 1.75])
    legend([L"c_L", L"c_D"])
    plt.show()

    figure()
    plot(alpha, c_l)
    plot(alpha, c_d)
    xlabel(L"\alpha")
    legend([L"c_L", L"c_D"])
    xlim([-15.0, 30.0])
    ylim([-1., 1.75])
    plt.show()

    figure()
    plot(alpha, c_l)
    plot(alpha, c_d)
    xlabel(L"\alpha")
    legend([L"c_L", L"c_D"])
    xlim([-15.0, 40.0])
    ylim([-0.05, 0.4])
    plt.show()
end