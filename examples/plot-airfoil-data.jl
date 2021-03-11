using DelimitedFiles
using PyPlot

# set airfoil data file directory and name
data_file_directory = "data/"
data_file_name = "DU21_A17.dat"

# import airfoil data
airfoil_data = readdlm(data_file_directory * data_file_name, skipstart=3)
alpha = airfoil_data[:,1]
cl = airfoil_data[:,2]
cd = airfoil_data[:,3]

figure()
plot(alpha, cl)
plot(alpha, cd)
xlabel(L"\alpha")
legend([L"c_L", L"c_D"])
plt.show()

figure()
plot(alpha, cl)
plot(alpha, cd)
xlabel(L"\alpha")
legend([L"c_L", L"c_D"])
xlim([-15.0, 40.0])
plt.show()

figure()
plot(alpha, cl)
plot(alpha, cd)
xlabel(L"\alpha")
legend([L"c_L", L"c_D"])
xlim([-15.0, 40.0])
ylim([-0.05, 0.3])
plt.show()