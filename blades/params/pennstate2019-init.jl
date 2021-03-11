# ================================================
# ============ MACHINE GENERATED FILE ============
# ================================================


#= 
DESCRIPTION: Defines rotor and blade parameters/files for the NREL 5MW wind turbine blades.
NOTES: This file is intended to be used by a blade analysis and/or optimization script.
=#

# define rotor parameters
B = 3   # number of blades
turbine = true
precone = 0.0*pi/180

# define hub, tip, and radial sections
Rhub = .0225
Rtip = .225
nsections = 17  # number of radial sections
r_spacing = (Rtip - Rhub)/nsections # distance between each radial section
r = collect(range(0, stop=1.0, length=nsections)).*(Rtip - Rhub - r_spacing) .+ Rhub .+ r_spacing/2 # evenly space the radial sections along the length of the blade
# r = [0.025570089, 0.03889375, 0.052217411, 0.064974107, 0.07858125, 0.091904911, 0.104945089, 0.11826875, 0.131592411, 0.144916071, 0.158239732, 0.171563393, 0.184603571, 0.197927232, 0.211250893, 0.224574554]

# get the airfoil data
Re_vec = [0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3] .* 1e6
airfoildata_directory = "airfoil-data/wortmann-fx63137/"
wortmann_fx_63137_filenames = Vector{String}(undef, length(Re_vec))
for i = 1:length(Re_vec)
    wortmann_fx_63137_filenames[i] = airfoildata_directory * "wortmann-fx63137-extrapolated-Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
    alpha_vec = readdlm(wortmann_fx_63137_filenames[i], skipstart=3)[:,1]
end
wortmann_fx_63137 = AlphaReAF(wortmann_fx_63137_filenames; radians=true)

# Define airfoil(s). This blade only has one airfoil that is used throughout the entire length.
naf = length(r)
aftypes = Array{AlphaReAF}(undef, naf)
aftypes[1] = wortmann_fx_63137

# create airfoil array
af_idx = Int.(ones(length(r)))
airfoils = aftypes[af_idx]

# chord and twist
chord = [0.05357108463239126, 0.05808740494056436, 0.06117913386002035, 0.06222132005487386, 0.060643989883516325, 0.057379768605734384, 0.05350977841323318, 0.04974564648499446, 0.046799, 0.044473324835019666, 0.04217231251394997, 0.04000789868646254, 0.03809201900222904, 0.03645076215246175, 0.0349442251691442, 0.033593395003744635, 0.03242794020975264]
theta = [0.747000919853573, 0.747000919853573, 0.6665218925794738, 0.5425370264613579, 0.43341104961073607, 0.34522798749755074, 0.28407186559174363, 0.2483125625215775, 0.2187483488049216, 0.1940536007100472, 0.17422831823695445, 0.15927250138564317, 0.14817870126853225, 0.13729253601751715, 0.1272780233423265, 0.11927739554180508, 0.11443288491479785]

# chord and twist spline parameters
cspline = [0.050974, 0.061282, 0.046799, 0.037665, 0.031924]
tspline = [0.747000919853573, 0.2719048441681966, 0.15465362501921753, 0.11355112113475108]

# optimal tip-speed ratio
tsr = 4.0

# optimal pitch angles
pitch = [0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659]