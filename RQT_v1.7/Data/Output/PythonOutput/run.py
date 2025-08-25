import sys
sys.path.append('/Users/iwasakihiroto/pyrate/Output/PythonOutput')

from SM import RGEsolver

##############################################
# First, create an instance of the RGEsolver #
##############################################

rge = RGEsolver('rge', tmin=0, tmax=20, initialScale=0)


##########################################################
# We fix the running scheme and initial conditions below #
##########################################################

# Running scheme :

rge.loops = {'GaugeCouplings': 3,
             'Yukawas': 2,
             'QuarticTerms': 2,
             'ScalarMasses': 2,
             'Vevs': 2}

# Gauge Couplings

rge.g1.initialValue = 0
rge.g2.initialValue = 0
rge.g3.initialValue = 0

# Yukawa Couplings

rge.Yu.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]

rge.Yd.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]

rge.Ye.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]


# Quartic Couplings

rge.lambda_.initialValue = 0

# Scalar Mass Couplings

rge.mu.initialValue = 0

# Vacuum-expectation Values

rge.v.initialValue = 0


############################
# Solve the system of RGEs #
############################

rge.solve(step = .05)

# Another way to call rge.solve() :
# rge.solve(Npoints = 500)

####################
# Plot the results #
####################

rge.plot(subPlots=True, printLoopLevel=True)


#############################################
# Possibly save the results for a later use #
#############################################

# Save the results in some file

# rge.save('rgeResults.save')

# Later, load the rge object with :

# rge = RGEsolver.load('rgeResults.save')

