# Size of variable arrays:
sizeAlgebraic = 7
sizeStates = 3
sizeConstants = 20
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_VOI = "time in component main (ms)"
    legend_states[0] = "V in component main (mV)"
    legend_states[1] = "h in component main (dimensionless)"
    legend_states[2] = "n in component main (dimensionless)"
    legend_algebraic[6] = "i_Na in component main (dimensionless)"
    legend_algebraic[5] = "i_K in component main (dimensionless)"
    legend_constants[0] = "G_Na in component main (dimensionless)"
    legend_constants[1] = "E_Na in component main (dimensionless)"
    legend_constants[2] = "East in component main (dimensionless)"
    legend_constants[3] = "Edag in component main (dimensionless)"
    legend_constants[4] = "g_2 in component main (dimensionless)"
    legend_constants[5] = "F_n in component main (dimensionless)"
    legend_constants[6] = "k_1 in component main (dimensionless)"
    legend_constants[7] = "k_2 in component main (dimensionless)"
    legend_constants[8] = "k_3 in component main (dimensionless)"
    legend_constants[9] = "E_1 in component main (dimensionless)"
    legend_constants[10] = "E_2 in component main (dimensionless)"
    legend_constants[11] = "E_3 in component main (dimensionless)"
    legend_constants[12] = "F_h in component main (dimensionless)"
    legend_constants[13] = "epsilon in component main (dimensionless)"
    legend_constants[14] = "epsilon_2 in component main (dimensionless)"
    legend_algebraic[3] = "G in component main (dimensionless)"
    legend_algebraic[0] = "H_1 in component main (dimensionless)"
    legend_algebraic[1] = "H_2 in component main (dimensionless)"
    legend_algebraic[2] = "H_3 in component main (dimensionless)"
    legend_algebraic[4] = "Istim in component main (dimensionless)"
    legend_constants[15] = "IstimStart in component main (dimensionless)"
    legend_constants[16] = "IstimEnd in component main (dimensionless)"
    legend_constants[17] = "IstimAmplitude in component main (dimensionless)"
    legend_constants[18] = "IstimPeriod in component main (dimensionless)"
    legend_constants[19] = "IstimPulseDuration in component main (dimensionless)"
    legend_rates[2] = "d/dt n in component main (dimensionless)"
    legend_rates[1] = "d/dt h in component main (dimensionless)"
    legend_rates[0] = "d/dt V in component main (mV)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -93.3333
    states[1] = 1.0
    states[2] = 0.0
    constants[0] = 33.33333
    constants[1] = 40.0
    constants[2] = -15.0
    constants[3] = -80.0
    constants[4] = -9.0
    constants[5] = 0.0037037
    constants[6] = 0.075
    constants[7] = 0.04
    constants[8] = 0.1
    constants[9] = -93.3333
    constants[10] = -55.0
    constants[11] = 1.0
    constants[12] = 0.5
    constants[13] = 1.0
    constants[14] = 1.0
    constants[15] = 0
    constants[16] = 50000
    constants[17] = 80.0
    constants[18] = 1000
    constants[19] = 1.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = custom_piecewise([less_equal(states[0] , constants[3]), 1.00000 , True, 0.000000])
    rates[1] = (constants[12]*(algebraic[1]-states[1]))/constants[13]
    algebraic[2] = custom_piecewise([greater_equal(states[0] , constants[3]), 1.00000 , True, 0.000000])
    rates[2] = constants[14]*constants[5]*(algebraic[2]-states[2])
    algebraic[0] = custom_piecewise([greater_equal(states[0] , constants[2]), 1.00000 , True, 0.000000])
    algebraic[6] = (constants[0]*states[1]*(constants[1]-states[0])*algebraic[0])/constants[13]
    algebraic[3] = custom_piecewise([less(states[0] , constants[3]), constants[6]*(constants[9]-states[0]) , greater_equal(states[0] , constants[2]), constants[8]*(constants[11]-states[0]) , True, constants[7]*(states[0]-constants[10])])
    algebraic[5] = constants[4]*algebraic[2]*(states[2]**4.00000)+algebraic[3]
    algebraic[4] = custom_piecewise([greater_equal(VOI , constants[15]) & less_equal(VOI , constants[16]) & less_equal((VOI-constants[15])-(floor(((VOI-constants[15])/constants[18])))*constants[18] , constants[19]), constants[17] , True, 0.000000])
    rates[0] = algebraic[4]+algebraic[6]+algebraic[5]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = custom_piecewise([less_equal(states[0] , constants[3]), 1.00000 , True, 0.000000])
    algebraic[2] = custom_piecewise([greater_equal(states[0] , constants[3]), 1.00000 , True, 0.000000])
    algebraic[0] = custom_piecewise([greater_equal(states[0] , constants[2]), 1.00000 , True, 0.000000])
    algebraic[6] = (constants[0]*states[1]*(constants[1]-states[0])*algebraic[0])/constants[13]
    algebraic[3] = custom_piecewise([less(states[0] , constants[3]), constants[6]*(constants[9]-states[0]) , greater_equal(states[0] , constants[2]), constants[8]*(constants[11]-states[0]) , True, constants[7]*(states[0]-constants[10])])
    algebraic[5] = constants[4]*algebraic[2]*(states[2]**4.00000)+algebraic[3]
    algebraic[4] = custom_piecewise([greater_equal(VOI , constants[15]) & less_equal(VOI , constants[16]) & less_equal((VOI-constants[15])-(floor(((VOI-constants[15])/constants[18])))*constants[18] , constants[19]), constants[17] , True, 0.000000])
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 10000, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-006, rtol=1e-006, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)

