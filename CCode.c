/*
   There are a total of 7 entries in the algebraic variable array.
   There are a total of 3 entries in each of the rate and state variable arrays.
   There are a total of 20 entries in the constant variable array.
 */
/*
 * VOI is time in component main (ms).
 * STATES[0] is V in component main (mV).
 * STATES[1] is h in component main (dimensionless).
 * STATES[2] is n in component main (dimensionless).
 * ALGEBRAIC[6] is i_Na in component main (dimensionless).
 * ALGEBRAIC[5] is i_K in component main (dimensionless).
 * CONSTANTS[0] is G_Na in component main (dimensionless).
 * CONSTANTS[1] is E_Na in component main (dimensionless).
 * CONSTANTS[2] is East in component main (dimensionless).
 * CONSTANTS[3] is Edag in component main (dimensionless).
 * CONSTANTS[4] is g_2 in component main (dimensionless).
 * CONSTANTS[5] is F_n in component main (dimensionless).
 * CONSTANTS[6] is k_1 in component main (dimensionless).
 * CONSTANTS[7] is k_2 in component main (dimensionless).
 * CONSTANTS[8] is k_3 in component main (dimensionless).
 * CONSTANTS[9] is E_1 in component main (dimensionless).
 * CONSTANTS[10] is E_2 in component main (dimensionless).
 * CONSTANTS[11] is E_3 in component main (dimensionless).
 * CONSTANTS[12] is F_h in component main (dimensionless).
 * CONSTANTS[13] is epsilon in component main (dimensionless).
 * CONSTANTS[14] is epsilon_2 in component main (dimensionless).
 * ALGEBRAIC[3] is G in component main (dimensionless).
 * ALGEBRAIC[0] is H_1 in component main (dimensionless).
 * ALGEBRAIC[1] is H_2 in component main (dimensionless).
 * ALGEBRAIC[2] is H_3 in component main (dimensionless).
 * ALGEBRAIC[4] is Istim in component main (dimensionless).
 * CONSTANTS[15] is IstimStart in component main (dimensionless).
 * CONSTANTS[16] is IstimEnd in component main (dimensionless).
 * CONSTANTS[17] is IstimAmplitude in component main (dimensionless).
 * CONSTANTS[18] is IstimPeriod in component main (dimensionless).
 * CONSTANTS[19] is IstimPulseDuration in component main (dimensionless).
 * RATES[2] is d/dt n in component main (dimensionless).
 * RATES[1] is d/dt h in component main (dimensionless).
 * RATES[0] is d/dt V in component main (mV).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = -93.3333;
STATES[1] = 1.0;
STATES[2] = 0.0;
CONSTANTS[0] = 33.33333;
CONSTANTS[1] = 40.0;
CONSTANTS[2] = -15.0;
CONSTANTS[3] = -80.0;
CONSTANTS[4] = -9.0;
CONSTANTS[5] = 0.0037037;
CONSTANTS[6] = 0.075;
CONSTANTS[7] = 0.04;
CONSTANTS[8] = 0.1;
CONSTANTS[9] = -93.3333;
CONSTANTS[10] = -55.0;
CONSTANTS[11] = 1.0;
CONSTANTS[12] = 0.5;
CONSTANTS[13] = 1.0;
CONSTANTS[14] = 1.0;
CONSTANTS[15] = 0;
CONSTANTS[16] = 50000;
CONSTANTS[17] = 80.0;
CONSTANTS[18] = 1000;
CONSTANTS[19] = 1.0;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[1] = (STATES[0]<=CONSTANTS[3] ? 1.00000 : 0.000000);
RATES[1] = ( CONSTANTS[12]*(ALGEBRAIC[1] - STATES[1]))/CONSTANTS[13];
ALGEBRAIC[2] = (STATES[0]>=CONSTANTS[3] ? 1.00000 : 0.000000);
RATES[2] =  CONSTANTS[14]*CONSTANTS[5]*(ALGEBRAIC[2] - STATES[2]);
ALGEBRAIC[0] = (STATES[0]>=CONSTANTS[2] ? 1.00000 : 0.000000);
ALGEBRAIC[6] = ( CONSTANTS[0]*STATES[1]*(CONSTANTS[1] - STATES[0])*ALGEBRAIC[0])/CONSTANTS[13];
ALGEBRAIC[3] = (STATES[0]<CONSTANTS[3] ?  CONSTANTS[6]*(CONSTANTS[9] - STATES[0]) : STATES[0]>=CONSTANTS[2] ?  CONSTANTS[8]*(CONSTANTS[11] - STATES[0]) :  CONSTANTS[7]*(STATES[0] - CONSTANTS[10]));
ALGEBRAIC[5] =  CONSTANTS[4]*ALGEBRAIC[2]*(pow(STATES[2], 4.00000))+ALGEBRAIC[3];
ALGEBRAIC[4] = (VOI>=CONSTANTS[15]&&VOI<=CONSTANTS[16]&&(VOI - CONSTANTS[15]) -  (floor(((VOI - CONSTANTS[15])/CONSTANTS[18])))*CONSTANTS[18]<=CONSTANTS[19] ? CONSTANTS[17] : 0.000000);
RATES[0] = ALGEBRAIC[4]+ALGEBRAIC[6]+ALGEBRAIC[5];
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[1] = (STATES[0]<=CONSTANTS[3] ? 1.00000 : 0.000000);
ALGEBRAIC[2] = (STATES[0]>=CONSTANTS[3] ? 1.00000 : 0.000000);
ALGEBRAIC[0] = (STATES[0]>=CONSTANTS[2] ? 1.00000 : 0.000000);
ALGEBRAIC[6] = ( CONSTANTS[0]*STATES[1]*(CONSTANTS[1] - STATES[0])*ALGEBRAIC[0])/CONSTANTS[13];
ALGEBRAIC[3] = (STATES[0]<CONSTANTS[3] ?  CONSTANTS[6]*(CONSTANTS[9] - STATES[0]) : STATES[0]>=CONSTANTS[2] ?  CONSTANTS[8]*(CONSTANTS[11] - STATES[0]) :  CONSTANTS[7]*(STATES[0] - CONSTANTS[10]));
ALGEBRAIC[5] =  CONSTANTS[4]*ALGEBRAIC[2]*(pow(STATES[2], 4.00000))+ALGEBRAIC[3];
ALGEBRAIC[4] = (VOI>=CONSTANTS[15]&&VOI<=CONSTANTS[16]&&(VOI - CONSTANTS[15]) -  (floor(((VOI - CONSTANTS[15])/CONSTANTS[18])))*CONSTANTS[18]<=CONSTANTS[19] ? CONSTANTS[17] : 0.000000);
}

