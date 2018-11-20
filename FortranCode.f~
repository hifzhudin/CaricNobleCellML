C
C There are a total of 7 entries in the algebraic variable array.
C There are a total of 3 entries in each of the rate and state variable arrays.
C There are a total of 20 entries in the constant variable array.
C
C
C VOI is time in component main (ms).
C STATES(0) is V in component main (mV).
C STATES(1) is h in component main (dimensionless).
C STATES(2) is n in component main (dimensionless).
C ALGBRC(6) is i_Na in component main (dimensionless).
C ALGBRC(5) is i_K in component main (dimensionless).
C CONSTS(0) is G_Na in component main (dimensionless).
C CONSTS(1) is E_Na in component main (dimensionless).
C CONSTS(2) is East in component main (dimensionless).
C CONSTS(3) is Edag in component main (dimensionless).
C CONSTS(4) is g_2 in component main (dimensionless).
C CONSTS(5) is F_n in component main (dimensionless).
C CONSTS(6) is k_1 in component main (dimensionless).
C CONSTS(7) is k_2 in component main (dimensionless).
C CONSTS(8) is k_3 in component main (dimensionless).
C CONSTS(9) is E_1 in component main (dimensionless).
C CONSTS(10) is E_2 in component main (dimensionless).
C CONSTS(11) is E_3 in component main (dimensionless).
C CONSTS(12) is F_h in component main (dimensionless).
C CONSTS(13) is epsilon in component main (dimensionless).
C CONSTS(14) is epsilon_2 in component main (dimensionless).
C ALGBRC(3) is G in component main (dimensionless).
C ALGBRC(0) is H_1 in component main (dimensionless).
C ALGBRC(1) is H_2 in component main (dimensionless).
C ALGBRC(2) is H_3 in component main (dimensionless).
C ALGBRC(4) is Istim in component main (dimensionless).
C CONSTS(15) is IstimStart in component main (dimensionless).
C CONSTS(16) is IstimEnd in component main (dimensionless).
C CONSTS(17) is IstimAmplitude in component main (dimensionless).
C CONSTS(18) is IstimPeriod in component main (dimensionless).
C CONSTS(19) is IstimPulseDuration in component main (dimensionless).
C RATES(2) is d/dt n in component main (dimensionless).
C RATES(1) is d/dt h in component main (dimensionless).
C RATES(0) is d/dt V in component main (mV).
C

SUBROUTINE initConsts(CONSTS, RATES, STATES)
      REAL CONSTS(*), RATES(*), STATES(*)
      STATES(0) = -93.3333
      STATES(1) = 1.0
      STATES(2) = 0.0
      CONSTS(0) = 33.33333
      CONSTS(1) = 40.0
      CONSTS(2) = -15.0
      CONSTS(3) = -80.0
      CONSTS(4) = -9.0
      CONSTS(5) = 0.0037037
      CONSTS(6) = 0.075
      CONSTS(7) = 0.04
      CONSTS(8) = 0.1
      CONSTS(9) = -93.3333
      CONSTS(10) = -55.0
      CONSTS(11) = 1.0
      CONSTS(12) = 0.5
      CONSTS(13) = 1.0
      CONSTS(14) = 1.0
      CONSTS(15) = 0
      CONSTS(16) = 50000
      CONSTS(17) = 80.0
      CONSTS(18) = 1000
      CONSTS(19) = 1.0
      RETURN
      END
      SUBROUTINE computeRates(VOI, CONSTS,  RATES, STATES, ALGBRC)
      REAL VOI, CONSTS(*), RATES(*), STATES(*), ALGBRC(*)
      ALGBRC(1) = TERNRY(STATES(0).LE.CONSTS(3), 1.00000, 0.000000)
      RATES(1) = ( CONSTS(12)*(ALGBRC(1) - STATES(1)))/CONSTS(13)
      ALGBRC(2) = TERNRY(STATES(0).GE.CONSTS(3), 1.00000, 0.000000)
      RATES(2) =  CONSTS(14)*CONSTS(5)*(ALGBRC(2) - STATES(2))
      ALGBRC(0) = TERNRY(STATES(0).GE.CONSTS(2), 1.00000, 0.000000)
      ALGBRC(6) = ( CONSTS(0)*STATES(1)*(CONSTS(1) - STATES(0))*ALGBRC(0))/CONSTS(13)
      ALGBRC(3) = TERNRY(STATES(0).LT.CONSTS(3),  CONSTS(6)*(CONSTS(9) - STATES(0)), TERNRY(STATES(0).GE.CONSTS(2),  CONSTS(8)*(CONSTS(11) - STATES(0)),  CONSTS(7)*(STATES(0) - CONSTS(10)))
      ALGBRC(5) =  CONSTS(4)*ALGBRC(2)*(STATES(2) ** 4.00000)+ALGBRC(3)
      ALGBRC(4) = TERNRY(VOI.GE.CONSTS(15).AND.VOI.LE.CONSTS(16).AND.(VOI - CONSTS(15)) -  (INT(((VOI - CONSTS(15))/CONSTS(18))))*CONSTS(18).LE.CONSTS(19), CONSTS(17), 0.000000)
      RATES(0) = ALGBRC(4)+ALGBRC(6)+ALGBRC(5)
      RETURN
      END
      SUBROUTINE computeVariables(VOI, CONSTS, RATES, STATES, ALGBRC)
      REAL VOI, CONSTS(*), RATES(*), STATES(*), ALGBRC(*)
      ALGBRC(1) = TERNRY(STATES(0).LE.CONSTS(3), 1.00000, 0.000000)
      ALGBRC(2) = TERNRY(STATES(0).GE.CONSTS(3), 1.00000, 0.000000)
      ALGBRC(0) = TERNRY(STATES(0).GE.CONSTS(2), 1.00000, 0.000000)
      ALGBRC(6) = ( CONSTS(0)*STATES(1)*(CONSTS(1) - STATES(0))*ALGBRC(0))/CONSTS(13)
      ALGBRC(3) = TERNRY(STATES(0).LT.CONSTS(3),  CONSTS(6)*(CONSTS(9) - STATES(0)), TERNRY(STATES(0).GE.CONSTS(2),  CONSTS(8)*(CONSTS(11) - STATES(0)),  CONSTS(7)*(STATES(0) - CONSTS(10)))
      ALGBRC(5) =  CONSTS(4)*ALGBRC(2)*(STATES(2) ** 4.00000)+ALGBRC(3)
      ALGBRC(4) = TERNRY(VOI.GE.CONSTS(15).AND.VOI.LE.CONSTS(16).AND.(VOI - CONSTS(15)) -  (INT(((VOI - CONSTS(15))/CONSTS(18))))*CONSTS(18).LE.CONSTS(19), CONSTS(17), 0.000000)
      RETURN
      END
      REAL FUNCTION TERNRY(TEST, VALA, VALB)
      LOGICAL TEST
      REAL VALA, VALB
      IF (TEST) THEN
        TERNRY = VALA
      ELSE
        TERNRY = VALB
      ENDIF
      RETURN
      END

