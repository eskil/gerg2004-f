      PROGRAM GERG2004
C===============================================================================
C     "The partial GERG-2004 implementation used to calculate the
C     compressibility factors of oxygen, helium, nitrogen and their
C     mixtures in the gas phase was inspired by open source code written
C     by David de Marneffe in December 2007. David can be contacted at
C     daviddemarneffe@yahoo.com."
C
C     Ported to Fortran by Eskil Heyn Olsen (eskil at eskil dot org).
C
C===============================================================================
      IMPLICIT NONE
      REAL Fraction_Oxygen, Fraction_Helium
      REAL Pressure, Temperature, Delta
      CHARACTER System*10
      REAL Z_VALUE, Z

C===============================================================================
C     Load input
      OPEN (UNIT=7, FILE='input', STATUS='UNKNOWN', ACCESS='SEQUENTIAL',
     *     FORM='FORMATTED')
      READ (7,*) Fraction_Oxygen, Fraction_Helium,
     *     Pressure, Temperature, Delta, System

C     Frob input
      IF (System .EQ. 'imperial') THEN
         Pressure = Pressure / 14.5037738
         Temperature = (Temperature - 32) / 1.8
      END IF

C     Convert pressure to absolute
      Pressure = Pressure + 1.01325
C     Convert temperature to kelvin
      Temperature = Temperature + 273.15

      Z = Z_VALUE(Fraction_Oxygen, Fraction_Helium, Pressure, 
     *     Temperature, Delta)

      WRITE (*,500) Fraction_Oxygen, Fraction_Helium
      WRITE (*,501) Temperature
      WRITE (*,502) Pressure
      WRITE (*,503) Z

C===============================================================================
 500  FORMAT (F5.2, '/', F5.2)
 501  FORMAT (F6.1, 'k')
 502  FORMAT (F6.1, 'bar(a)')
 503  FORMAT ('Z = ', F20.16)

C===============================================================================
      CLOSE (UNIT = 7, STATUS = 'KEEP')
      END

C===============================================================================
      FUNCTION Z_VALUE(fO2, fHe, Pressure, Temperature, MaxDelta)
      REAL fO2, fHe
      REAL Pressure, Temperature, MaxDelta
      REAL Z_VALUE

      REAL R, Tau, Dens, Dens_Red, Delta, Z_Old, ZDelta
      
      R = 0.08314472
      Tau = TEMP_R(fO2, fHe) / Temperature
      Dens = 1.0
      Dens_Red = DENS_R(fO2, fHe);

C     Initial guess
      Z_VALUE = 1.0

      Delta = 0.0
      WRITE (*, 601) 'Pressure', Pressure
      WRITE (*, 601) 'Z', Z_VALUE
      WRITE (*, 601) 'R', R
      WRITE (*, 601) 'Temp', Temperature
      WRITE (*, 601) 'Dens_R', Dens_Red
      WRITE (*, 601) 'Tau', Tau
      WRITE (*, 602)

      DO WHILE (.TRUE.)
C        calculate molar density (mole per liter)
         Dens = Pressure / (Z_VALUE * R * Temperature)
         WRITE (*, 601) 'Dens', Dens

C        calculate reduced density
         Delta = Dens / Dens_Red
         WRITE (*, 601) 'Delta', Delta

C        remember previous value
         Z_Old = Z_VALUE

C        calculate new Z
         Z_VALUE = 1.0 + Delta * 
     *        DALPHAR_DDELTA_MIX(fO2, fHe, Delta, Tau)
         WRITE (*, 601) 'Z', Z_VALUE

         ZDelta = ABS (Z_VALUE - Z_Old)
         WRITE (*, 601) 'ZDelta', ZDelta
         WRITE (*, 601) 'MaxDelta', MaxDelta
         IF (ZDelta .LT. MaxDelta) THEN
            EXIT
         END IF
         WRITE (*, 602)
      END DO

 601  FORMAT (A20, F16.12)
 602  FORMAT ()

      RETURN
      END


C===============================================================================
C     Checked
      FUNCTION TEMP_R(fO2, fHe)
      REAL fO2, fHe, fN2
      REAL Temp_cN2, Temp_cO2, Temp_cHe
      REAL BTN2O2, BTN2He, BTO2He
      REAL GTN2O2, GTN2He, GTO2He
      REAL TEMP_R

      fN2 = 1.0 - fO2 - fHe
      Temp_cN2 = 126.192
      Temp_cO2 = 154.595
      Temp_cHe = 5.1953
      BTN2O2 = 0.997190589
      BTN2He = 0.692868765
      BTO2He = 1.0
      GTN2O2 = 0.995157044
      GTN2He = 1.47183158
      GTO2He = 1.0
      
      IF ((fO2 .NE. 1.0).AND.(fHe .NE. 1.0).AND.(fN2 .NE. 1.0)) THEN
         TEMP_R = fN2**2.0 * temp_cN2 + fO2**2.0 * temp_cO2 + 
     *        fHe**2.0 * temp_cHe
     *        + 2.0 * fN2 * fO2 * BTN2O2 * GTN2O2 * 
     *        (fN2 + fO2) / (BTN2O2**2.0 * fN2 + fO2)
     *        * (temp_cN2 * temp_cO2)**0.5
     *        + 2.0 * fN2 * fHe * BTN2He * GTN2He * 
     *        (fN2 + fHe) / (BTN2He**2.0 * fN2 + fHe)
     *        * (temp_cN2 * temp_cHe)**0.5
     *        + 2.0 * fO2 * fHe * BTO2He * GTO2He * 
     *        (fO2 + fHe) / (BTO2He**2.0 * fO2 + fHe)
     *        * (temp_cO2 * temp_cHe)**0.5

      ELSE
         TEMP_R = fN2**2 * Temp_cN2 + 
     *        fO2**2 * Temp_cO2 + fHe**2 * Temp_cHe
      END IF
      RETURN
      END


C===============================================================================
C     Checked
      FUNCTION DENS_R(fO2, fHe)
      REAL fO2, fHe, fN2
      REAL DENS_R
      REAL dens_cN2, dens_cO2, dens_cHe
      REAL BmuN2O2, BmuN2He, BmuO2He
      REAL GmuN2O2, GmuN2He, GmuO2He
      REAL InvDensR
      
      fN2 = 1.0 - fO2 - fHe
      dens_cN2 = 11.1839;
      dens_cO2 = 13.63;
      dens_cHe = 17.399;
      BmuN2O2 = 0.99952177
      BmuN2He = 0.969501055
      BmuO2He = 1.0
      GmuN2O2 = 0.997082328
      GmuN2He = 0.932629867
      GmuO2He = 1.0
        
      IF ((fO2 .NE. 1.0).AND.(fHe .NE. 1.0).AND.(fN2 .NE. 1.0)) THEN
         InvDensR = 
     *        fN2**2 / dens_cN2 + fO2**2 / dens_cO2 + fHe**2 /
     *        dens_cHe 
     *        + 0.25 * fN2 * fO2 * BmuN2O2 * GmuN2O2 * (fN2 + fO2) / 
     *        (BmuN2O2**2 * fN2 + fO2) 
     *        * (1.0 / dens_cN2**(1.0 / 3.0) + 
     *        1.0 / dens_cO2**(1.0 / 3.0))**3 
     *        + 0.25 * fN2 * fHe * BmuN2He * GmuN2He * (fN2 + fHe) / 
     *        (BmuN2He**2 * fN2 + fHe) 
     *        * (1.0 / dens_cN2**(1.0 / 3.0) + 
     *        1.0 / dens_cHe**(1.0 / 3.0))**3 
     *        + 0.25 * fO2 * fHe * BmuO2He * GmuO2He * (fO2 + fHe) / 
     *        (BmuO2He**2 * fO2 + fHe) 
     *        * (1.0 / dens_cO2**(1.0 / 3.0) + 
     *        1.0 / dens_cHe**(1.0 / 3.0))**3
      ELSE
         InvDensR = fN2**2 / dens_cN2 + fO2**2 / dens_cO2 + 
     *        fHe**2 / dens_cHe
      END IF 

      DENS_R = 1.0 / InvDensR

      RETURN
      END


C===============================================================================
      FUNCTION DALPHAR_DDELTA_MIX(fO2, fHe, Delta, Tau)
      REAL fO2, fHe, fN2
      REAL Delta, Tau
      REAL DALPHAR_DDELTA_MIX

      fN2 = 1.0 - fO2 - fHe

      DALPHAR_DDELTA_MIX = fO2 * DALPHAR_DDELTA_O2 (Delta, Tau) +
     *     fHe * DALPHAR_DDELTA_He (Delta, Tau) +
     *     fN2 * DALPHAR_DDELTA_N2 (Delta, Tau)


      RETURN
      END


C===============================================================================
      FUNCTION DALPHAR_DDELTA_O2(Delta, Tau)
      REAL Delta, Tau
      REAL DALPHAR_DDELTA_O2

      DALPHAR_DDELTA_O2 = 0.88878286369701 * tau**0.25
     *     - 2.4879433312148 * tau**1.125
     *     + 0.59750190775886 * tau**1.5
     *     + 9.6501817061881E-03 * 2.0 * delta * tau**1.375
     *     + 0.07197042871277 * 3.0 * delta**2.0 * tau**0.25
     *     + 2.2337443000195E-04 * 7.0 * delta**6.0 * tau**0.875
     *     + 0.18558686391474 * delta * (2.0 - delta) * tau**0.625 
     *     * EXP (-delta)
     *     - 0.03812936803576 * delta**4.0 * (5.0 - delta) * tau**1.75 
     *     * EXP (-delta)
     *     - 0.15352245383006 * (1.0 - 2.0 * delta**2.0) * tau**3.625 
     *     * EXP (-delta**2.0)
     *     - 0.026726814910919 * delta**3.0 * (4.0 - 2.0 * delta**2.0) 
     *     * tau**3.625 * EXP (-delta**2.0)
     *     - 0.025675298677127 * delta**2.0 * (3.0 - 3.0 * delta**3.0) 
     *     * tau**14.5 * EXP (-delta**3.0)
     *     + 9.5714302123668E-03 * delta**3.0 * (4.0 - 3.0 * delta**3.0) 
     *     * tau**12.0 * EXP (-delta**3.0)


      RETURN
      END


C===============================================================================
      FUNCTION DALPHAR_DDELTA_He(Delta, Tau)
      REAL Delta, Tau
      REAL DALPHAR_DDELTA_O2

      DALPHAR_DDELTA_He = -0.45579024006737
     *     + 1.2516390754925 * tau**0.125
     *     - 1.5438231650621 * tau**0.75
     *     + 0.020467489707221 * 4.0 * delta**3.0 * tau
     *     - 0.34476212380781 * (1.0 - delta) * tau**0.75 * EXP (-delta)
     *     - 0.020858459512787 * delta**2.0 * (3.0 - delta) 
     *     * tau**2.625 * EXP (-delta)
     *     + 0.016227414711778 * delta**4.0 * (5.0 - delta) 
     *     * tau**0.125 * EXP (-delta)
     *     - 0.057471818200892 * delta**4.0 * (5.0 - delta) 
     *     * tau**1.25 * EXP (-delta)
     *     + 0.019462416430715 * delta**4.0 * (5.0 - delta) 
     *     * tau**2.0 * EXP (-delta)
     *     - 0.03329568012302 * delta * (2.0 - 2.0 * delta**2.0) 
     *     * tau * EXP (-delta**2.0)
     *     - 0.010863577372367 * (1.0 - 3.0 * delta**3.0) * tau**4.5 
     *     * EXP (-delta**3.0)
     *     - 0.022173365245954 * delta * (2.0 - 3.0 * delta**3.0) 
     *     * tau**5.0 * EXP (-delta**3.0)

      RETURN
      END


C===============================================================================
      FUNCTION DALPHAR_DDELTA_N2(Delta, Tau)
      REAL Delta, Tau
      REAL DALPHAR_DDELTA_N2

      DALPHAR_DDELTA_N2 = 0.59889711801201 * tau**0.125
     *     - 1.6941557480731 * tau**1.125
     *     + 0.24579736191718 * 2.0 * delta * tau**0.375
     *     - 0.23722456755175 * 2.0 * delta * tau**1.125
     *     + 0.017954918715141 * 4.0 * delta**3.0 * tau**0.625
     *     + 0.014592875720215 * 4.0 * delta**3 * tau**1.5
     *     + 0.10008065936206 * (1.0 - delta) * tau**0.625 
     *     * EXP (-delta)
     *     + 0.73157115385532 * (1.0 - delta) * tau**2.625 
     *     * EXP (-delta)
     *     - 0.88372272336366 * (1.0 - delta) * tau**2.75 * EXP (-delta)
     *     + 0.31887660246708 * delta * (2.0 - delta) * tau**2.125 
     *     * EXP (-delta)
     *     + 0.20766491728799 * delta**2.0 * (3.0 - delta) * tau**2.0 
     *     * EXP (-delta)
     *     - 0.019379315454158 * delta**5.0 * (6.0 - delta) 
     *     * tau**1.75 * EXP (-delta)
     *     - 0.16936641554983 * delta * (2.0 - 2.0 * delta**2.0) 
     *     * tau**4.5 * EXP (-delta**2.0)
     *     + 0.13546846041701 * delta**2.0 * (3.0 - 2.0 * delta**2.0) 
     *     * tau**4.75 * EXP (-delta**2.0)
     *     - 0.033066712095307 * delta**2.0 * (3.0 - 2.0 * delta**2.0) 
     *     * tau**5.0 * EXP (-delta**2.0)
     *     - 0.06069081701855 * delta**3.0 * (4.0 - 2.0 * delta**2.0) 
     *     * tau**4.0 * EXP (-delta**2.0)
     *     + 0.012797548292871 * delta**3.0 * (4.0 - 2.0 * delta**2.0) 
     *     * tau**4.5 * EXP (-delta**2.0)
     *     + 5.8743664107299E-03 * delta * (2.0 - 3.0 * delta**3.0) 
     *     * tau**7.5 * EXP (-delta**3.0)
     *     - 0.018451951971969 * delta**2.0 * (3.0 - 3.0 * delta**3.0) 
     *     * tau**14.0 * EXP (-delta**3.0)
     *     + 4.7226622042472E-03 * delta**3.0 * (4.0 - 3.0 * delta**3.0) 
     *     * tau**11.5 * EXP (-delta**3.0)
     *     - 5.2024079680599E-03 * delta**4.0 * (5.0 - 6.0 * delta**6.0) 
     *     * tau**26.0 * EXP (-delta**6.0)
     *     + 0.043563505956635 * delta**5.0 * (6.0 - 6.0 * delta**6.0) 
     *     * tau**28.0 * EXP (-delta**6.0)
     *     - 0.036251690750939 * delta**5.0 * (6.0 - 6.0 * delta**6.0) 
     *     * tau**30.0 * EXP (-delta**6.0)
     *     - 2.8974026866543E-03 * delta**6.0 * (7.0 - 6.0 * delta**6.0) 
     *     * tau**16.0 * EXP (-delta**6.0)

      RETURN
      END


