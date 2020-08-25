

MODULE SUBSTRATE_CONTRACTION

!USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt

IMPLICIT NONE

CONTAINS



!---Subroutine Substrate Contraction----------------------------------

SUBROUTINE SubstrateDeformation(TS, XP0Ip, YP0Ip, XPIp, YPIp)



	USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XP0Ip, YP0Ip

	REAL(KIND = DP), INTENT(OUT) :: XPIp, YPIp




	IF (TS < TimeStepEquil) THEN

		XPIp = XP0Ip + MaxStrain * XP0Ip
		YPIp = YP0Ip - PoissonRatio * MaxStrain * YP0Ip

	ELSE IF (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt >= 0.0_DP) THEN

		XPIp = XP0Ip + (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt) * XP0Ip
		YPIp = YP0Ip - PoissonRatio * (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt) * YP0Ip

	ELSE

		XPIp = XP0Ip
		YPIp = YP0Ip

	END IF


END SUBROUTINE SubstrateDeformation



!---Subroutine Mapping of Coordinates onto Unstrained Substrate-------


SUBROUTINE Mapping(TS, XI,YI,XX,YY)


	USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XI, YI

	REAL(KIND = DP), INTENT(OUT) :: XX, YY



	IF (TS < TimeStepEquil) THEN

		XX = XI/(1.0_DP + MaxStrain)
		YY = YI/(1.0_DP - PoissonRatio * MaxStrain)

	ELSE IF (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt >= 0.0_DP) THEN

		XX = XI/(1.0_DP + (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt))
		YY = YI/(1.0_DP - PoissonRatio * (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt))

	ELSE

		XX = XI
		YY = YI

	END IF


END SUBROUTINE Mapping



END MODULE SUBSTRATE_CONTRACTION



!---------------------------------------------------------------------------------------------------


MODULE SUBSTRATE_STRETCH

!USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt

IMPLICIT NONE

CONTAINS



!---Subroutine Substrate Contraction----------------------------------

SUBROUTINE SubstrateDeformation(TS, XP0Ip, YP0Ip, XPIp, YPIp)



	USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XP0Ip, YP0Ip

	REAL(KIND = DP), INTENT(OUT) :: XPIp, YPIp




	IF (TS < TimeStepEquil) THEN

		XPIp = XP0Ip
		YPIp = YP0Ip

	ELSE IF (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt >= 0.0_DP) THEN

		XPIp = XP0Ip + StrainRate*(TS - TimeStepEquil)*dt * XP0Ip
		YPIp = YP0Ip - PoissonRatio * StrainRate*(TS - TimeStepEquil)*dt * YP0Ip

	ELSE

		XPIp = XP0Ip + MaxStrain * XP0Ip
		YPIp = YP0Ip - PoissonRatio * MaxStrain * YP0Ip

	END IF


END SUBROUTINE SubstrateDeformation



!---Subroutine Mapping of Coordinates onto Unstrained Substrate-------


SUBROUTINE Mapping(TS, XI,YI,XX,YY)


	USE PARAMETERS, ONLY : Range15, DP, TimeStepEquil, MaxStrain, PoissonRatio, StrainRate, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XI, YI

	REAL(KIND = DP), INTENT(OUT) :: XX, YY



	IF (TS < TimeStepEquil) THEN

		XX = XI
		YY = YI

	ELSE IF (MaxStrain - StrainRate*(TS - TimeStepEquil)*dt >= 0.0_DP) THEN

		XX = XI/(1.0_DP + StrainRate*(TS - TimeStepEquil)*dt)
		YY = YI/(1.0_DP - PoissonRatio * StrainRate*(TS - TimeStepEquil)*dt)

	ELSE

		XX = XI/(1.0_DP + MaxStrain)
		YY = YI/(1.0_DP - PoissonRatio * MaxStrain)

	END IF


END SUBROUTINE Mapping



END MODULE SUBSTRATE_STRETCH


!---------------------------------------------------------------------------------------------------



MODULE SUBSTRATE_PERIODIC_STRETCH

!USE PARAMETERS, ONLY : Range15, DP, pi, TimeStepEquil, MaxStrain, Freq, PoissonRatio, dt

IMPLICIT NONE

CONTAINS



!---Subroutine Substrate Contraction----------------------------------

SUBROUTINE SubstrateDeformation(TS, XP0Ip, YP0Ip, XPIp, YPIp)



	USE PARAMETERS, ONLY : Range15, DP, pi, MaxStrain, Freq, PoissonRatio, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XP0Ip, YP0Ip

	REAL(KIND = DP), INTENT(OUT) :: XPIp, YPIp




	XPIp = XP0Ip + MaxStrain * XP0Ip * (0.5_DP - 0.5_DP * DCOS(2.0_DP*pi*Freq*TS*dt))
	YPIp = YP0Ip - PoissonRatio * MaxStrain * YP0Ip * (0.5_DP - 0.5_DP * DCOS(2.0_DP*pi*Freq*TS*dt))


END SUBROUTINE SubstrateDeformation



!---Subroutine Mapping of Coordinates onto Unstrained Substrate-------


SUBROUTINE Mapping(TS, XI,YI,XX,YY)


	USE PARAMETERS, ONLY : Range15, DP, pi, MaxStrain, Freq, PoissonRatio, dt


	INTEGER(KIND = Range15), INTENT(IN) :: TS

	REAL(KIND = DP), INTENT(IN) :: XI, YI

	REAL(KIND = DP), INTENT(OUT) :: XX, YY




	XX = XI/(1.0_DP + MaxStrain * (0.5_DP - 0.5_DP * DCOS(2.0_DP*pi*Freq*TS*dt)))
	YY = YI/(1.0_DP - PoissonRatio * MaxStrain * (0.5_DP - 0.5_DP * DCOS(2.0_DP*pi*Freq*TS*dt)))


END SUBROUTINE Mapping



END MODULE SUBSTRATE_PERIODIC_STRETCH
