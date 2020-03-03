MODULE PARAMETERS

IMPLICIT NONE


! Unit: um, sec, pN



INTEGER, PARAMETER ::		Range15 = SELECTED_INT_KIND(15), &
				DP = SELECTED_REAL_KIND(14), &
				NumAssay = 1, &
				NumMol = 1, &
				NumBeads = 13.0, &
				NumSpecies = 2, &		! 2Motors
!				MaxAreaNum = 1000, &
				MaxNumIteration=1E6, &
				AreaRenewDiv = 2E4, &
				OutPutDiv=2E4, &
				MaxInteractinNum = 1, & !Maximum Number of MTs which a complex can interact with
				seed = 213 !97845 !445 !78448 !657 !2215 !78448 !7484541 !4397 ! 242, 6990, 52804, 2072544

INTEGER(KIND = Range15), PARAMETER :: NumTimeStep = 10.0E6, TimeStepEquil = 1E9, MaxNumParticles = 1E6, TimeForceON = 1E9

INTEGER, PARAMETER :: MaxAreaNum = NumTimeStep/AreaRenewDiv + 1

REAL(KIND = DP), PARAMETER :: 	pi = 3.14159265358979_DP, &
				kBT = 0.0041_DP, &
				dt=5.0E-7_DP, &
				Tol=1.0E-6_DP, &
				BondLength=0.25_DP, &
				EI=0.073_DP, & !Actin Filament (Gittes et al., JCB 1993)
				Motor_Density = 3000.0_DP, &
				k = 300.0_DP,& 
!
				Stepsize = 1.0E-2_DP, & ! was 1.0E-2
				k_a = 40.0_DP,&
				k_d0 = 350.0_DP,&
				k_t = 2.0_DP,& ! /ATP[uM]
				ATP = 2000.0_DP,& ![uM]
				k_hp = 100.0_DP,&
				k_hm = 10.0_DP,&
				delta_x = -1.86E-3_DP,&
!
				F_Motor_Detach1 = 9.2_DP,&
				F_Motor_Detach2 = 10.0_DP,&
				CaptureRadius=0.020_DP, &
				Mergin4List = BondLength*1.5_DP,&
				track_width = 4.0_DP ,&
!
				Species1Ratio = 1.0_DP, &

!---External Force----------------------------------------------------
				ExtForceDensity0 = 0.0_DP, & !pN/um
				ExtForceMag = 0.0_DP, &
!---Track Surface-----------------------------------------------------
				XLimit = 500.0_DP, &

!---Track Surface-----------------------------------------------------

			!Periodic Stretch & Contraction
!				MaxStrain = 0.2_DP, PoissonRatio = 0.4_DP, &
!				MaxStrain = 0.75_DP, PoissonRatio = 0.2_DP, &
!				MaxStrain = 1.4_DP, PoissonRatio = 0.15_DP, &
				Freq = 1.0_DP, & !(Hz)

			!Contraction
				MaxStrain = 0.0_DP, &
!				MaxStrain = 0.2_DP, &
!				MaxStrain = 0.5_DP, &
!				MaxStrain = 1.3_DP, &
				StrainRate = 0.05_DP, &
				PoissonRatio = 0.5_DP, &

				HorizontalLength = DBLE(NumBeads-1)*BondLength + 1.0_DP, &
				VerticalLength = 0.5*DBLE(NumBeads-1)*BondLength + 1.0_DP



END MODULE PARAMETERS
