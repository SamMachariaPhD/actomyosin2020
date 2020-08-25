PROGRAM BEADRODPOLYMER

! Modified from MotilityAssayMain_v29.f90


USE PARAMETERS
USE mtmod

!USE PLANAR_TRACK_CONFINEMENT
USE X_CHANNEL_CONFINEMENT


!USE UNIFORM_FORCE_X
USE END_DRAG



IMPLICIT NONE

INTEGER	:: seed_Assay, I_Assay, I, J, II, JJ, JJJ, J_Contact, AreaCounter, I_Area, EraseCounter, I_Erase, &
	NumIteration, &
	OutputfileCounter, &
	Iint, CounterBuffer, ISpecies

INTEGER(KIND = Range15):: TS, Ip, NumParticles

INTEGER(KIND = Range15), DIMENSION(MaxAreaNum) :: AddedParticleNum

INTEGER, DIMENSION(MaxNumParticles) :: MotorSpecies		! 2Motors	1=Myosin, 2=Dead Motor (or "Zombie" in Dan's email)

INTEGER, DIMENSION(NumMol*NumBeads, MaxNumParticles) :: CandidateList

INTEGER, DIMENSION(MaxNumParticles) :: Release_ADP		! 2Motors

LOGICAL :: StateConstraint, StatConfinement, StatConfinementParticle

REAL(KIND = DP) ::	IniAngle, &
			UR1, UR2, UR3, &
			NormRand4Particls, &
			Intercept_Contact, X_Contact, Y_Contact, Z_Contact, DisHead2Tail, ProjHead2Tail_MT, &
			gamma_Bead, D_Bead, gamma_Particle, &
			v_Motor, Intercept, &
			XCM, YCM

REAL(KIND = DP), DIMENSION(MaxAreaNum) :: AreaOriginX, AreaOriginY, AreaOriginUx, AreaOriginUy

REAL(KIND = DP), DIMENSION(NumBeads) :: XI, YI, ZI, XI_temp, YI_temp, ZI_temp, &
			FIx, FIy, FIz, &
			NormRandVector4Beads, XINormRandVector4Beads, YINormRandVector4Beads, ZINormRandVector4Beads

REAL(KIND = DP), DIMENSION(MaxNumParticles) :: XP, YP, ZP, &
!			FPx, FPy, FPz, &										!v27
			NormRandVector4Particles, &
			XPBuffer, YPBuffer, ZPBuffer, &
			Elongation

REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles) :: ContactState, TempContact

REAL(KIND = DP), DIMENSION(MaxNumParticles) :: F_Motor_X, F_Motor_Y, F_Motor_Z

REAL(KIND = DP), DIMENSION(NumBeads, NumMol) :: X, Y, Z, Fx, Fy, Fz, Force_Bead_X, Force_Bead_Y, Force_Bead_Z, &
			XNormRandVector4Beads, YNormRandVector4Beads, ZNormRandVector4Beads

CHARACTER(LEN=30) :: OutFileName, OutFileNameF, OutFileNameP



seed_Assay = seed

!---Calculating Parameters

gamma_Bead = 3.0_DP*pi*0.001_DP*BondLength/DLOG(BondLength/0.006_DP)
D_Bead = kBT/gamma_Bead

!gamma_Particle = kBT/D_Particle




!Output chamber boundary in vtk format--------------------------------
!CALL OutputBoundary


OPEN(10,FILE='InitialCondition.txt')


!---Repeat Assay------------------------------------------------------

Loop_Assay: DO I_Assay=1, NumAssay



!---Assay Initial Condition-------------------------------------------


CALL sgrnd(seed_Assay)

TS = 0
NumParticles = 0


AreaCounter = 1
EraseCounter = 0

OutputfileCounter = 0

XP=-100.0_DP
YP=5.0_DP



!UR1 = grnd()
!IniAngle = 2.0_DP*pi*UR1
IniAngle = 0.0_DP*pi



WRITE(10,'(I3,I7,F15.10)') I_Assay, seed_Assay, IniAngle


seed_Assay = seed_Assay+1


!OPEN OUTPUT FILE-----------------------------------------------------


WRITE(OutFileName,'(A,I3.3,A)') 'Conformation_A', I_Assay, '.txt'
OPEN(11,FILE=OutFileName)

WRITE(OutFileName,'(A,I3.3,A)') 'MTPlusEndPositions_A', I_Assay, '.txt'
OPEN(15,FILE=OutFileName)

WRITE(OutFileName,'(A,I3.3,A)') 'TipXY_A', I_Assay, '.txt'
OPEN(16,FILE=OutFileName)

WRITE(OutFileName,'(A,I3.3,A)') 'KinesinTrajectories_A', I_Assay, '.txt'
OPEN(17,FILE=OutFileName)

!WRITE(OutFileName,'(A,I3.3,A)') 'KinesinElongation', I_Assay, '.txt'
!OPEN(18,FILE=OutFileName)

!WRITE(OutFileName,'(A,I3.3,A)') 'KinesinProjectionMT', I_Assay, '.txt'
!OPEN(19,FILE=OutFileName)

WRITE(OutFileName,'(A,I3.3,A)') 'KinesinHeadTail_A', I_Assay, '.txt'
OPEN(20,FILE=OutFileName)

WRITE(OutFileName,'(A,I3.3,A)') 'AllKinesin_A', I_Assay, '.txt'			!v25
OPEN(21,FILE=OutFileName)							!v25

WRITE(OutFileName,'(A,I3.3,A)') 'KinesinForce_A', I_Assay, '.txt'		!v27
OPEN(22,FILE=OutFileName)							!v27

WRITE(OutFileName,'(A,I3.3,A)') 'AreaEraseCounter_A', I_Assay, '.txt'		!v28
OPEN(23,FILE=OutFileName)


!---Initial Conformation of Polymer Chain-----------------------------

DO I=1, NumMol

	X(NumBeads,I) = DCOS(IniAngle)
	Y(NumBeads,I) = DSIN(IniAngle)
	Z = 0.0125_DP

	DO J=NumBeads-1, 1, -1

		X(J,I) = X(J+1,I) + BondLength*DCOS(IniAngle)
		Y(J,I) = Y(J+1,I) + BondLength*DSIN(IniAngle)

	END DO
END DO





!---Initial Particle Location-----------------------------------------

XCM=SUM(X(:,1))/DBLE(NumBeads)		!v11
YCM=SUM(Y(:,1))/DBLE(NumBeads)		!v11

AreaOriginUx(AreaCounter) = DCOS(IniAngle)
AreaOriginUy(AreaCounter) = DSIN(IniAngle)

AreaOriginX(AreaCounter) = XCM		!v11
AreaOriginY(AreaCounter) = YCM		!v11



DO Ip=1, INT(Motor_Density*HorizontalLength*VerticalLength, Range15)		!v11

	UR1 = grnd()
	UR2 = grnd()

	XP(Ip) = XCM + HorizontalLength * (UR1 - 0.5_DP) * AreaOriginUx(AreaCounter) - VerticalLength * (UR2 - 0.5_DP) * AreaOriginUy(AreaCounter)		!v11
	YP(Ip) = YCM + HorizontalLength * (UR1 - 0.5_DP) * AreaOriginUy(AreaCounter) + VerticalLength * (UR2 - 0.5_DP) * AreaOriginUx(AreaCounter)		!v11
	WRITE(21,'(3F15.10)') XP(Ip), YP(Ip)									!v29


	UR1 = grnd()					! 2Motors
	IF (UR1 <= Species1Ratio) THEN			! 2Motors

		MotorSpecies(Ip) = 1			! 2Motors
		UR2 = grnd()				! 2Motors
		IF(UR2 <= 0.091)THEN			! 2Motors
			ContactState(:,Ip) = -1.0_DP	! 2Motors
		ELSE					! 2Motors
			ContactState(:,Ip) = 0.0_DP	! 2Motors
		END IF					! 2Motors

	ELSE						! 2Motors

		MotorSpecies(Ip) = 2			! 2Motors
		ContactState(:,Ip) = 0.0_DP		! 2Motors

	END IF						! 2Motors





END DO

NumParticles = NumParticles + INT(Motor_Density*HorizontalLength*VerticalLength, Range15)		!v11

!AreaCounter = 1
AddedParticleNum(AreaCounter) = INT(Motor_Density*HorizontalLength*VerticalLength, Range15)		!v11





Release_ADP = 0





CandidateList(:,1:NumParticles) = 1
CandidateList(:,(NumParticles+1):MaxNumParticles) = 0





!Filament-Particle Contact--------------------------------------------

TempContact = 0.0_DP															!8/16
CALL CheckingMotorFilamentContact(X, Y, Z, XP, YP, ZP, NumParticles, CandidateList, ContactState, TempContact)		!8/16

!CALL MotorBinding(NumParticles, ContactState, TempContact)								!8/16


!---Motor Binding at Initial Condition--------------------------------

DO Ip=1, NumParticles

	DO II=1, NumMol

		IF ((TempContact(II, Ip) >= 1.0_DP) .AND. (TempContact(II, Ip) <= DBLE(NumBeads))) THEN

			UR1 = grnd()
			IF(UR1 <= 0.1_DP)THEN   ! TS=0 ��1���̃~�I�V���̂݌���
				ContactState(II, Ip) = TempContact(II, Ip)
				TempContact(II, Ip) = 0.0_DP
			END IF

		END IF

	END DO

END DO





!---Calculating Forces upon Beads, Moters and Particles---------------


CALL CalculateForceMotor(X, Y, Z, XP, YP, ZP, ContactState, F_Motor_X, F_Motor_Y, F_Motor_Z, NumParticles, Elongation)						!8/16

DO Ip=1, NumParticles
	CALL MotorForcedDetachment(F_Motor_X(Ip), F_Motor_Y(Ip), F_Motor_Z(Ip), ContactState(:,Ip), Elongation(Ip), MotorSpecies(Ip))			!v9
END DO

CALL CalculateForceBead(ContactState, Force_Bead_X, Force_Bead_Y, Force_Bead_Z, F_Motor_X, F_Motor_Y, F_Motor_Z,NumParticles)				!8/16





!---Calculating Forces upon Beads-------------------------------------

Fx=Force_Bending(X) + Force_Bead_X + ExtF_X(TS, X, Y, Z)
Fy=Force_Bending(Y) + Force_Bead_Y + ExtF_Y(TS, X, Y, Z)
Fz=Force_Bending(Z) + Force_Bead_Z + ExtF_Z(TS, X, Y, Z)





!---Output Polymer Initial Conformation-------------------------------

DO J=1, NumBeads

	WRITE(11,'(I10)', ADVANCE = "NO") TS	!v4

	DO I=1, NumMol-1
		WRITE(11,'(3F25.15)', ADVANCE = "NO") X(J,I), Y(J,I), Z(J,I)
	END DO
	WRITE(11,'(3F25.15)') X(J,NumMol), Y(J,NumMol), Z(J,NumMol)
END DO





!---Output Initial MT Plus End Location-------------------------------

DO I=1, NumMol-1
	WRITE(15,'(3F14.6)', ADVANCE = "NO") X(NumBeads,I), Y(NumBeads,I), Z(NumBeads,I)
END DO
WRITE(15,'(3F14.6)') X(NumBeads,NumMol), Y(NumBeads,NumMol), Z(NumBeads,NumMol)





!---Output Initial Leading Tip XY-Location----------------------------

DO I=1, NumMol-1
	WRITE(16,'(I10, 2F14.6)', ADVANCE = "NO") TS, X(1,I), Y(1,I)	!v10
END DO
WRITE(16,'(I10, 2F14.6)') TS, X(1,NumMol), Y(1,NumMol)			!v10





!---Output Initial Kinesin Location on MT-----------------------------

DO Ip=1, NumParticles-1
	WRITE(17,'(F14.6)', ADVANCE = "NO") ContactState(1,Ip)*BondLength
END DO
WRITE(17,'(F14.6)') ContactState(1,NumParticles)*BondLength





!---Output Initial Kinesin Elongation & Projection--------------------

DO Ip=1, NumParticles

	IF (ContactState(1,Ip) >= 1.0_DP) THEN

		J_Contact = INT(ContactState(1,Ip))
		Intercept_Contact = ContactState(1,Ip) - DBLE(J_Contact)

		IF (J_Contact >= NumBeads) THEN

			X_Contact = X(J_Contact,1) + Intercept_Contact * (X(J_Contact,1) - X(J_Contact-1,1))
			Y_Contact = Y(J_Contact,1) + Intercept_Contact * (Y(J_Contact,1) - Y(J_Contact-1,1))
			Z_Contact = Z(J_Contact,1) + Intercept_Contact * (Z(J_Contact,1) - Z(J_Contact-1,1))

		ELSE

			X_Contact = X(J_Contact,1) + Intercept_Contact * (X(J_Contact+1,1) - X(J_Contact,1))
			Y_Contact = Y(J_Contact,1) + Intercept_Contact * (Y(J_Contact+1,1) - Y(J_Contact,1))
			Z_Contact = Z(J_Contact,1) + Intercept_Contact * (Z(J_Contact+1,1) - Z(J_Contact,1))

		END IF

!		DisHead2Tail = (X_Contact - XP(Ip))*(X_Contact - XP(Ip)) + &
!				(Y_Contact - YP(Ip))*(Y_Contact - YP(Ip)) + &
!				(Z_Contact - ZP(Ip))*(Z_Contact - ZP(Ip))
!		DisHead2Tail = DSQRT(DisHead2Tail)


!		ProjHead2Tail_MT = (X_Contact - XP(Ip))*(X(J_Contact+1,1) - X(J_Contact,1)) + &
!					(Y_Contact - YP(Ip))*(Y(J_Contact+1,1) - Y(J_Contact,1)) + &
!					(Z_Contact - ZP(Ip))*(Z(J_Contact+1,1) - Z(J_Contact,1))
!		ProjHead2Tail_MT = ProjHead2Tail_MT/BondLength

		WRITE(20,'(2I10, 7F14.6, I2)') TS, Ip, ContactState(1,Ip), X_Contact, Y_Contact, Z_Contact, XP(Ip), YP(Ip), ZP(Ip), MotorSpecies(Ip)

!	ELSE

!		DisHead2Tail = 0.0_DP
!		ProjHead2Tail_MT = 0.0_DP

	END IF


!	WRITE(18,'(F14.6)', ADVANCE = "NO") DisHead2Tail
!	WRITE(19,'(F14.6)', ADVANCE = "NO") ProjHead2Tail_MT


END DO
!WRITE(18,*)
!WRITE(19,*)



!---Output Initial Kinesin Force--------------------------------------

DO Ip=1, NumParticles

	IF (ContactState(1,Ip) >= 1.0_DP) THEN

		WRITE(22,'(2I10, 7F14.6)') TS, Ip, ContactState(1,Ip), F_Motor_X(Ip), F_Motor_Y(Ip), F_Motor_Z(Ip)	!v27 I=1

	END IF

END DO



!---Output Area & Erase Counters--------------------------------------

WRITE(23,*) TS, AreaCounter, EraseCounter						!v28



!---Output Polymer Conformation in vtk format-------------------------

	WRITE(OutFileNameF,'(A,I3.3,A,I7.7,A)') 'Filament_A', I_Assay, 'T', OutputfileCounter, '.vtk'

	OPEN(33,FILE=OutFileNameF)

	WRITE(33,'(A)') "# vtk DataFile Version 2.0"
	WRITE(33,'(A, I7)') "Filament: OutputfileCounter = ", OutputfileCounter
	WRITE(33,'(A)') "ASCII"
	WRITE(33,'(A)') "DATASET UNSTRUCTURED_GRID"
	WRITE(33,'(A, I10, A)') "POINTS", NumMol*NumBeads, " float"

	DO I=1, NumMol
		DO J=1, NumBeads
			WRITE(33,'(3F10.5)') X(J,I), Y(J,I), Z(J,I)
		END DO
	END DO

	WRITE(33,'(A, I10, I10)') "CELLS", NumMol, (1+NumBeads)*NumMol
	DO I=1, NumMol
		!WRITE(33,'(I3, I10)') NumBeads, (J*I-1, J=1, NumBeads)
		WRITE(33,*) NumBeads, (J+(I-1)*NumBeads-1, J=1, NumBeads)
	END DO

	WRITE(33,'(A, I10)') "CELL_TYPES ", NumMol
	DO I=1, NumMol
		WRITE(33,'(I3)') 4
	END DO


	CLOSE(33)



!---Output MT Plus End Locations in vtk format------------------------

	WRITE(OutFileNameP,'(A,I3.3,A,I7.7,A)') 'MTPlusEnd_A', I_Assay, 'T', OutputfileCounter, '.vtk'

	OPEN(35,FILE=OutFileNameP)

	WRITE(35,'(A)') "# vtk DataFile Version 2.0"
	WRITE(35,'(A, I7)') "MT Plus Ends: OutputfileCounter = ", OutputfileCounter
	WRITE(35,'(A)') "ASCII"
	WRITE(35,'(A)') "DATASET UNSTRUCTURED_GRID"
	WRITE(35,'(A, I10, A)') "POINTS", NumMol, " float"

	DO I=1, NumMol
		WRITE(35,'(3F10.5)') X(NumBeads,I), Y(NumBeads,I), Z(NumBeads,I)
	END DO

	WRITE(35,'(A, I10, I10)') "CELLS", NumMol, 2*NumMol

	DO I=1, NumMol
		WRITE(35,'(I3, I10)') 1, I-1
	END DO

	WRITE(35,'(A, I10)') "CELL_TYPES ", NumMol
	DO I=1, NumMol
		WRITE(35,'(I3)') 1
	END DO


	CLOSE(35)



!---Output Contact States of particles--------------------------------

!!!---Assuming that NumMol = 1

	WRITE(OutFileName,'(A,I3.3,A,I7.7,A)') 'ContactStates_A', I_Assay, 'T', OutputfileCounter, '.txt'

	OPEN(36,FILE=OutFileName)

	DO Ip=1, NumParticles
		WRITE(36,'(I7.7,I3.1,F10.5)', ADVANCE = "YES") Ip, MotorSpecies(Ip), ContactState(1,Ip)
	END DO

	CLOSE(36)



!---Output Interacting Motor Locations in vtk Format for Each Species

DO ISpecies = 1, NumSpecies

		WRITE(OutFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'IntSpecie', ISpecies, '_A', I_Assay, 'T', OutputfileCounter, '.vtk'

		OPEN(38,FILE=OutFileName)

		CounterBuffer = 0
		DO Ip=1, NumParticles

			IF ((MotorSpecies(Ip)==ISpecies) .AND. (SUM(ContactState(:,Ip))>1)) THEN

				CounterBuffer = CounterBuffer + 1
				XPBuffer(CounterBuffer) = XP(Ip)
				YPBuffer(CounterBuffer) = YP(Ip)
				ZPBuffer(CounterBuffer) = ZP(Ip)

			END IF

		END DO

		WRITE(38,'(A)') "# vtk DataFile Version 2.0"
		WRITE(38,'(A, I7)') "Particle: OutputfileCounter = ", OutputfileCounter
		WRITE(38,'(A)') "ASCII"
		WRITE(38,'(A)') "DATASET UNSTRUCTURED_GRID"
		WRITE(38,'(A, I10, A)') "POINTS", CounterBuffer, " float"

		DO Ip=1, CounterBuffer
			WRITE(38,'(3F10.5)') XPBuffer(Ip), YPBuffer(Ip), ZPBuffer(Ip)
		END DO

		WRITE(38,'(A, I10, I10)') "CELLS", CounterBuffer, 2*CounterBuffer

		DO Ip=1, CounterBuffer
			WRITE(38,'(I3, I10)') 1, Ip-1
		END DO

		WRITE(38,'(A, I10)') "CELL_TYPES ", CounterBuffer
		DO Ip=1, CounterBuffer
			WRITE(38,'(I3)') 1
		END DO


		CLOSE(38)

	END DO



!---Output Particle Locations in vtk format for each species


	DO ISpecies = 1, NumSpecies

		WRITE(OutFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'MotorSpecie', ISpecies, '_A', I_Assay, 'T', OutputfileCounter, '.vtk'

		OPEN(38,FILE=OutFileName)

		CounterBuffer = 0
		DO Ip=1, NumParticles

			IF (MotorSpecies(Ip)==ISpecies) THEN

				CounterBuffer = CounterBuffer + 1
				XPBuffer(CounterBuffer) = XP(Ip)
				YPBuffer(CounterBuffer) = YP(Ip)
				ZPBuffer(CounterBuffer) = ZP(Ip)

			END IF

		END DO

		WRITE(38,'(A)') "# vtk DataFile Version 2.0"
		WRITE(38,'(A, I7)') "Particle: OutputfileCounter = ", OutputfileCounter
		WRITE(38,'(A)') "ASCII"
		WRITE(38,'(A)') "DATASET UNSTRUCTURED_GRID"
		WRITE(38,'(A, I10, A)') "POINTS", CounterBuffer, " float"

		DO Ip=1, CounterBuffer
			WRITE(38,'(3F10.5)') XPBuffer(Ip), YPBuffer(Ip), ZPBuffer(Ip)
		END DO

		WRITE(38,'(A, I10, I10)') "CELLS", CounterBuffer, 2*CounterBuffer

		DO Ip=1, CounterBuffer
			WRITE(38,'(I3, I10)') 1, Ip-1
		END DO

		WRITE(38,'(A, I10)') "CELL_TYPES ", CounterBuffer
		DO Ip=1, CounterBuffer
			WRITE(38,'(I3)') 1
		END DO


		CLOSE(38)

	END DO


!---------------------------------------------------------------------





!---Time Evolution----------------------------------------------------

Loop_Time: DO TS=1, NumTimeStep





!---Generating RND Num Matrix-----------------------------------------

DO I=1,NumMol
	CALL NormRNDVector(NumBeads, NormRandVector4Beads)

	DO J=1,NumBeads
		XNormRandVector4Beads(J,I) = NormRandVector4Beads(J)
	END DO
END DO

DO I=1,NumMol
	CALL NormRNDVector(NumBeads, NormRandVector4Beads)

	DO J=1,NumBeads
		YNormRandVector4Beads(J,I) = NormRandVector4Beads(J)
	END DO
END DO

DO I=1,NumMol
	CALL NormRNDVector(NumBeads, NormRandVector4Beads)

	DO J=1,NumBeads
		ZNormRandVector4Beads(J,I) = NormRandVector4Beads(J)
	END DO
END DO





Loop_Mol: DO I=1, NumMol


	DO J=1, NumBeads
		XI(J)=X(J,I)
	END DO
	DO J=1, NumBeads
		YI(J)=Y(J,I)
	END DO
	DO J=1, NumBeads
		ZI(J)=Z(J,I)
	END DO

	DO J=1, NumBeads
		FIx(J)=Fx(J,I)
	END DO
	DO J=1, NumBeads
		FIy(J)=Fy(J,I)
	END DO
	DO J=1, NumBeads
		FIz(J)=Fz(J,I)
	END DO

	DO J=1, NumBeads
		XINormRandVector4Beads(J)=XNormRandVector4Beads(J,I)
	END DO
	DO J=1, NumBeads
		YINormRandVector4Beads(J)=YNormRandVector4Beads(J,I)
	END DO
	DO J=1, NumBeads
		ZINormRandVector4Beads(J)=ZNormRandVector4Beads(J,I)
	END DO


!---Unconstrained Movements of Beads


	CALL NormRNDVector(NumBeads, NormRandVector4Beads)
	XI_temp = XI + FIx/gamma_Bead*dt + XINormRandVector4Beads * DSQRT(2.0_DP * D_Bead * dt)

	CALL NormRNDVector(NumBeads, NormRandVector4Beads)
	YI_temp = YI + FIy/gamma_Bead*dt + YINormRandVector4Beads * DSQRT(2.0_DP * D_Bead * dt)

	CALL NormRNDVector(NumBeads, NormRandVector4Beads)
	ZI_temp = ZI + FIz/gamma_Bead*dt + ZINormRandVector4Beads * DSQRT(2.0_DP * D_Bead * dt)



!---Iteration until Constraint and Confinement are Statisfied

	NumIteration = 0
	StateConstraint = .FALSE.
	StatConfinement = .FALSE.

	IterationConstraintConfinement: DO

!		IF ((StateConstraint) .OR. (NumIteration > MaxNumIteration)) EXIT IterationConstraintConfinement

		IF (StateConstraint .AND. StatConfinement) EXIT IterationConstraintConfinement

		IF (NumIteration > MaxNumIteration) THEN
			WRITE(*,*) "Too many iterations!"
			EXIT IterationConstraintConfinement
		END IF

		NumIteration = NumIteration + 1



!---Confinement (Filament)--------------------------------------------


		StatConfinement = .TRUE.


		Loop_Bead_Confinement: DO J=1, NumBeads

			CALL Confinement(XI_temp(J), YI_temp(J), ZI_temp(J), XI(J), YI(J), ZI(J), StatConfinement)

		END DO Loop_Bead_Confinement



!---Constraint--------------------------------------------------------

		StateConstraint = .TRUE.

!		DX = 0.0_DP
!		DY = 0.0_DP
!		DZ = 0.0_DP

		Loop_Bond_Constraint: DO J=1, NumBeads-1

			CALL Constraint(XI_temp(J), YI_temp(J), ZI_temp(J), XI_temp(J+1), YI_temp(J+1), ZI_temp(J+1), &
					XI(J), YI(J), ZI(J), XI(J+1), YI(J+1), ZI(J+1), StateConstraint)

		END DO Loop_Bond_Constraint



!---------------------------------------------------------------------

	END DO IterationConstraintConfinement



	DO J=1,NumBeads
		X(J,I) = XI_temp(J)
	END DO
	DO J=1,NumBeads
		Y(J,I) = YI_temp(J)
	END DO
	DO J=1,NumBeads
		Z(J,I) = ZI_temp(J)
	END DO


END DO Loop_Mol





!---Motor State Conversion--------------------------------------------

DO Ip=1, NumParticles															! 2Motors
	IF (MotorSpecies(Ip) == 1) THEN													! 2Motors
		CALL MotorStateConv(X, Y, Z, F_Motor_X(Ip), F_Motor_Y(Ip), F_Motor_Z(Ip), ContactState(:,Ip), Release_ADP(Ip))		! 2Motors
	END IF																! 2Motors
END DO																	! 2Motors





!---Renew Motor Population--------------------------------------------

IF (MOD(TS,AreaRenewDiv)==0) THEN

	CALL RenewMotorPopulation(X, Y, XP, YP, ZP, ContactState, CandidateList, AddedParticleNum, AreaCounter, EraseCounter, AreaOriginX, AreaOriginY, AreaOriginUx, AreaOriginUy, NumParticles, MotorSpecies)

END IF





!Filament-Particle Contact--------------------------------------------

TempContact = 0.0_DP
CALL CheckingMotorFilamentContact(X, Y, Z, XP, YP, ZP, NumParticles, CandidateList, ContactState, TempContact)		!8/16
CALL MotorBinding(NumParticles, ContactState, TempContact)								!8/16

DO Ip=1, NumParticles													! 2Motors
	IF (MotorSpecies(Ip) == 1) THEN											! 2Motors
		CALL MotorMove(ContactState(:,Ip), TempContact(:,Ip))							! 2Motors
	ELSE IF (MotorSpecies(Ip) == 2) THEN                                ! v07 on 2020/01/22
        CALL MotorStuck(ContactState(:,Ip), TempContact(:,Ip))
    END IF														! 2Motors
END DO															! 2Motors




!---Calculating Forces upon Beads, Moters and Particles---------------


CALL CalculateForceMotor(X, Y, Z, XP, YP, ZP, ContactState, F_Motor_X, F_Motor_Y, F_Motor_Z, NumParticles, Elongation)						!8/16

DO Ip=1, NumParticles
	CALL MotorForcedDetachment(F_Motor_X(Ip), F_Motor_Y(Ip), F_Motor_Z(Ip), ContactState(:,Ip), Elongation(Ip), MotorSpecies(Ip))			!v9
END DO

CALL CalculateForceBead(ContactState, Force_Bead_X, Force_Bead_Y, Force_Bead_Z, F_Motor_X, F_Motor_Y, F_Motor_Z,NumParticles)				!8/16





!---Calculating Forces upon Beads-------------------------------------

Fx=Force_Bending(X) + Force_Bead_X + ExtF_X(TS, X, Y, Z)
Fy=Force_Bending(Y) + Force_Bead_Y + ExtF_Y(TS, X, Y, Z)
Fz=Force_Bending(Z) + Force_Bead_Z + ExtF_Z(TS, X, Y, Z)





IF (MOD(TS,OutPutDiv)==0) THEN
!IF ((TS >= 5E6) .AND. (TS <= 6E6) .AND. (MOD(TS,OutPutDiv)==0)) THEN

	OutputfileCounter = OutputfileCounter + 1



!---Output Polymer Conformation---------------------------------------

	DO J=1, NumBeads

		WRITE(11,'(I10)', ADVANCE = "NO") TS	!v4

		DO I=1, NumMol-1
			WRITE(11,'(3F25.15)', ADVANCE = "NO") X(J,I), Y(J,I), Z(J,I)
		END DO
		WRITE(11,'(3F25.15)') X(J,NumMol), Y(J,NumMol), Z(J,NumMol)
	END DO




	DO I=1, NumMol-1
		WRITE(15,'(3F14.6)', ADVANCE = "NO") X(NumBeads,I), Y(NumBeads,I), Z(NumBeads,I)
	END DO
	WRITE(15,'(3F14.6)') X(NumBeads,NumMol), Y(NumBeads,NumMol), Z(NumBeads,NumMol)




!---Output Leading Tip XY-Location------------------------------------

	DO I=1, NumMol-1
		WRITE(16,'(I10, 2F14.6)', ADVANCE = "NO") TS, X(1,I), Y(1,I)	!v10
	END DO
	WRITE(16,'(I10, 2F14.6)') TS, X(1,NumMol), Y(1,NumMol)			!v10





!---Output Initial Kinesin Location on MT-----------------------------

	DO Ip=1, NumParticles-1
		WRITE(17,'(F14.6)', ADVANCE = "NO") ContactState(1,Ip)*BondLength
	END DO
	WRITE(17,'(F14.6)') ContactState(1,NumParticles)*BondLength





!---Output Kinesin Elongation & Projection----------------------------

DO Ip=1, NumParticles

	IF (ContactState(1,Ip) >= 1.0_DP) THEN

		J_Contact = INT(ContactState(1,Ip))
		Intercept_Contact = ContactState(1,Ip) - DBLE(J_Contact)

		IF (J_Contact >= NumBeads) THEN

			X_Contact = X(J_Contact,1) + Intercept_Contact * (X(J_Contact,1) - X(J_Contact-1,1))
			Y_Contact = Y(J_Contact,1) + Intercept_Contact * (Y(J_Contact,1) - Y(J_Contact-1,1))
			Z_Contact = Z(J_Contact,1) + Intercept_Contact * (Z(J_Contact,1) - Z(J_Contact-1,1))

		ELSE

			X_Contact = X(J_Contact,1) + Intercept_Contact * (X(J_Contact+1,1) - X(J_Contact,1))
			Y_Contact = Y(J_Contact,1) + Intercept_Contact * (Y(J_Contact+1,1) - Y(J_Contact,1))
			Z_Contact = Z(J_Contact,1) + Intercept_Contact * (Z(J_Contact+1,1) - Z(J_Contact,1))

		END IF

!		DisHead2Tail = (X_Contact - XP(Ip))*(X_Contact - XP(Ip)) + &
!				(Y_Contact - YP(Ip))*(Y_Contact - YP(Ip)) + &
!				(Z_Contact - ZP(Ip))*(Z_Contact - ZP(Ip))
!		DisHead2Tail = DSQRT(DisHead2Tail)


!		ProjHead2Tail_MT = (X_Contact - XP(Ip))*(X(J_Contact+1,1) - X(J_Contact,1)) + &
!					(Y_Contact - YP(Ip))*(Y(J_Contact+1,1) - Y(J_Contact,1)) + &
!					(Z_Contact - ZP(Ip))*(Z(J_Contact+1,1) - Z(J_Contact,1))

		WRITE(20,'(2I10, 7F14.6, I2)') TS, Ip, ContactState(1,Ip), X_Contact, Y_Contact, Z_Contact, XP(Ip), YP(Ip), ZP(Ip), MotorSpecies(Ip)

!	ELSE

!		DisHead2Tail = 0.0_DP
!		ProjHead2Tail_MT = 0.0_DP

	END IF


!	WRITE(18,'(F14.6)', ADVANCE = "NO") DisHead2Tail
!	WRITE(19,'(F14.6)', ADVANCE = "NO") ProjHead2Tail_MT


END DO
!WRITE(18,*)
!WRITE(19,*)



!---Output Initial Kinesin Force--------------------------------------

DO Ip=1, NumParticles

	IF (ContactState(1,Ip) >= 1.0_DP) THEN

		WRITE(22,'(2I10, 7F14.6)') TS, Ip, ContactState(1,Ip), F_Motor_X(Ip), F_Motor_Y(Ip), F_Motor_Z(Ip)	!v27 I=1

	END IF

END DO



!---Output Area & Erase Counters--------------------------------------

WRITE(23,*) TS, AreaCounter, EraseCounter						!v28



!---Output Polymer Conformation in vtk format-------------------------

	WRITE(OutFileNameF,'(A,I3.3,A,I7.7,A)') 'Filament_A', I_Assay, 'T', OutputfileCounter, '.vtk'

	OPEN(33,FILE=OutFileNameF)

	WRITE(33,'(A)') "# vtk DataFile Version 2.0"
	WRITE(33,'(A, I7)') "Filament: OutputfileCounter = ", OutputfileCounter
	WRITE(33,'(A)') "ASCII"
	WRITE(33,'(A)') "DATASET UNSTRUCTURED_GRID"
	WRITE(33,'(A, I10, A)') "POINTS", NumMol*NumBeads, " float"

	DO I=1, NumMol
		DO J=1, NumBeads
			WRITE(33,'(3F10.5)') X(J,I), Y(J,I), Z(J,I)
		END DO
	END DO

	WRITE(33,'(A, I10, I10)') "CELLS", NumMol, (1+NumBeads)*NumMol
	DO I=1, NumMol
		!WRITE(33,'(I3, I10)') NumBeads, (J*I-1, J=1, NumBeads)
		WRITE(33,*) NumBeads, (J+(I-1)*NumBeads-1, J=1, NumBeads)
	END DO

	WRITE(33,'(A, I10)') "CELL_TYPES ", NumMol
	DO I=1, NumMol
		WRITE(33,'(I3)') 4
	END DO


	CLOSE(33)



!---Output MT Plus End Locations in vtk format------------------------

	WRITE(OutFileNameP,'(A,I3.3,A,I7.7,A)') 'MTPlusEnd_A', I_Assay, 'T', OutputfileCounter, '.vtk'

	OPEN(35,FILE=OutFileNameP)

	WRITE(35,'(A)') "# vtk DataFile Version 2.0"
	WRITE(35,'(A, I7)') "MT Plus Ends: OutputfileCounter = ", OutputfileCounter
	WRITE(35,'(A)') "ASCII"
	WRITE(35,'(A)') "DATASET UNSTRUCTURED_GRID"
	WRITE(35,'(A, I10, A)') "POINTS", NumMol, " float"

	DO I=1, NumMol
		WRITE(35,'(3F10.5)') X(NumBeads,I), Y(NumBeads,I), Z(NumBeads,I)
	END DO

	WRITE(35,'(A, I10, I10)') "CELLS", NumMol, 2*NumMol

	DO I=1, NumMol
		WRITE(35,'(I3, I10)') 1, I-1
	END DO

	WRITE(35,'(A, I10)') "CELL_TYPES ", NumMol
	DO I=1, NumMol
		WRITE(35,'(I3)') 1
	END DO


	CLOSE(35)



!---Output Contact States of particles--------------------------------

!!!---Assuming that NumMol = 1

	WRITE(OutFileName,'(A,I3.3,A,I7.7,A)') 'ContactStates_A', I_Assay, 'T', OutputfileCounter, '.txt'

	OPEN(36,FILE=OutFileName)

	DO Ip=1, NumParticles
		WRITE(36,'(I7.7,I3.1,F10.5)', ADVANCE = "YES") Ip, MotorSpecies(Ip), ContactState(1,Ip)
	END DO

	CLOSE(36)



!---Output Interacting Motor Locations in vtk Format for Each Species

DO ISpecies = 1, NumSpecies

		WRITE(OutFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'IntSpecie', ISpecies, '_A', I_Assay, 'T', OutputfileCounter, '.vtk'

		OPEN(38,FILE=OutFileName)

		CounterBuffer = 0
		DO Ip=1, NumParticles

			IF ((MotorSpecies(Ip)==ISpecies) .AND. (SUM(ContactState(:,Ip))>1)) THEN

				CounterBuffer = CounterBuffer + 1
				XPBuffer(CounterBuffer) = XP(Ip)
				YPBuffer(CounterBuffer) = YP(Ip)
				ZPBuffer(CounterBuffer) = ZP(Ip)

			END IF

		END DO

		WRITE(38,'(A)') "# vtk DataFile Version 2.0"
		WRITE(38,'(A, I7)') "Particle: OutputfileCounter = ", OutputfileCounter
		WRITE(38,'(A)') "ASCII"
		WRITE(38,'(A)') "DATASET UNSTRUCTURED_GRID"
		WRITE(38,'(A, I10, A)') "POINTS", CounterBuffer, " float"

		DO Ip=1, CounterBuffer
			WRITE(38,'(3F10.5)') XPBuffer(Ip), YPBuffer(Ip), ZPBuffer(Ip)
		END DO

		WRITE(38,'(A, I10, I10)') "CELLS", CounterBuffer, 2*CounterBuffer

		DO Ip=1, CounterBuffer
			WRITE(38,'(I3, I10)') 1, Ip-1
		END DO

		WRITE(38,'(A, I10)') "CELL_TYPES ", CounterBuffer
		DO Ip=1, CounterBuffer
			WRITE(38,'(I3)') 1
		END DO


		CLOSE(38)

	END DO



!---Output Particle Locations in vtk format for each species


	DO ISpecies = 1, NumSpecies

		WRITE(OutFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'MotorSpecie', ISpecies, '_A', I_Assay, 'T', OutputfileCounter, '.vtk'

		OPEN(38,FILE=OutFileName)

		CounterBuffer = 0
		DO Ip=1, NumParticles

			IF (MotorSpecies(Ip)==ISpecies) THEN

				CounterBuffer = CounterBuffer + 1
				XPBuffer(CounterBuffer) = XP(Ip)
				YPBuffer(CounterBuffer) = YP(Ip)
				ZPBuffer(CounterBuffer) = ZP(Ip)

			END IF

		END DO

		WRITE(38,'(A)') "# vtk DataFile Version 2.0"
		WRITE(38,'(A, I7)') "Particle: OutputfileCounter = ", OutputfileCounter
		WRITE(38,'(A)') "ASCII"
		WRITE(38,'(A)') "DATASET UNSTRUCTURED_GRID"
		WRITE(38,'(A, I10, A)') "POINTS", CounterBuffer, " float"

		DO Ip=1, CounterBuffer
			WRITE(38,'(3F10.5)') XPBuffer(Ip), YPBuffer(Ip), ZPBuffer(Ip)
		END DO

		WRITE(38,'(A, I10, I10)') "CELLS", CounterBuffer, 2*CounterBuffer

		DO Ip=1, CounterBuffer
			WRITE(38,'(I3, I10)') 1, Ip-1
		END DO

		WRITE(38,'(A, I10)') "CELL_TYPES ", CounterBuffer
		DO Ip=1, CounterBuffer
			WRITE(38,'(I3)') 1
		END DO


		CLOSE(38)

	END DO


!---------------------------------------------------------------------


END IF


IF (X(1,1) > XLimit) EXIT Loop_Time	!v12


END DO Loop_Time



!---Output All Particle Locations in vtk format-----------------------


!	WRITE(OutFileNameP,'(A,I7.7,A)') 'Particle', OutputfileCounter, '.vtk'

!	OPEN(21,FILE=OutFileNameP)

!	WRITE(21,'(A)') "# vtk DataFile Version 2.0"
!	WRITE(21,'(A, I7)') "Particle: OutputfileCounter = ", OutputfileCounter
!	WRITE(21,'(A)') "ASCII"
!	WRITE(21,'(A)') "DATASET UNSTRUCTURED_GRID"
!	WRITE(21,'(A, I10, A)') "POINTS", NumParticles, " float"

!	DO Ip=1, NumParticles
!		WRITE(21,'(3F10.5)') XP0(Ip), YP0(Ip), ZP(Ip)
!	END DO

!	WRITE(21,'(A, I10, I10)') "CELLS", NumParticles, 2*NumParticles

!	DO Ip=1, NumParticles
!		WRITE(21,'(I3, I10)') 1, Ip-1
!	END DO

!	WRITE(21,'(A, I10)') "CELL_TYPES ", NumParticles
!	DO Ip=1, NumParticles
!		WRITE(21,'(I3)') 1
!	END DO


!	CLOSE(21)





!CLOSE(10)
CLOSE(11)
CLOSE(15)
CLOSE(16)
CLOSE(17)
!CLOSE(18)
!CLOSE(19)
CLOSE(20)











END DO Loop_Assay

CLOSE(10)




CONTAINS



!---Subroutine NormRNDNum---------------------------------------------
!Generate Normal Random Number
!---------------------------------------------------------------------

SUBROUTINE NormRNDNum(NRN)

	USE PARAMETERS, ONLY : DP
	USE mtmod, ONLY : grnd

	REAL(KIND = DP), INTENT(OUT) :: NRN

	REAL(KIND = DP) :: UR1, UR2



	UR1 = grnd()
	IF (UR1 == 0.0_DP) UR1 = grnd()
	UR2 = grnd()

	NRN = DSQRT(-2.0_DP*DLOG(UR1)) * DCOS(2.0_DP*pi*(UR2))


END SUBROUTINE NormRNDNum



!---Subroutine NormRNDVector-----------------------------------------
!Generate N-Dimensinal Normal Random Vector
!---------------------------------------------------------------------

SUBROUTINE NormRNDVector(N_Dim, NRV)

	USE PARAMETERS, ONLY : DP
	USE mtmod, ONLY : grnd

	INTEGER, INTENT(IN) :: N_Dim
	REAL(KIND = DP), DIMENSION(N_Dim), INTENT(OUT) :: NRV

	INTEGER :: I
	REAL(KIND = DP) :: UR1, UR2, NR



	DO I=1, N_Dim

		UR1 = grnd()
		IF (UR1 == 0.0_DP) UR1 = grnd()
		UR2 = grnd()

		NR = DSQRT(-2.0_DP*DLOG(UR1)) * DCOS(2.0_DP*pi*(UR2))

		NRV(I) = NR

	END DO

END SUBROUTINE NormRNDVector





!---Subroutine MakeList---------------------------------------------

SUBROUTINE MakeList(X, Y, Z, XP, YP, ZP, ContactState, CandidateList)

	USE PARAMETERS, ONLY : DP, MaxNumParticles, NumMol, NumBeads, BondLength, Mergin4List


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(IN) :: XP, YP, ZP

	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(IN) :: ContactState

	INTEGER, DIMENSION(NumMol*NumBeads, MaxNumParticles), INTENT(OUT) :: CandidateList

	REAL(KIND = DP) :: Intercept, SqDistance



	CandidateList = 0



	DO Ip=1, NumParticles

		Loop_Filament_MakeList: DO II=1, NumMol

			IF (ContactState(II, Ip) >= 1.0_DP) CYCLE Loop_Filament_MakeList

			DO JJ=1, NumBeads-1

				Intercept = (XP(Ip) - X(JJ,II))*(X(JJ+1,II) - X(JJ,II)) + &
					(YP(Ip) - Y(JJ,II))*(Y(JJ+1,II) - Y(JJ,II)) + &
					(ZP(Ip) - Z(JJ,II))*(Z(JJ+1,II) - Z(JJ,II))
				Intercept = Intercept/BondLength

				IF ((Intercept > BondLength + Mergin4List) .OR. (Intercept < 0.0_DP - Mergin4List)) CYCLE

				SqDistance = (XP(Ip) - X(JJ,II))*(XP(Ip) - X(JJ,II)) + &
					(YP(Ip) - Y(JJ,II))*(YP(Ip) - Y(JJ,II)) + &
					(ZP(Ip) - Z(JJ,II))*(ZP(Ip) - Z(JJ,II)) - Intercept**2

				IF (SqDistance <= Mergin4List**2) THEN

					CandidateList(NumBeads*(II-1)+JJ, Ip) = 1

				END IF

			END DO

		END DO Loop_Filament_MakeList

	END DO



END SUBROUTINE MakeList




!---Subroutine CheckingMotorFilamentContact-------------------------------------

SUBROUTINE CheckingMotorFilamentContact(X, Y, Z, XP, YP, ZP, NumParticles, CandidateList, ContactState, TempContact)

	USE PARAMETERS, ONLY : DP, MaxNumParticles, NumMol, NumBeads, BondLength, Mergin4List,Stepsize


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(IN) :: XP, YP, ZP

	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(INOUT) :: ContactState, TempContact

	INTEGER(KIND = Range15), INTENT(IN):: NumParticles

	INTEGER, DIMENSION(NumMol*NumBeads, MaxNumParticles), INTENT(IN) :: CandidateList

	REAL(KIND = DP) :: Intercept, SqDistance

	INTEGER(KIND = Range15) :: Ip
    
	INTEGER :: II, JJ



!!!---Assuming that NumMol = 1


!$omp parallel
!$omp do private(JJ,Intercept,SqDistance)
	Loop_Particle_ContactCheck: DO Ip=1, NumParticles

		IF (SUM(CandidateList(:,Ip))==0) CYCLE Loop_Particle_ContactCheck

		IF (ContactState(1, Ip) >= 1.0_DP) CYCLE Loop_Particle_ContactCheck

		IF (ContactState(1, Ip) == -1.0_DP) CYCLE Loop_Particle_ContactCheck


		DO JJ=1, NumBeads-1

			IF (CandidateList(JJ, Ip)==1) THEN

				Intercept = (XP(Ip) - X(JJ,1))*(X(JJ+1,1) - X(JJ,1)) + &
							(YP(Ip) - Y(JJ,1))*(Y(JJ+1,1) - Y(JJ,1)) + &
							(ZP(Ip) - Z(JJ,1))*(Z(JJ+1,1) - Z(JJ,1))
				Intercept = Intercept/BondLength


				IF ((Intercept > BondLength) .OR. (Intercept < 0.0_DP)) CYCLE


				SqDistance = (XP(Ip) - X(JJ,1))*(XP(Ip) - X(JJ,1)) + &
							(YP(Ip) - Y(JJ,1))*(YP(Ip) - Y(JJ,1)) + &
							(ZP(Ip) - Z(JJ,1))*(ZP(Ip) - Z(JJ,1)) - Intercept**2

				IF (SqDistance <= CaptureRadius**2) THEN
					TempContact(1, Ip) = DBLE(JJ) + Intercept/BondLength
				END IF

			END IF

		END DO


	END DO Loop_Particle_ContactCheck
!$omp end do
!$omp end parallel


END SUBROUTINE CheckingMotorFilamentContact



!--Motor Binding----------------------------------------------------------------

SUBROUTINE MotorBinding(NumParticles, ContactState, TempContact)

	USE PARAMETERS , ONLY : DP, NumMol, dt, k_a

	USE mtmod, ONLY : grnd
    
	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(INOUT) :: ContactState, TempContact

	INTEGER(KIND = Range15), INTENT(IN):: NumParticles

	REAL(KIND = DP), DIMENSION(MaxNumParticles) :: URArray
    
	INTEGER(KIND = Range15) :: Ip
    
	INTEGER :: II



	DO Ip=1, NumParticles
		URArray(Ip) = grnd()
	END DO



	Loop_Motor_ContactCheck: DO Ip=1, NumParticles

		Loop_Filament_ContactCheck: DO II=1, NumMol

			IF ((TempContact(II, Ip) >= 1.0_DP ) .AND. (INT(ContactState(II, Ip)) == 0)) THEN

				IF(URArray(Ip) > k_a*dt)THEN  
                			TempContact(II, Ip) = 0.0_DP
				END IF

			END IF

		END DO Loop_Filament_ContactCheck

	END DO Loop_Motor_ContactCheck



END SUBROUTINE MotorBinding



!--Motor Move-------------------------------------------------------------------

SUBROUTINE MotorMove(ContactState_Ip, TempContact_Ip)				! 2Motors

	USE PARAMETERS , ONLY : DP, NumMol, dt, k_a

	USE mtmod, ONLY : grnd

	REAL(KIND = DP), DIMENSION(NumMol), INTENT(INOUT) :: ContactState_Ip, TempContact_Ip


	INTEGER :: II



	Loop_Filament_ContactCheck: DO II=1, NumMol

		IF (TempContact_Ip(II) >= 1.0_DP )THEN

               		ContactState_Ip(II) = TempContact_Ip(II) + Stepsize/BondLength

		END IF

	END DO Loop_Filament_ContactCheck



END SUBROUTINE MotorMove



!--Motor Move-------------------------------------------------------------------

SUBROUTINE MotorStuck(ContactState_Ip, TempContact_Ip)				! 2Motors

	USE PARAMETERS , ONLY : DP, NumMol, dt, k_a

	USE mtmod, ONLY : grnd

	REAL(KIND = DP), DIMENSION(NumMol), INTENT(INOUT) :: ContactState_Ip, TempContact_Ip


	INTEGER :: II



	Loop_Filament_ContactCheck: DO II=1, NumMol

		IF (TempContact_Ip(II) >= 1.0_DP )THEN

               		ContactState_Ip(II) = TempContact_Ip(II)

		END IF

	END DO Loop_Filament_ContactCheck



END SUBROUTINE MotorStuck



!---Subroutine Calculate Force Actin on Motor-----------------------------------

SUBROUTINE CalculateForceMotor(X, Y, Z, XP, YP, ZP, ContactState, F_Motor_X, F_Motor_Y, F_Motor_Z, NumParticles, Elongation)


	USE PARAMETERS, ONLY : DP, MaxNumParticles, NumMol, CaptureRadius, k



	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(IN) :: XP, YP, ZP

	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(IN) :: ContactState

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(OUT) :: F_Motor_X, F_Motor_Y, F_Motor_Z

	INTEGER(KIND = Range15), INTENT(IN):: NumParticles

	REAL(KIND = DP) :: Intercept_Contact, X_Contact, Y_Contact, Z_Contact

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(OUT) :: Elongation

	INTEGER(KIND = Range15) :: Ip

	INTEGER :: II, J_Contact



	F_Motor_X = 0.0_DP
	F_Motor_Y = 0.0_DP
	F_Motor_Z = 0.0_DP


!$omp parallel
!$omp do private(J_Contact,Intercept_Contact,X_Contact,Y_Contact,Z_Contact)
	DO Ip=1, NumParticles
		IF (ContactState(1,Ip) >= 1.0_DP) THEN

			J_Contact = INT(ContactState(1,Ip))
			Intercept_Contact = ContactState(1,Ip) - DBLE(J_Contact)

			IF (J_Contact >= NumBeads) THEN
				X_Contact = X(NumBeads,1) + Intercept_Contact * (X(NumBeads,1) - X(NumBeads-1,1))
				Y_Contact = Y(NumBeads,1) + Intercept_Contact * (Y(NumBeads,1) - Y(NumBeads-1,1))
				Z_Contact = Z(NumBeads,1) + Intercept_Contact * (Z(NumBeads,1) - Z(NumBeads-1,1))
			ELSE
				X_Contact = X(J_Contact,1) + Intercept_Contact * (X(J_Contact+1,1) - X(J_Contact,1))
				Y_Contact = Y(J_Contact,1) + Intercept_Contact * (Y(J_Contact+1,1) - Y(J_Contact,1))
				Z_Contact = Z(J_Contact,1) + Intercept_Contact * (Z(J_Contact+1,1) - Z(J_Contact,1))
			END IF

			Elongation(Ip) = (XP(Ip) - X_Contact)**2 + (YP(Ip) - Y_Contact)**2 + (ZP(Ip) - Z_Contact)**2
			Elongation(Ip) = DSQRT(Elongation(Ip))

			F_Motor_X(Ip) = k * Elongation(Ip) * (XP(Ip) - X_Contact)/Elongation(Ip)
			F_Motor_Y(Ip) = k * Elongation(Ip) * (YP(Ip) - Y_Contact)/Elongation(Ip)
			F_Motor_Z(Ip) = k * Elongation(Ip) * (ZP(Ip) - Z_Contact)/Elongation(Ip)

		END IF
	END DO
!$omp end do
!$omp end parallel



END SUBROUTINE CalculateForceMotor



!---Subroutine Forced Detachment of Motor---------------------------------------

SUBROUTINE MotorForcedDetachment(F_Motor_X_Ip, F_Motor_Y_Ip, F_Motor_Z_Ip, ContactState_Ip, Elongation_Ip, MotorSpecies_Ip)


	USE PARAMETERS, ONLY : DP, MaxNumParticles, NumMol, k, F_Motor_Detach1, F_Motor_Detach2


	INTEGER, INTENT(IN) :: MotorSpecies_Ip
	
	REAL(KIND = DP), INTENT(INOUT) :: F_Motor_X_Ip, F_Motor_Y_Ip, F_Motor_Z_Ip

	REAL(KIND = DP), DIMENSION(NumMol), INTENT(INOUT) :: ContactState_Ip

	REAL(KIND = DP), INTENT(IN) :: Elongation_Ip


	INTEGER :: II



	IF( ((MotorSpecies_Ip==1) .AND. (k*Elongation_Ip >= F_Motor_Detach1)) .OR. &
			((MotorSpecies_Ip==2) .AND. (k*Elongation_Ip >= F_Motor_Detach2)) )THEN

		F_Motor_X_Ip = 0.0_DP
		F_Motor_Y_Ip = 0.0_DP
		F_Motor_Z_Ip = 0.0_DP

		Scan_Contact_Filament: DO II=1, NumMol
			IF (ContactState_Ip(II) >= 1.0_DP) THEN
				ContactState_Ip(II) = -1.0_DP
			END IF
		END DO Scan_Contact_Filament

	END IF



END SUBROUTINE MotorForcedDetachment



!---Subroutine Calculating Force Acting on Bead---------------------------------

SUBROUTINE CalculateForceBead(ContactState, Force_Bead_X, Force_Bead_Y, Force_Bead_Z, F_Motor_X, F_Motor_Y, F_Motor_Z,NumParticles)


	USE PARAMETERS, ONLY : DP, MaxNumParticles, NumMol





	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(IN) :: ContactState

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(OUT) :: Force_Bead_X, Force_Bead_Y, Force_Bead_Z

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(IN) :: F_Motor_X, F_Motor_Y, F_Motor_Z

	INTEGER(KIND = Range15), INTENT(IN):: NumParticles

	REAL(KIND = DP) :: Intercept_Contact

	INTEGER(KIND = Range15) :: Ip

	INTEGER :: II, J_Contact



	Force_Bead_X = 0.0_DP
	Force_Bead_Y = 0.0_DP
	Force_Bead_Z = 0.0_DP



	DO Ip=1, NumParticles

		Scan_Contact_Filament: DO II=1, NumMol

			IF (ContactState(II,Ip) >= 1.0_DP) THEN

				J_Contact = INT(ContactState(II,Ip))
				Intercept_Contact = ContactState(II,Ip) - DBLE(J_Contact)

				IF (J_Contact >= NumBeads) THEN

					Force_Bead_X(NumBeads,II) = Force_Bead_X(NumBeads,II) + F_Motor_X(Ip)
					Force_Bead_Y(NumBeads,II) = Force_Bead_Y(NumBeads,II) + F_Motor_Y(Ip)
					Force_Bead_Z(NumBeads,II) = Force_Bead_Z(NumBeads,II) + F_Motor_Z(Ip)

				ELSE

					Force_Bead_X(J_Contact,II) = Force_Bead_X(J_Contact,II) + (1.0_DP - Intercept_Contact)*F_Motor_X(Ip)
					Force_Bead_Y(J_Contact,II) = Force_Bead_Y(J_Contact,II) + (1.0_DP - Intercept_Contact)*F_Motor_Y(Ip)
					Force_Bead_Z(J_Contact,II) = Force_Bead_Z(J_Contact,II) + (1.0_DP - Intercept_Contact)*F_Motor_Z(Ip)

					Force_Bead_X(J_Contact+1,II) = Force_Bead_X(J_Contact+1,II) + Intercept_Contact*F_Motor_X(Ip)
					Force_Bead_Y(J_Contact+1,II) = Force_Bead_Y(J_Contact+1,II) + Intercept_Contact*F_Motor_Y(Ip)
					Force_Bead_Z(J_Contact+1,II) = Force_Bead_Z(J_Contact+1,II) + Intercept_Contact*F_Motor_Z(Ip)

				END IF

			END IF

		END DO Scan_Contact_Filament

	END DO



END SUBROUTINE CalculateForceBead





!---Constraint--------------------------------------------------------

SUBROUTINE Constraint(XItempJ0, YItempJ0, ZItempJ0, XItempJ1, YItempJ1, ZItempJ1, XIJ0, YIJ0, ZIJ0, XIJ1, YIJ1, ZIJ1, StateConstraint)


	USE PARAMETERS, ONLY : DP, Tol, BondLength



	REAL(KIND = DP), INTENT(INOUT) :: XItempJ0, XItempJ1, YItempJ0, YItempJ1, ZItempJ0, ZItempJ1

	REAL(KIND = DP), INTENT(IN) :: XIJ0, XIJ1, YIJ0, YIJ1, ZIJ0, ZIJ1

	LOGICAL, INTENT(INOUT) :: StateConstraint

	REAL(KIND = DP) :: XI_tempAB, YI_tempAB, ZI_tempAB, XIAB, YIAB, ZIAB, BeadsDisSq, DiffSq, gAB, DX, DY, DZ


!	XItempJ0 = XI_temp(J), XItempJ1 = XI_temp(J+1), ...
!	XIJ0 = XI(J), XIJ1 = XI(J+1), ...


	XI_tempAB = XItempJ1 - XItempJ0
	YI_tempAB = YItempJ1 - YItempJ0
	ZI_tempAB = ZItempJ1 - ZItempJ0

	XIAB = XIJ1 - XIJ0
	YIAB = YIJ1 - YIJ0
	ZIAB = ZIJ1 - ZIJ0

	BeadsDisSq = XI_tempAB**2 + YI_tempAB**2 + ZI_tempAB**2
	DiffSq = BondLength**2 - BeadsDisSq

!	WRITE(*,*) BeadsDisSq, DiffSq

	IF (DABS(DiffSq) > 2.0_DP*Tol*BondLength**2) THEN

		StateConstraint = .FALSE.

		gAB = DiffSq / 4.0_DP / (XI_tempAB*XIAB + YI_tempAB*YIAB + ZI_tempAB*ZIAB)

		DX = gAB*XIAB
		DY = gAB*YIAB
		DZ = gAB*ZIAB

		XItempJ0 = XItempJ0 - DX
		YItempJ0 = YItempJ0 - DY
		ZItempJ0 = ZItempJ0 - DZ

		XItempJ1 = XItempJ1 + DX
		YItempJ1 = YItempJ1 + DY
		ZItempJ1 = ZItempJ1 + DZ

	END IF


END SUBROUTINE Constraint






!---Subroutine RenewMotorPopulation-----------------------------------

SUBROUTINE RenewMotorPopulation(X, Y, XP, YP, ZP, ContactState, CandidateList, AddedParticleNum, AreaCounter, EraseCounter, AreaOriginX, AreaOriginY, AreaOriginUx, AreaOriginUy, NumParticles, MotorSpecies)


	USE PARAMETERS, ONLY : Range15, DP, MaxNumParticles, NumMol, NumBeads, BondLength, HorizontalLength, VerticalLength, Motor_Density


	INTEGER, DIMENSION(MaxNumParticles), INTENT(INOUT) :: MotorSpecies		! 2Motors	1=Myosin, 2=Dead Motor (or "Zombie" in Dan's email)

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y

	REAL(KIND = DP), DIMENSION(MaxNumParticles), INTENT(INOUT) :: XP, YP, ZP

	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles), INTENT(INOUT) :: ContactState

	INTEGER, DIMENSION(NumMol*NumBeads, MaxNumParticles), INTENT(INOUT) :: CandidateList

	INTEGER(KIND = Range15), DIMENSION(MaxAreaNum), INTENT(INOUT) :: AddedParticleNum

	INTEGER, INTENT(INOUT)	:: AreaCounter, EraseCounter

	REAL(KIND = DP), DIMENSION(MaxAreaNum), INTENT(INOUT) :: AreaOriginX, AreaOriginY, AreaOriginUx, AreaOriginUy

	INTEGER(KIND = Range15), INTENT(INOUT) :: NumParticles



	INTEGER(KIND = Range15):: Ip, AddedParticleCounter

	INTEGER :: JJ, JJJ

	REAL(KIND = DP) :: XCM, YCM, XP_New, YP_New, UR1, UR2, UR3

	LOGICAL :: NearBoundary, InNewArea, OutOldArea





!---Erase Old Motor Region--------------------------------------------



	OutOldArea = .TRUE.

	DO JJJ=1, NumBeads
		IF ((DABS((X(JJJ,1)-AreaOriginX(EraseCounter+1))*AreaOriginUx(EraseCounter+1) + (Y(JJJ,1)-AreaOriginY(EraseCounter+1))*AreaOriginUy(EraseCounter+1)) <= 0.5_DP*HorizontalLength + 0.5_DP) .AND. &
		(DABS(-(X(JJJ,1)-AreaOriginX(EraseCounter+1))*AreaOriginUy(EraseCounter+1) + (Y(JJJ,1)-AreaOriginY(EraseCounter+1))*AreaOriginUx(EraseCounter+1)) <= 0.5_DP*VerticalLength + 0.5_DP)) THEN

			OutOldArea = .FALSE.

		END IF
	END DO


	IF (OutOldArea) THEN

		EraseCounter = EraseCounter + 1


!---Renew CandidateList-----------------------------------------------		!v25

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			CandidateList(:,Ip) = CandidateList(:,Ip + AddedParticleNum(EraseCounter))
		END DO



!---Renew ContactState------------------------------------------------		!v25

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			ContactState(:,Ip) = ContactState(:,Ip + AddedParticleNum(EraseCounter))
		END DO



!---Renew XP, YP, ZP--------------------------------------------------		!v25

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			XP(Ip) = XP(Ip + AddedParticleNum(EraseCounter))
		END DO

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			YP(Ip) = YP(Ip + AddedParticleNum(EraseCounter))
		END DO

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			ZP(Ip) = ZP(Ip + AddedParticleNum(EraseCounter))
		END DO


!---Renew MotorSpecies------------------------------------------------		!v25

		DO Ip=1, NumParticles - AddedParticleNum(EraseCounter)
			MotorSpecies(Ip) = MotorSpecies(Ip + AddedParticleNum(EraseCounter))
		END DO


		NumParticles = NumParticles - AddedParticleNum(EraseCounter)		!v25


	END IF





!---Checking if beads are near boundary-------------------------------

XCM=SUM(X(:,1))/DBLE(NumBeads)		!v11
YCM=SUM(Y(:,1))/DBLE(NumBeads)		!v11

!IF ((DABS((XCM-AreaOriginX(AreaCounter))*AreaOriginUx(AreaCounter) + (YCM-AreaOriginY(AreaCounter))*AreaOriginUy(AreaCounter)) >= 0.5_DP*HorizontalLength - 0.5_DP*DBLE(NumBeads-1)*BondLength - 1.0_DP) .OR. &			!v7
!	(DABS(-(XCM-AreaOriginX(AreaCounter))*AreaOriginUy(AreaCounter) + (YCM-AreaOriginY(AreaCounter))*AreaOriginUx(AreaCounter)) >= 0.5_DP*VerticalLength - 0.5_DP*DBLE(NumBeads-1)*BondLength - 1.0_DP)) THEN	!v7


NearBoundary = .FALSE.

DO JJ=1, NumBeads

	IF (DABS((X(JJ,1)-AreaOriginX(AreaCounter))*AreaOriginUx(AreaCounter) + (Y(JJ,1)-AreaOriginY(AreaCounter))*AreaOriginUy(AreaCounter)) >= 0.5_DP*HorizontalLength - 0.5_DP*DBLE(NumBeads-1)*BondLength - 1.0_DP) NearBoundary = .TRUE.

	IF (DABS(-(X(JJ,1)-AreaOriginX(AreaCounter))*AreaOriginUy(AreaCounter) + (Y(JJ,1)-AreaOriginY(AreaCounter))*AreaOriginUx(AreaCounter)) >= 0.5_DP*VerticalLength - 0.5_DP*DBLE(NumBeads-1)*BondLength - 1.0_DP) NearBoundary = .TRUE.

END DO


IF (NearBoundary == .TRUE.) THEN

!---Generate New Motor Region-----------------------------------------


!	WRITE(*,*) "CheckPoint1", TS

	AreaCounter = AreaCounter + 1


	AreaOriginX(AreaCounter) = XCM
	AreaOriginY(AreaCounter) = YCM
	AreaOriginUx(AreaCounter) = (X(1,1) - X(2,1))/DSQRT((X(1,1) - X(2,1))*(X(1,1) - X(2,1)) + (Y(1,1) - Y(2,1))*(Y(1,1) - Y(2,1)))
	AreaOriginUy(AreaCounter) = (Y(1,1) - Y(2,1))/DSQRT((X(1,1) - X(2,1))*(X(1,1) - X(2,1)) + (Y(1,1) - Y(2,1))*(Y(1,1) - Y(2,1)))


	AddedParticleCounter = 0
	DO Ip=1, INT(Motor_Density*HorizontalLength*VerticalLength, Range15)


		UR1 = grnd()
		UR2 = grnd()

!		XP(Ip) = HorizontalLength*(UR1 - 0.5_DP)
!		YP(Ip) = VerticalLength*(UR2 - 0.5_DP)

		XP_New = XCM + HorizontalLength * (UR1 - 0.5_DP) * AreaOriginUx(AreaCounter) - VerticalLength * (UR2 - 0.5_DP) * AreaOriginUy(AreaCounter)
		YP_New = YCM + HorizontalLength * (UR1 - 0.5_DP) * AreaOriginUy(AreaCounter) + VerticalLength * (UR2 - 0.5_DP) * AreaOriginUx(AreaCounter)

		InNewArea = .TRUE.
		IF (AreaCounter == 1) THEN																				!v11

			IF ((DABS((XP_New-AreaOriginX(AreaCounter))*AreaOriginUx(AreaCounter) + (YP_New-AreaOriginY(AreaCounter))*AreaOriginUy(AreaCounter)) <= 0.5_DP*HorizontalLength) .AND. &	!v11
			(DABS(-(XP_New-AreaOriginX(AreaCounter))*AreaOriginUy(AreaCounter) + (YP_New-AreaOriginY(AreaCounter))*AreaOriginUx(AreaCounter)) <= 0.5_DP*VerticalLength)) THEN		!v11

				InNewArea = .FALSE.																			!v11

			END IF																						!v11

		ELSE																							!v11
			DO I_Area = EraseCounter+1, AreaCounter-1																	!v11

				IF ((DABS((XP_New-AreaOriginX(I_Area))*AreaOriginUx(I_Area) + (YP_New-AreaOriginY(I_Area))*AreaOriginUy(I_Area)) <= 0.5_DP*HorizontalLength) .AND. &			!v11
				(DABS(-(XP_New-AreaOriginX(I_Area))*AreaOriginUy(I_Area) + (YP_New-AreaOriginY(I_Area))*AreaOriginUx(I_Area)) <= 0.5_DP*VerticalLength)) THEN				!v11

					InNewArea = .FALSE.																		!v11

				END IF																					!v11

			END DO																						!v11
		END IF																							!v11

		IF (InNewArea == .TRUE.) THEN

!			WRITE(*,*) "CheckPoint2"

			AddedParticleCounter = AddedParticleCounter + 1

			XP(NumParticles + AddedParticleCounter) = XP_New
			YP(NumParticles + AddedParticleCounter) = YP_New
			ZP(NumParticles + AddedParticleCounter) = 0.0_DP
			CandidateList(:,NumParticles + AddedParticleCounter) = 1	!v25

			UR1 = grnd()									! 2Motors
			IF (UR1 <= Species1Ratio) THEN							! 2Motors

				MotorSpecies(NumParticles + AddedParticleCounter) = 1			! 2Motors
				UR2 = grnd()								! 2Motors
				IF(UR2 <= 0.091)THEN							! 2Motors
					ContactState(:,NumParticles + AddedParticleCounter) = -1.0_DP	! 2Motors
				ELSE									! 2Motors
					ContactState(:,NumParticles + AddedParticleCounter) = 0.0_DP	! 2Motors
				END IF									! 2Motors
			ELSE										! 2Motors
				MotorSpecies(NumParticles + AddedParticleCounter) = 2			! 2Motors
				ContactState(:,NumParticles + AddedParticleCounter) = 0.0_DP		! 2Motors
			END IF										! 2Motors



!			WRITE(21,'(3F15.10)') XP_New, YP_New				!v29

		END IF



	END DO

	AddedParticleNum(AreaCounter) = AddedParticleCounter

	NumParticles = NumParticles + AddedParticleCounter



END IF



END SUBROUTINE RenewMotorPopulation





!---Motor State Conversion------------------------------------------------------

SUBROUTINE MotorStateConv(X, Y, Z, F_Motor_X_Ip, F_Motor_Y_Ip, F_Motor_Z_Ip, ContactState_Ip, Release_ADP_Ip)		! 2Motors

	USE PARAMETERS
	USE mtmod, ONLY : grnd

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	INTEGER, INTENT(INOUT) :: Release_ADP_Ip

!	INTEGER(KIND = Range15), INTENT(IN):: NumParticles

	REAL(KIND = DP), DIMENSION(NumMol), INTENT(INOUT) :: ContactState_Ip

	REAL(KIND = DP), DIMENSION(NumMol), INTENT(IN) :: F_Motor_X_Ip, F_Motor_Y_Ip, F_Motor_Z_Ip

	REAL(KIND = DP) :: UR, F_Motor_Tangent, k_d

	INTEGER(KIND = Range15) :: Ip

	INTEGER :: II, J_Contact



	DO II=1, NumMol

		J_Contact = INT(ContactState_Ip(II))

		UR = grnd()

		SELECT CASE(J_Contact)

		CASE (-1)

			IF(UR <= k_hp*dt)THEN
				ContactState_Ip(II) = 0.0_DP
			END IF

		CASE (0)

			IF(UR <= k_hm*dt)THEN
				ContactState_Ip(II) = -1.0_DP
			END IF

		CASE (1:)

			IF (Release_ADP_Ip == 0) THEN

				IF (J_Contact >= NumBeads) THEN
					F_Motor_Tangent = (F_Motor_X_Ip(II)*(X(J_Contact,II) - X(J_Contact-1,II))+ &
						F_Motor_Y_Ip(II)*(Y(J_Contact,II) - Y(J_Contact-1,II))+ &
						F_Motor_Z_Ip(II)*(Z(J_Contact,II) - Z(J_Contact-1,II)))/BondLength
				ELSE
					F_Motor_Tangent = (F_Motor_X_Ip(II)*(X(J_Contact+1,II) - X(J_Contact,II))+ &
						F_Motor_Y_Ip(II)*(Y(J_Contact+1,II) - Y(J_Contact,II))+ &
						F_Motor_Z_Ip(II)*(Z(J_Contact+1,II) - Z(J_Contact,II)))/BondLength
				END IF

				k_d = k_d0*DEXP(-F_Motor_Tangent*delta_x/KBT)

				IF(UR <= k_d*dt)THEN
					Release_ADP_Ip = 1
				END IF

			ELSE

				IF(UR <= k_t*ATP*dt)THEN 
					Release_ADP_Ip = 0
					ContactState_Ip(II) = -1.0_DP
				END IF

			END IF

		END SELECT

	END DO



END SUBROUTINE MotorStateConv





!---Bending Force Calculation-----------------------------------------

FUNCTION Force_Bending(Coordinate)


	USE PARAMETERS, ONLY : DP, NumMol, NumBeads, BondLength, EI


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol) :: Force_Bending

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: Coordinate

	INTEGER :: I, J

	REAL(KIND = DP) :: F


	Force_Bending = 0.0_DP

	DO I=1, NumMol
		DO J=2, NumBeads-1

			F = Coordinate(J+1,I) - 2.0_DP*Coordinate(J,I) + Coordinate(J-1,I)

			Force_Bending(J-1,I) = Force_Bending(J-1,I) + (-1.0_DP)*F
			Force_Bending(J,I) = Force_Bending(J,I) + 2.0_DP*F
			Force_Bending(J+1,I) = Force_Bending(J+1,I) + (-1.0_DP)*F

		END DO
	END DO

	Force_Bending = Force_Bending*EI/(BondLength**3)


END FUNCTION Force_Bending











!---Position of Bound Particle----------------------------------------

FUNCTION BoundParticlePosition(Ip, Coordinate, ContactState)


	USE PARAMETERS, ONLY : DP, NumMol, NumBeads, MaxNumParticles


	REAL(KIND = DP) :: BoundParticlePosition, Intercept_Contact

	INTEGER :: II, BonndCounter, J_Contact

	INTEGER, INTENT(IN) :: Ip

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: Coordinate

	REAL(KIND = DP), DIMENSION(NumMol, MaxNumParticles) :: ContactState



	BoundParticlePosition = 0.0_DP
	BonndCounter = 0
	DO II=1, NumMol

		IF (ContactState(II,Ip) > 1.0_DP) THEN

			BonndCounter = BonndCounter + 1

			J_Contact = INT(ContactState(II,Ip))
			Intercept_Contact = ContactState(II,Ip) - DBLE(J_Contact)

			BoundParticlePosition = BoundParticlePosition + &
				Coordinate(J_Contact,II) + Intercept_Contact * (Coordinate(J_Contact+1,II) - Coordinate(J_Contact,II))

		END IF

	END DO
	BoundParticlePosition = BoundParticlePosition/BonndCounter


END FUNCTION BoundParticlePosition


END PROGRAM BEADRODPOLYMER
