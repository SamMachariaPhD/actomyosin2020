
subroutine CalDragF(X1,Y1,X2,Y2,FX,FY,FZ)

USE PARAMETERS, ONLY : DP, ExtForceMag

real(kind = DP), intent(in) :: X1,Y1,X2,Y2 ! input: X1, Y1: The trailing bead; X2, Y2: The second trailing bead 
real(kind = DP), intent(out) :: FX,FY,FZ ! output resolved forces
real(kind = DP) :: Mag

Mag = sqrt((X1-X2)**2+(Y1-Y2)**2)
FX = (ExtForceMag/Mag)*(X1-X2)
FY = (ExtForceMag/Mag)*(Y1-Y2)
FZ = 0.0_DP

end subroutine CalDragF


