#! /bin/bash

#---------------------------------------------------------------------
ulimit -s unlimited

./a.out
#
#aa=1
#while test $aa -le 1
#do
#
#	if test $aa -le 9
#	then
#	aa=00$aa
#	elif test $aa -le 99
#	then
#	aa=0$aa
#	fi
#
#	tar cfz Filament_A$aa.tar.gz Filament_A$aa*.vtk
#	tar cfz MabikiFilament_A$aa.tar.gz Filament_A$aa*0.vtk
#
#	tar cfz MTPlusEnd_A$aa.tar.gz MTPlusEnd_A$aa*.vtk
#	tar cfz MabikiMTPlusEnd_A$aa.tar.gz MTPlusEnd_A$aa*0.vtk
#
#	tar cfz Particle_A$aa.tar.gz Particle_A$aa*.vtk
#	tar cfz MabikiParticle_A$aa.tar.gz Particle_A$aa*0.vtk
#
#	tar cfz MotorSpecie1_A$aa.tar.gz MotorSpecie1_A$aa*.vtk
#	tar cfz MabikiMotorSpecie1_A$aa.tar.gz MotorSpecie1_A$aa*0.vtk
#
#	tar cfz MotorSpecie2_A$aa.tar.gz MotorSpecie2_A$aa*.vtk
#	tar cfz MabikiMotorSpecie2_A$aa.tar.gz MotorSpecie2_A$aa*0.vtk
#
#	tar cfz IntParticle_A$aa.tar.gz IntParticle_A$aa*.vtk
#	tar cfz MabikiIntParticleA$aa.tar.gz IntParticle_A$aa*0.vtk
#
#
#	tar cfz ContactStates_A$aa.tar.gz ContactStates_A$aa*
#
#	gzip KinesinTrajectories_A$aa.txt
#
#	gzip KinesinHeadTail_A$aa.txt
#
#
#	echo $aa
#	aa=$(expr $aa + 1)
#
#done
#rm Filament_A*.vtk
#rm MTPlusEnd_A*.vtk
#rm Particle_A*.vtk
#rm MotorSpecie1_A*.vtk
#rm MotorSpecie2_A*.vtk
#rm IntParticle_A*.vtk
#rm ContactStates_A*.txt
##rm AllKinesin_A*.vtk
#
##---------------------------------------------------------------------
#