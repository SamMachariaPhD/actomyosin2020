#! /bin/bash

ifort mt.f90 MotilityAssayActin2MotorsParameters_v6.f90 MotilityAssayConfinements_v2.f90 MotilityAssaySubstrateDeformation_v2.f90 MotilityAssayForceForceFunctions_v3.f90 MotilityAssayActin2MotorsMain_v9.f90

#-openmp
