#! /bin/bash

./clean
cp codes/cg/ILU.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2ILU.csv

./clean
cp codes/cg/SSOR.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2SSOR.csv

./clean
cp codes/cg/trilinos/AMG.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILAMG.csv

./clean
cp codes/cg/trilinos/AMGMueLu.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILAMGMueLu.csv

./clean
cp codes/cg/trilinos/Chebyshev.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILChebyshev.csv

./clean
cp codes/cg/trilinos/ILU.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILILU.csv

./clean
cp codes/cg/trilinos/ILUT.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILILUT.csv

./clean
cp codes/cg/trilinos/SSOR.cc mark.cc
cmake ./
make release run >> ../../../performancedata/TRILSSOR.csv

./clean
cp codes/others/Bicgstab.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2BiCGStab.csv

./clean
cp codes/others/GMRES.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2GMRES.csv

./clean
cp codes/others/MINRES.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2MINRES.csv

./clean
cp codes/others/QMRS.cc mark.cc
cmake ./
make release run >> ../../../performancedata/D2QMRS.csv

./clean
cp codes/others/SDUMFPACK.cc
cmake ./
make release run >> ../../../performancedata/D2SDUMFPACK.csv


