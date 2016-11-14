cd ~/geant4/DDM-build 
rm -rf ./*
cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.2.p02/x86_64-slc6-gcc49-opt/lib64/ ~/ddm/DDM
make -j4
