# ddm
Directional Dark Matter research project

Gabriel Penn and Will Sugg (Bristol MSci final year project students 2016/17)

## Running the simulation
This guide outlines the basic steps you will need to take to run the DDM simulation once you have a local clone of this repository. (Please note that this guide assumes you are working on soolin.dice.priv, and may not work on other machines.)

### Environment setup (once per login)

A shell script is provided to set up the necessary Geant4 and ROOT environment using libraries from cvmfs.

From the ddm directory:
```
source init/init_geant4_plus_ROOT.sh
```

### Building the simulation (once, unless changes are made to the source code)

You may wish to use the shell script provided (ddm/DDM/DDM_build.sh) to build the simulation, but this assumes that your clone of the ddm repository is located in your home directory, and that you have created the expected build directory (~/geant4/DDM-build).

Otherwise, execute the following commands in the directory of your choice:

```
mkdir [desired build directory]
cd [desired build directory]
cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.2.p02/x86_64-slc6-gcc49-opt/lib64/ [path to ddm]/DDM
make -j4
```


## Modifying the simulation

The sensitive volume is currently set to contain gaseous argon. The following parts of the simulation are specific to argon, and would have to be modified if the sensitive material is to be changed :

  - The generation of drift electrons, which uses the ionisation energy of argon
  - The equation used to calculate the electron drift velocity, obtained from [1]
  - The equation used to calculate the yield of scintillation photons generated per drift electron, obtained from [2]

references
