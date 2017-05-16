# ddm
Directional Dark Matter research project

Gabriel Penn and Will Sugg (Bristol MSci final year project students 2016/17)

## Running the simulation
This guide outlines the basic steps you will need to take to run the DDM simulation.

### Environment setup (once per login)




## Modifying the simulation

The sensitive volume is currently set to contain gaseous argon. The following parts of the simulation are specific to argon, and would have to be modified if the sensitive material is to be changed :

  - The generation of drift electrons, which uses the ionisation energy of argon
  - The equation used to calculate the electron drift velocity, obtained from [1]
  - The equation used to calculate the yield of scintillation photons generated per drift electron, obtained from [2]
