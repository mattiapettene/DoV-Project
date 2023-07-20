# Double track vehicle model

## Authors

 - Francesco Biral       <francesco.biral@unitn.it>
 - Edaordo Pagot         <edoardo.pagot@unitn.it>
 - -Mattia Piccinini     <mattia.piccinini@unitn.it>
 - GastoneRosati Papini  <gastone.rosatipapini@unitn.it>
 - Sebastiano Taddei      <sebastiano.taddei@unitn.it>

## Description
Double track model with Simulink structure including ABS and ESP control blocks.
Simple sensors model block also included.
The model is parameterised and can be customised for different type of vehicles

## Revision history

### [1.1.0] - 2022-05-18
First version
#### Added
 - vehicle animation: an animation of the vehicle dynamics has been added at the end of the simulation using VehiCool library
 - FSAE data
#### Changed
 - adapted tyre formulas to match convention used for tyre forces used in fitting task. 
#### Removed

### [1.0.1] - 2022-05-18
First version
#### Added
 - added Fz saturation.
 - initial condition is now parametric in initial speed
 - changed the braking torque parameters to adapt to a GP2 model
#### Changed
 - corrected typos in the tyre forces formulas: wrong longitudinal slip variable used in Fx__rl and Fx__fl 

#### Removed

### [1.0.0] - 2022-05-18
First version
#### Added
#### Changed
#### Removed

