# Dynamics of Vehicles - Course project
#### Team 6: Consalvi Natale, Pettene Mattia, Zumerle Matteo
<br>

## Assignment 1 - Tyre fitting
In [this](https://github.com/MatteoZumerle/DoV-Project/tree/main/Assignements/01-tyre_fitting) folder there are all the material related to the first assignment of the course, which aim is to fit the Pacejka MF96 tyre model coefficients on a dataset of measured forces for a F-SAE tyre. Using [Maple](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/01-tyre_fitting/TyreModel-MF96.mw) the coefficients needed to fit the data have been extracted from the Magic Formula equations and the MATLAB functions collected in [tyre_lib](https://github.com/MatteoZumerle/DoV-Project/tree/main/Assignements/01-tyre_fitting/tyre_lib) folder have been generated. The whole analysis have been carried out in [MATLAB](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/01-tyre_fitting/main_tyre_data_analysis.m), starting from two different [dataset](https://github.com/MatteoZumerle/DoV-Project/tree/main/Assignements/01-tyre_fitting/dataset). At the end, all the coefficents required by the MF96 tyre model have been obtained and collected [here](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/01-tyre_fitting/tyre_coeffs_team6.mat). In the final [report](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/01-tyre_fitting/Final_report_ass1_team6.pdf) there are all the plots and considerations related to this assignment.
<br> 

## Assignment 2 - Steady-state handling analysis
In [this](https://github.com/MatteoZumerle/DoV-Project/tree/main/Assignements/02-SS_handling) folder there are all the material related to the second assignment of the course, which aim is the analysis of the steady-state behaviour of a F-SAE formula type car, given its double track model and loading the tyre coefficients of the first assignment. Two different test have been carried out, modifying the [simulink](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/02-SS_handling/DT_model-1_1_0/VehicleModel_DoubleTrack/Vehicle_Model_2Track.slx) file:
- speed ramp test
- steer ramp test

The type of the test can be selected through the variable ```switch_test_type``` in the [main](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/02-SS_handling/DT_model-1_1_0/VehicleModel_DoubleTrack/mainVehicleModel_2Track.m) file. For each test, have been [computed](https://github.com/MatteoZumerle/DoV-Project/blob/main/Assignements/02-SS_handling/DT_model-1_1_0/VehicleModel_DoubleTrack/Utilities/dataAnalysis.m) the main parameters and indices, such as: lateral load transfer, normalized axle characteristics, handling diagram, undesteering gradient (fitted and theoretical), yaw gain, beta gain. Finally, the influence of suspensions stiffness, camber and toe angle have been [inspected](https://github.com/MatteoZumerle/DoV-Project/tree/main/Assignements/02-SS_handling/DT_model-1_1_0/VehicleModel_DoubleTrack/Utilities).

<br>

All the considerations and the analysis results for both the first and second assignments are collected in the final [presentation](https://github.com/MatteoZumerle/DoV-Project/blob/main/Final_presentation_team_6_DoV.pdf).
