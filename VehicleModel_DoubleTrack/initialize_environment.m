% -----------------------
%% Initialization
% -----------------------

close all; clear all; clc;
format long;

% Set LaTeX as default interpreter for axis labels, ticks and legends
set( 0,     'defaulttextinterpreter',          'latex' );
set( groot, 'defaultAxesTickLabelInterpreter', 'latex' );
set( groot, 'defaultLegendInterpreter',        'latex' );

set( 0, 'DefaultFigureWindowStyle', 'docked' );
set( 0, 'defaultAxesFontSize',      18 );
set( 0, 'DefaultLegendFontSize',    18 );

addpath( genpath( 'Lib_Vehicle' ) );
addpath( genpath( 'Lib_Tyre' ) );
addpath( genpath( 'Longit_Controller' ) );
addpath( genpath( 'Lateral_Controllers' ) );
addpath( genpath( 'Scenario' ) );
addpath( genpath( 'Utilities' ) );
addpath( 'Clothoids/matlab' );

% --------------------
%% Open Simulink model
% --------------------

% Check if the Simulink model is already opened. Otherwise open it
openModels = find_system( 'SearchDepth', 0 );
if ( isempty( find( strcmp( openModels, 'Vehicle_Model_2Track_contABS_EMPTY' ), 1 ) ) )
    load_system( 'Vehicle_Model_2Track_contABS_EMPTY.slx' );
    open_system( 'Vehicle_Model_2Track_contABS_EMPTY' );
end