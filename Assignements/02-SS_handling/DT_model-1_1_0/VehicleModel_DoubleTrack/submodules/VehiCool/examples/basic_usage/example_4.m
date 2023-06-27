%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%
%% Prepare the workspace
%------------------------------------------------------------------------------%

% Clear the workspace
clear; close all; clc;

% Add the path to the VehiCool library
addpath( genpath( '../../' ) );

%------------------------------------------------------------------------------%
%% Load the data
%------------------------------------------------------------------------------%

load( '../data/adria.mat' );     % racetrack data
load( '../data/state_rf0.mat' ); % state of the main reference frame
load( '../data/time_sim.mat' );  % time vector of the simulation

%------------------------------------------------------------------------------%
%% Create the scenario
%------------------------------------------------------------------------------%

% Create a racetrack using the RaceTrack class. This class enables you to render
% any (we hope) flat racetrack. The only thing you need to do is to provide the
% left and right margins of the track. The margins are defined as a matrix
% where each row is a point of the margin. The first column is the x coordinate,
% the second column is the y coordinate and the third column is the z coordinate
% of the point.
track = RaceTrack( adria.left_margin, adria.right_margin );

% Create the main reference frame using the STLObject class. This class enables
% you to render any STL object (which by default is a 3D reference frame). The
% only thing you need to do is to provide the state of the object as a N x 6
% matrix where each row is the state of the object at a given time. The first
% three columns are the position of the object, the last three columns are the
% orientation of the object.
rf0_T = makehgtform( 'scale', 0.05 ); % initial transformation matrix
rf0   = STLObject( state_rf0, 'InitTrans', rf0_T );

% Create the camera using the FixedCamera class. This class enables you to
% render a camera that is fixed in a given position and that is pointing in the
% direction of the target. The only thing you need to do is to provide the
% reference frame of the camera and the position of the camera as a 1 x 3
% vector (by default the camera is placed in the origin of the world reference
% frame).
camera = FixedCamera( rf0, 'FixedPosition', [10, 10, 10] );

% Create the scenario using the VehiCool class. This is the main class of the
% entire library.
scen = VehiCool();

% Add the objects to the scenario
scen.set_track( track );     % set the racetrack, we use “set” because every
                             % scenario has only one racetrack
scen.add_camera( camera );   % add the camera, which in theory could be more
                             % than one
scen.add_root_object( rf0 ); % add the main reference frame as a root object,
                             % which in theory could be more than one

% Note
% The terminology “root object” is used to indicate an object that is not
% attached to any other object. In other words, it is an object that is not the
% child of any other object. In VehiCool, the root objects are the objects that
% are attached to the world frame. There is no need to add child objects to the
% scenario, because they are automatically added when you add a root object
% (doing so will result in an error).

%------------------------------------------------------------------------------%
%% Animate the scenario
%------------------------------------------------------------------------------%

% Animate the scenario using the animate method. This method takes as input the
% final time of the simulation. The scenario will be animated from time 0 to
% time time_sim(end).
scen.animate( time_sim(end) );

%------------------------------------------------------------------------------%
%% End of file
%------------------------------------------------------------------------------%

% User:
% I must admit, this is not the most exciting example. But creating a camera
% is not the most exciting thing to do, either. The point is that VehiCool
% is able to do that in an incredibly easy way, and that's what matters.
%
% VehiCool developers:
% Wait, weren't we supposed to be the ones to say that? Well, thank you user.
