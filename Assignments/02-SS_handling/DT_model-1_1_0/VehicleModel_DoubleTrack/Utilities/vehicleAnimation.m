%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

function vehicleAnimation( model_sim, vehicle_data, Ts )
% Temporarily undock the figures
set(0,'DefaultFigureWindowStyle','normal');

%------------------------------------------------------------------------------%
%% Create the scenario
%------------------------------------------------------------------------------%

% Extract the data
states        = model_sim.states;
time_sim      = states.x.Time;
state_vehicle = [states.x.Data, states.y.Data, 0 .* states.x.Data, ...
                 states.psi.Data, 0 .* states.psi.Data, 0 .* states.psi.Data]; % [x, y, z, yaw, roll, pitch]

% Construct the right and left margins using clothoids and the vehicle
% trajectory
CL_centre = ClothoidList();
CL_centre.build_G1( state_vehicle(:, 1), state_vehicle(:, 2) );

sample_dist  = linspace( 0, CL_centre.length(), floor( CL_centre.length() / 0.1 ) );
track_pad    = vehicle_data.vehicle.Wf + 1;
left_margin  = CL_centre.eval( sample_dist, track_pad )';
right_margin = CL_centre.eval( sample_dist, -track_pad )';

% Create a racetrack using the RaceTrack class. This class enables you to render
% any (we hope) flat racetrack. The only thing you need to do is to provide the
% left and right margins of the track. The margins are defined as a matrix
% where each row is a point of the margin. The first column is the x coordinate,
% the second column is the y coordinate and the third column is the z coordinate
% of the point.
track = RaceTrack( left_margin, right_margin );

% Create the vehicle using the STLObject class. This class enables
% you to render any STL object (which by default is a 3D reference frame). The
% only thing you need to do is to provide the state of the object as a N x 6
% matrix where each row is the state of the object at a given time. The first
% three columns are the position of the object, the last three columns are the
% orientation of the object.
vehicle_T = makehgtform( 'scale', 0.25, 'translate', [0, 0, 2.5] ); % initial transformation matrix
vehicle   = STLObject( state_vehicle, 'InitTrans', vehicle_T, ...
                       'STLPath', which( 'Renault_5_Rallye_Edition_Chassis.stl' ) );

% Create the camera using the FollowerCamera class. This class enables you to
% render a camera that follows a given object. The only thing you need to do is
% to provide the object to follow.
camera = FollowerCamera( vehicle );

% Create the scenario using the VehiCool class. This is the main class of the
% entire library.
scen = VehiCool();

% Add the objects to the scenario
scen.set_track( track );     % set the racetrack, we use “set” because every
                             % scenario has only one racetrack
scen.add_camera( camera );   % add the camera, which in theory could be more
                             % than one
scen.add_root_object( vehicle ); % add the vehicle as a root object,
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
scen.animate( time_sim(end), 'SampleTime', Ts, 'ShowProgress',  true, 'FrameRate', 24);

% Dock back the figures
set(0,'DefaultFigureWindowStyle','docked');

end
