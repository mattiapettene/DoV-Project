%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

classdef FollowerCamera < BaseCamera
    %% FollowerCamera
    % This class creates a simple follower camera object for the simulation. It
    % is a subclass of BaseCamera, please refer to the documentation of
    % BaseCamera for more information.
    %
    % Properties
    % ----------
    %  - ax                -> axis handle.
    %  - state             -> camera state. It is a N x 3 matrix, where N is the
    %                         number of samples and each row contains the camera
    %                         position [x, y, z].
    %  - target            -> handle of the target. It has to be a BaseObject
    %                         object.
    %  - projection        -> camera projection. Options are 'perspective' and
    %                         'orthographic'. Default is 'perspective'.
    %  - view_angle        -> camera viewing angle [deg]. Default is 40 deg.
    %  - relative_position -> camera relative position w.r.t the target. Default
    %                         is [-10, 0, 6] (i.e., [x, y, z]).
    %
    % Methods
    % -------
    %  - FollowerCamera( target, varargin ) -> constructor.
    %  - set_target( target )               -> set camera target.
    %  - get_target_state()                       -> get camera current target.
    %  - create_state()                     -> create the camera state.
    %  - plot( ax, varargin )               -> set the camera.
    %  - update( varargin )                 -> update camera state and target.
    %
    % Usage
    % -----
    %  - obj = FollowerCamera( state, target, varargin )
    %  - obj.set_target( target )
    %  - out = obj.get_target_state()
    %  - out = obj.create_state()
    %  - obj.plot( ax, varargin )
    %  - obj.update( varargin )
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        relative_position % camera relative position w.r.t the target

    end

    %% Methods
    methods

        % Constructor
        function obj = FollowerCamera( target, varargin )
            % This function creates the camera object for the simulation.
            %
            % Arguments
            % ----------
            %  - target             -> handle of the target. It has to be a
            %                          BaseObject object.
            %  - 'Delay'            -> camera delay in following the target. It
            %                          is the ratio between the time step of the
            %                          simulation and the time constant of the
            %                          camera dynamics, which is assumed to be a
            %                          first order system. Default is 1 (i.e.,
            %                          no delay in following the target).
            %  - 'Projection'       -> camera projection. Options are
            %                          'perspective' and 'orthographic'. Default
            %                          is 'perspective'.
            %  - 'ViewAngle'        -> camera viewing angle [deg]. Default is
            %                          40.
            %  - 'RelativePosition' -> camera relative position w.r.t the
            %                          target. Default is [-10, 0, 6] (i.e.,
            %                          [x, y, z]).
            %
            %
            % Outputs
            % -------
            %  - obj -> camera object.
            %
            % Usage
            % -----
            %  - obj = FollowerCamera( target, varargin )
            %

            % Parse the inputs
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter( p, 'RelativePosition', [-10, 0, 6], @(x) isnumeric(x) && isvector(x) && length(x) == 3 );
            parse( p, varargin{:} );

            % Call the superclass constructor
            obj@BaseCamera( target, p.Unmatched );

            % Set the camera relative position
            obj.relative_position = [p.Results.RelativePosition, 1];

        end

        % Create the camera state
        function out = create_state( obj )
            % This function creates the camera state.
            %
            % Usage
            % -----
            %  - obj.create_state()
            %

            % Get the target current transformation matrix without the default
            % transformation
            target_trans = obj.target.get_transform() / ...
                           obj.target.default_transform;

            % Compute the camera position
            camera_state = target_trans * obj.relative_position';

            % Extract the camera position
            out = camera_state(1:3);

        end

    end

end
