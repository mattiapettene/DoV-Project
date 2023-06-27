%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

classdef BaseCamera < handle
    %% BaseCamera
    % This class creates the camera object for the simulation.
    %
    % This is an abstract class. It is meant to be inherited by other classes
    % that represent specific cameras. The inherited classes must implement the
    % 'create_state' method.
    %
    % Properties
    % ----------
    %  - ax           -> axis handle.
    %  - state        -> camera state. It is a 1 x 3 matrix containing the
    %                    camera current position [x, y, z].
    %  - target       -> handle of the target. It has to be a BaseObject object.
    %  - target_state -> target state. It is a 1 x 3 matrix containing the
    %                    target current position [x, y, z].
    %  - delay        -> camera delay in following the target. It represents the
    %                    ratio between the time step of the simulation and the
    %                    time constant of the camera dynamics, which is assumed
    %                    to be a first order system. Default is 1 (i.e., no
    %                    delay in following the target).
    %  - projection   -> camera projection. Options are 'perspective' and
    %                    'orthographic'. Default is 'perspective'.
    %  - view_angle   -> camera viewing angle [deg]. Default is 40 deg.
    %
    % Methods - abstract
    % ------------------
    %  - out = create_state( obj ) -> create the state of the camera.
    %
    % Methods
    % -------
    %  - BaseCamera( target, varargin )           -> constructor.
    %  - set_target( target )                     -> set camera target.
    %  - out = get_target_state()                       -> get camera current target.
    %  - set_view_angle( view_angle )             -> set camera view angle.
    %  - out = first_order_dynamics( state, ...
    %                                next_state ) -> compute the next of a
    %                                                system following the
    %                                                dynamics of a first order
    %                                                system.
    %  - plot( ax )                               -> set the camera.
    %  - update()                                 -> update camera state and
    %                                                target.
    %
    % Usage
    % -----
    %  - obj = BaseCamera( state, target, varargin )
    %  - obj.set_target( target )
    %  - out = obj.get_target_state()
    %  - obj.set_view_angle( view_angle )
    %  - out = obj.first_order_dynamics( state, next_state )
    %  - obj.plot( ax )
    %  - obj.update()
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        ax           % axis handle
        state        % camera state
        target       % handle of the target
        target_state % target state
        delay        % camera delay in following the target
        projection   % camera projection
        view_angle   % camera viewing angle

    end

    %% Methods - abstract
    methods (Abstract)

        out = create_state( obj ) % create the state of the camera

    end

    %% Methods
    methods

        % Constructor
        function obj = BaseCamera( target, varargin )
            % This function creates the camera object for the simulation.
            %
            % Arguments
            % ----------
            %  - target       -> handle of the target. It has to be a BaseObject
            %                    object.
            %  - 'Delay'      -> camera delay in following the target. It is the
            %                    ratio between the time step of the simulation
            %                    and the time constant of the camera dynamics,
            %                    which is assumed to be a first order system.
            %                    Default is 1 (i.e., no delay in following the
            %                    target).
            %  - 'Projection' -> camera projection. Options are 'perspective'
            %                    and 'orthographic'. Default is 'perspective'.
            %  - 'ViewAngle'  -> camera viewing angle in degrees. The value must
            %                    be between 0 and 180. Default is 40 deg.
            %
            %
            % Outputs
            % -------
            %  - obj -> camera object.
            %
            % Usage
            % -----
            %  - obj = BaseCamera( target, varargin )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'target', @(x) isa( x, 'BaseObject' ) );
            addParameter( p, 'Delay', 1, @(x) isnumeric( x ) && x > 0 );
            addParameter( p, 'Projection', 'perspective',  @(x) ischar( x ) && any( strcmp( x, {'perspective', 'orthographic'} ) ) );
            addParameter( p, 'ViewAngle', 40, @(x) isnumeric( x ) && x > 0 && x < 180 );
            parse( p, target, varargin{:} );

            % Set the object properties
            obj.target       = p.Results.target;
            obj.delay        = p.Results.Delay;
            obj.projection   = p.Results.Projection;
            obj.view_angle   = p.Results.ViewAngle;

        end

        % Set camera target
        function set_target( obj, target )
            % This function sets the camera target.
            %
            % Arguments
            % ----------
            %  - target -> handle of the target. It has to be a BaseObject
            %              object.
            %
            % Usage
            % -----
            %  - obj.set_target( target )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'target', @(x) isa( x, 'BaseObject' ) );
            parse( p, target );

            % Set the target
            obj.target = p.Results.target;

        end

        % Get the camera current target
        function out = get_target_state( obj )
            % This function gets the camera current target.
            %
            % Outputs
            % -------
            %  - out -> state of camera current target.
            %
            % Usage
            % -----
            %  - target = obj.get_target_state()
            %

            out = obj.target.get_state();

        end

        % Set the view angle
        function set_view_angle( obj, view_angle )
            % This function sets the camera view angle.
            %
            % Arguments
            % ----------
            %  - view_angle -> camera view angle [deg].
            %
            % Usage
            % -----
            %  - obj.set_view_angle( view_angle )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'view_angle', @(x) isnumeric( x ) && x > 0 && x < 180 );
            parse( p, view_angle );

            obj.view_angle = p.Results.view_angle;

        end

        % Compute the next of a system following the dynamics of a first order
        % system
        out = first_order_dynamics( obj, state, next_state );

        % Set the camera
        plot( obj, ax );

        % Update camera state and target
        update( obj );

    end

end
