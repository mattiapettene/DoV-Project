%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

classdef FreeCamera < BaseCamera
    %% FreeCamera
    % This class creates a simple camera object able to follow a predefined path
    % for the simulation. It is a subclass of BaseCamera, please refer to the
    % documentation of BaseCamera for more information.
    %
    % Properties
    % ----------
    %  - path -> camera path to follow. It is an N x 3 matrix where each row
    %            represents a point in the path (i.e., [x, y, z] coordinates).
    %  - idx  -> current index of the path.
    %
    % Methods
    % -------
    %  - FreeCamera( path, target, ...
    %                varargin )        -> constructor.
    %  - create_state()                -> create the camera state.
    %
    % Usage
    % -----
    %  - obj = FreeCamera( state, target, varargin )
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        path % camera path to follow
        idx  % current index of the path

    end

    %% Methods
    methods

        % Constructor
        function obj = FreeCamera( path, target, varargin )
            % This function creates the camera object for the simulation.
            %
            % Arguments
            % ----------
            %  - path            -> camera path to follow. It is an N x 3 matrix
            %                       where each row represents a point in the
            %                       path (i.e., [x, y, z] coordinates).
            %  - target          -> handle of the target. It has to be a
            %                       BaseObject object.
            %  - 'Projection'    -> camera projection. Options are 'perspective'
            %                       and 'orthographic'. Default is
            %                       'perspective'.
            %  - 'ViewAngle'     -> camera viewing angle [deg]. Default is
            %                       40.
            %  - 'InitialIndex'  -> initial index of the path. Default is 1.
            %
            %
            % Outputs
            % -------
            %  - obj -> camera object.
            %
            % Usage
            % -----
            %  - obj = FreeCamera( target, varargin )
            %

            % Parse the inputs
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired( p, 'path', @(x) isnumeric(x) && ismatrix(x) && size(x, 2) == 3 );
            addParameter( p, 'InitialIndex', 1, @(x) isnumeric(x) && isscalar(x) );
            parse( p, path, varargin{:} );

            % Call the superclass constructor
            obj@BaseCamera( target, p.Unmatched );

            % Set the camera parameters
            obj.path = p.Results.path;
            obj.idx  = p.Results.InitialIndex;

        end

        % Create the camera state
        function out = create_state( obj )
            % This function creates the camera state. In particular, it returns
            % the current position of the camera.
            %
            % Usage
            % -----
            %  - obj.create_state()
            %

            % Return the current position of the camera
            out = obj.path( obj.idx, : );

            % Update the index
            obj.idx = obj.idx + 1;

        end

    end

end
