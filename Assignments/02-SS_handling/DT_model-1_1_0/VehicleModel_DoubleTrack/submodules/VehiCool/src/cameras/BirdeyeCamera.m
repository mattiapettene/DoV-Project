%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

classdef BirdeyeCamera < BaseCamera
    %% BirdeyeCamera
    % This class creates a simple birdeye camera object that follows a target
    % for the simulation. It is a subclass of BaseCamera, please refer to the
    % documentation of BaseCamera for more information.
    %
    % Properties
    % ----------
    %  - padding -> padding of the camera [m] with respect to the bounding box
    %               of the target. It is a 1 x 2 vector, where the first element
    %               is the padding in the x direction and the second element is
    %               the padding in the y direction. Default is [1, 1].
    % Methods
    % -------
    %  - BirdeyeCamera( target, varargin ) -> constructor.
    %  - create_state()                    -> create the camera state.
    %
    % Usage
    % -----
    %  - obj = BirdeyeCamera( target, varargin )
    %  - obj.create_state()
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        padding % padding of the camera [m] with respect to the bounding box of
                % the target

    end

    %% Methods
    methods

        % Constructor
        function obj = BirdeyeCamera( target, varargin )
            % This function creates the camera object for the simulation.
            %
            % Arguments
            % ----------
            %  - target          -> handle of the target. It has to be a
            %                       BaseObject object.
            %  - 'Projection'    -> camera projection. Options are 'perspective'
            %                       and 'orthographic'. Default is
            %                       'orthographic'.
            %  - 'Padding'       -> padding of the camera [m] with respect to
            %                       the bounding box of the target. It is a
            %                       1 x 2 vector, where the first element is the
            %                       padding in the x direction and the second
            %                       element is the padding in the y direction.
            %                       Default is [1, 1].
            %
            % Outputs
            % -------
            %  - obj -> camera object.
            %
            % Usage
            % -----
            %  - obj = BirdeyeCamera( target, varargin )
            %

            % Parse the inputs
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter( p, 'Projection', 'orthographic', @(x) ischar(x) && ismember(x, {'perspective', 'orthographic'} ) );
            addParameter( p, 'Padding', [1, 1], @(x) isnumeric(x) && isvector(x) && length(x) == 2 );
            parse( p, varargin{:} );

            % Call the superclass constructor
            obj@BaseCamera( target, ...
                            'Projection', p.Results.Projection, ...
                            p.Unmatched );

            % Set the camera padding
            obj.padding = p.Results.Padding;

        end

        % Update the camera view angle to frame the target
        function update_view_angle( obj )
            % This function computes the camera view angle to frame the target.
            %
            % Usage
            % -----
            %  - obj.update_view_angle( target )
            %

            % Get the target bounding box
            bbox = obj.target.bounding_box();

            % Compute the camera height relative to the target
            camera_height = 10 + bbox(2, 3);  % 3 times the height of the target

            % Compute the target dimensions on the image plane from the eight
            % corners of the bounding box
            target_length = abs( bbox(1, 1) - bbox(2, 1) );
            target_width  = abs( bbox(1, 2) - bbox(2, 2) );

            % Add the padding to the target dimensions
            target_length = target_length + obj.padding(1);
            target_width  = target_width  + obj.padding(2);

            % Compute the maximum dimension of the target on the image plane
            max_target_dim = max( target_length, target_width );

            % Compute the camera view angle to frame the target
            view_angle = 2 * rad2deg( atan2( max_target_dim / 2, ...
                                             camera_height ) );

            % Set the camera view angle
            obj.set_view_angle( view_angle );

        end

        % Get the camera current target
        function out = get_target_state( obj )
            % This function gets the camera current target.
            %
            % Outputs
            % -------
            %  - out -> camera current target.
            %
            % Usage
            % -----
            %  - target = obj.get_target_state()
            %

            % Get the target's bounding box
            bbox = obj.target.bounding_box();

            % extract the target's center
            out = [ (bbox(1, 1) + bbox(2, 1)) / 2, ...
                    (bbox(1, 2) + bbox(2, 2)) / 2, ...
                    (bbox(1, 3) + bbox(2, 3)) / 2 ];

        end

        % Create the camera state
        function out = create_state( obj )
            % This function creates the camera state.
            %
            % Usage
            % -----
            %  - obj.create_state()
            %

            % Set the camera view angle to frame the target
            obj.update_view_angle();

            % Get the target position
            target_position = obj.get_target_state();

            % Compute the camera position
            out = [ target_position(1), ...
                    target_position(2), ...
                    10 + 2 * target_position(3) ];

        end

    end

end
