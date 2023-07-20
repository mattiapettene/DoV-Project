%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

classdef VehiCool < handle
    %% VehiCool
    % This class handles the creation of a scenario (i.e., a track, a set of
    % vehicles, etc.) as well as the animation of the scenario.
    %
    % Properties
    % ----------
    %  - track       -> track of the scenario.
    %  - camera      -> camera view of the scenario.
    %  - objects     -> cell array of root objects in the scenario.
    %  - sample_time -> sample time of the data.
    %  - frame_rate  -> frame rate of the animation.
    %  - palette     -> colour palette for every plot in the scene.
    %
    % Methods
    % -------
    %  - VehiCool()                     -> constructor.
    %  - set_track( track )             -> set the track of the scenario.
    %  - add_camera( camera )           -> add a camera to the scenario.
    %  - add_root_object( root_object ) -> add a root object to the scenario.
    %  - render( varargin )             -> render the scenario.
    %  - advance( simulation_step )     -> advance the scenario by one step.
    %  - animate( tf, varargin )        -> animate the scenario.
    %
    % Methods - static
    % ----------------
    %  - plot_objects( objects, ax )       -> plot the objects of the scenario.
    %  - update_objects( objects, ...
    %                    simulation_step ) -> update the objects of the
    %                                         scenario.
    %  - skip_to( objects, idx )           -> skip the scenario to the given
    %                                         index.
    %
    % Usage
    % -----
    %  - obj = VehiCool()
    %  - obj.set_track( track )
    %  - obj.add_camera( camera )
    %  - obj.add_root_object( object )
    %  - obj.render( varargin )
    %  - obj.advance( simulation_step )
    %  - obj.animate( tf, varargin )
    %  - VehiCool.plot_objects( objects, ax )
    %  - VehiCool.update_objects( objects, simulation_step )
    %  - VehiCool.skip_to( objects, idx )
    %

    %% Properties - all private
    properties (SetAccess = private, Hidden = true)

        track       % track of the scenario
        camera      % camera view of the scenario
        objects     % cell array of root objects in the scenario
        sample_time % sample time of the data
        frame_rate  % frame rate of the animation
        palette     % colour palette for every plot in the scene

    end

    %% Methods
    methods

        % Constructor
        function obj = VehiCool( varargin )
            % VehiCool constructor.
            %
            % Arguments
            % ---------
            %  - 'ColourPalette' -> the colour palette to apply to all objects
            %                       in the scene (where applicable). It can be
            %                       any of the MATLAB's defaults (e.g.,
            %                       'parula') or an M x 3 matrix of RGB
            %                       triplets.
            %
            % Outputs
            % -------
            %  - obj -> the scenario object.
            %
            % Usage
            % -----
            %  - obj = VehiCool()
            %

            % Parse the arguments
            p = inputParser;
            p.addParameter( 'ColourPalette', 'parula', @(x) ischar(x) || ismatrix(x) && size(x, 2) == 3 );
            p.parse( varargin{:} );

            % Set the colour palette
            obj.palette = p.Results.ColourPalette;

        end

        % Set the track
        function set_track( obj, track )
            % Set the track of the scenario.
            %
            % Arguments
            % ---------
            %  - track -> track to add to the scenario. It must be a subclass
            %             of BaseTrack.
            %
            % Usage
            % -----
            %  - obj.set_track( track )
            %

            % Parse the inputs
            p = inputParser;
            p.addRequired( 'track', @(x) isa(x, 'BaseTrack') );
            p.parse( track );

            % Set the track
            obj.track = p.Results.track;

        end

        % Add a camera
        function add_camera( obj, camera )
            % Add a camera to the scenario.
            %
            % Arguments
            % ---------
            %  - camera -> camera to add to the scenario. It must be a subclass
            %              of BaseCamera.
            %
            % Usage
            % -----
            %  - obj.add_camera( camera )
            %

            % Parse the inputs
            p = inputParser;
            p.addRequired( 'camera', @(x) isa(x, 'BaseCamera') );
            p.parse( camera );

            % Set the camera
            obj.camera = p.Results.camera;

        end

        % Add an object
        function add_root_object( obj, object )
            % Add a root object to the scenario.
            %
            % Arguments
            % ---------
            %  - object -> object to add to the scenario. It must be a subclass
            %              of BaseObject.
            %
            % Usage
            % -----
            %  - obj.add_root_object( object )
            %

            % Parse the inputs
            p = inputParser;
            p.addRequired( 'object', @(x) isa(x, 'BaseObject') );
            p.parse( object );

            % Add the object
            obj.objects{end + 1} = p.Results.object;

        end

        % Render the scenario
        [fig, ax] = render( obj, varargin );

        % Advance the scenario
        advance( obj, simulation_step );

        % Animate the scenario
        animate( obj, tf, varargin );

    end

    %% Methods - static
    methods (Static)

        % Plot the objects
        plot_objects( objects, ax );

        % Update the objects
        update_objects( objects, simulation_step );

        % Skip the scenario to the given index
        skip_to( objects, idx );

    end

end
