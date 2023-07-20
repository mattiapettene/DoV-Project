%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

classdef RingTrack < BaseTrack
    %% RingTrack
    % This class handles the creation of a ringtrack object. It contains all the
    % information about the map of the ringtrack. A ringtrack is a track
    % composed of concentric circular roads with a central straight road.
    %
    % Properties
    % ----------
    %  - ax    -> axes handle.
    %  - rings -> N x 1 cell array of ring data structures. Each ring is
    %             composed of a centreline, a left margin, a right margin, a
    %             left kerb and a right kerb.
    %
    % Methods
    % -------
    %  - RingTrack( radius, num_rings, varargin )  -> constructor.
    %  - create_ring( radius, ring_width, ...
    %                 kerb_width, elevation )      -> create the centrelines,
    %                                                 the margins, and the
    %                                                 kerbs.
    %  - create_track( radius, num_rings,      ...
    %                  ring_width, kerb_width, ...
    %                  delta_elev )                -> create the centrelines,
    %                                                 the margins, and the
    %                                                 kerbs.
    %  - plot( ax )                                -> plot the ringtrack.
    %
    % Usage
    % -----
    %  - obj = RingTrack( radius, num_rings, varargin )
    %  - obj.create_ring( radius, ring_width, kerb_width, elevation )
    %  - obj.create_track( radius, num_rings, ring_width, kerb_width, ...
    %                      delta_elev )
    %  - obj.plot( ax )
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        ax    % axes handle
        rings % N x 1 cell array of ring data structures

    end

    %% Methods
    methods

        % Constructor
        function obj = RingTrack( radius, num_rings, varargin )
            % RingTrack constructor.
            %
            % Arguments
            % ---------
            %  - radius      -> radius of the innermost ring.
            %  - num_rings   -> number of rings.
            %  - 'RingWidth' -> width of the road. Default: 3.5.
            %  - 'KerbWidth' -> width of the kerbs. Default: 0.4.
            %  - 'DeltaElev  -> difference in elevation between the left and
            %                   right kerbs. Default: 0.
            %
            % Outputs
            % -------
            %  - obj -> the ringtrack object.
            %
            % Usage
            % -----
            %  - obj = RingTrack( radius, num_rings, varargin )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'radius',  @isnumeric );
            addRequired( p, 'num_rings', @isnumeric );
            addParameter( p, 'RingWidth', 3.5, @isnumeric );
            addParameter( p, 'KerbWidth', 0.2, @isnumeric );
            addParameter( p, 'DeltaElev', 0, @isnumeric );
            parse( p, radius, num_rings, varargin{:} );

            % Create the centrelines, the margins, and the kerbs
            obj.create_track( p.Results.radius, p.Results.num_rings, ...
                              p.Results.RingWidth, p.Results.KerbWidth, ...
                              p.Results.DeltaElev );

        end

        % Create a ring
        function create_ring( obj, radius, ring_width, kerb_width, elevation )
            % Create a ring.
            %
            % Arguments
            % ---------
            %  - radius     -> radius of the ring.
            %  - ring_width -> width of the ring.
            %  - kerb_width -> width of the kerbs.
            %  - elevation  -> elevation of the ring.
            %                  [left_kerb_elev, right_kerb_elev].
            %
            % Usage
            % -----
            %  - obj.create_ring( radius, ring_width, kerb_width, elevation )
            %

            % Compute the line between the left and right kerbs
            eval_elev = @(x) elevation(1) + (elevation(2) - elevation(1)) / ...
                             (radius + 2 * kerb_width) * x;

            % Angles
            angles = linspace( 0, 2 * pi, 1000 )';

            % Base coordinates
            cos_angles = cos( angles );
            sin_angles = sin( angles );

            % Create the centreline
            centreline = [radius * cos_angles, radius * sin_angles, ...
                          zeros( size( angles ) )];

            % Create the margins
            left_margin  = centreline + ring_width / 2 *       ...
                           [cos_angles, sin_angles, zeros( size( angles ) )] ...
                           + [zeros( size( angles ) ), ...
                              zeros( size( angles ) ), ...
                              eval_elev( kerb_width + radius ) * ...
                              ones( size( angles ) )];
            right_margin = centreline - ring_width / 2 *       ...
                           [cos_angles, sin_angles, zeros( size( angles ) )] ...
                           + [zeros( size( angles ) ), ...
                              zeros( size( angles ) ), ...
                              eval_elev( kerb_width ) * ones( size( angles ) )];

            % Create the kerbs
            left_kerb  = centreline + (ring_width / 2 + kerb_width) * ...
                         [cos_angles, sin_angles, zeros( size( angles ) )] ...
                         + [zeros( size( angles ) ), ...
                            zeros( size( angles ) ), ...
                            elevation(2) * ones( size( angles ) )];
            right_kerb = centreline - (ring_width / 2 + kerb_width) * ...
                         [cos_angles, sin_angles, zeros( size( angles ) )] ...
                         + [zeros( size( angles ) ), ...
                            zeros( size( angles ) ), ...
                            elevation(1) * ones( size( angles ) )];

            % Create the ring data structure
            ring.centreline   = centreline;
            ring.left_margin  = left_margin;
            ring.right_margin = right_margin;
            ring.left_kerb    = left_kerb;
            ring.right_kerb   = right_kerb;

            % Append the ring to the rings cell array
            obj.rings{end + 1} = ring;

        end

        % Create the centrelines, the margins, and the kerbs
        function create_track( obj, radius, num_rings, ring_width, ...
                               kerb_width, delta_elev )
            % Create the margins and kerbs.
            %
            % Arguments
            % ---------
            %  - radius     -> radius of the innermost ring.
            %  - num_rings  -> number of rings.
            %  - ring_width -> width of a ring.
            %  - kerb_width -> width of the kerb.
            %  - delta_elev -> difference in elevation between the left and
            %                  right kerbs.
            %
            % Usage
            % -----
            %  - obj.create_track( radius, num_rings, ring_width, ...
            %                      kerb_width, delta_elev )
            %

            % Create the rings
            for i = 1:num_rings
                radius_incr = (i - 1) * (ring_width + 2 * kerb_width);
                elevation   = [(i - 1) * delta_elev, i * delta_elev];
                obj.create_ring( radius + radius_incr, ring_width, ...
                                 kerb_width, elevation );
            end

        end

        % Plot the ringtrack
        function plot( obj, ax )
            % Plot the ringtrack.
            %
            % Arguments
            % ---------
            %  - ax -> axes handle.
            %
            % Usage
            % -----
            %  - obj.plot( ax )
            %

            % Set the axes handle
            obj.ax = ax;

            % Get the edges of the ringtrack
            x_min = min( [obj.rings{end}.left_kerb(:, 1); ...
                          obj.rings{end}.right_kerb(:, 1)] );
            x_max = max( [obj.rings{end}.left_kerb(:, 1); ...
                          obj.rings{end}.right_kerb(:, 1)] );
            y_min = min( [obj.rings{end}.left_kerb(:, 2); ...
                          obj.rings{end}.right_kerb(:, 2)] );
            y_max = max( [obj.rings{end}.left_kerb(:, 2); ...
                          obj.rings{end}.right_kerb(:, 2)] );

            % Plot the ground with a 100 m padding and a height of -1e-3 m
            ground_pad = 100;
            patch( ax, ...
                   [x_min - ground_pad, x_max + ground_pad,  ...
                    x_max + ground_pad, x_min - ground_pad], ...
                   [y_min - ground_pad, y_min - ground_pad,  ...
                    y_max + ground_pad, y_max + ground_pad], ...
                   zeros( 1, 4 ) - 1e-3,                     ...
                   rgb( 'ForestGreen' ),                     ...
                   'EdgeColor', 'none',                      ...
                   'FaceLighting', 'none',                   ...
                   'FaceAlpha', 1.0 );

            % Plot the track
            colour_flag = true;
            sample_skip = 9;
            for i = 1:sample_skip:(size( obj.rings{1}.left_margin, 1 ) - sample_skip)

                % Alternate the colour of the kerbs
                if colour_flag
                    left_colour  = rgb( 'White' );
                    right_colour = rgb( 'White' );
                    colour_flag  = false;
                else
                    left_colour  = rgb( 'DarkSlateGray' );
                    right_colour = rgb( 'DarkSlateGray' );
                    colour_flag  = true;
                end

                % For each ring
                for j = 1:length( obj.rings )

                    % Plot the left kerb
                    patch( ax,                                            ...
                           [obj.rings{j}.left_margin(i, 1),               ...
                            obj.rings{j}.left_margin(i + sample_skip, 1), ...
                            obj.rings{j}.left_kerb(i + sample_skip, 1),   ...
                            obj.rings{j}.left_kerb(i, 1)],                ...
                           [obj.rings{j}.left_margin(i, 2),               ...
                            obj.rings{j}.left_margin(i + sample_skip, 2), ...
                            obj.rings{j}.left_kerb(i + sample_skip, 2),   ...
                            obj.rings{j}.left_kerb(i, 2)],                ...
                           [obj.rings{j}.left_margin(i, 3),               ...
                            obj.rings{j}.left_margin(i + sample_skip, 3), ...
                            obj.rings{j}.left_kerb(i + sample_skip, 3),   ...
                            obj.rings{j}.left_kerb(i, 3)],                ...
                           left_colour,                                   ...
                           'EdgeColor', 'none',                           ...
                           'FaceLighting', 'none',                        ...
                           'FaceAlpha', 1.0 );

                    % Plot the track in between the kerbs
                    patch( ax,                                             ...
                           [obj.rings{j}.left_margin(i, 1),                ...
                            obj.rings{j}.left_margin(i + sample_skip, 1),  ...
                            obj.rings{j}.right_margin(i + sample_skip, 1), ...
                            obj.rings{j}.right_margin(i, 1)],              ...
                           [obj.rings{j}.left_margin(i, 2),                ...
                            obj.rings{j}.left_margin(i + sample_skip, 2),  ...
                            obj.rings{j}.right_margin(i + sample_skip, 2), ...
                            obj.rings{j}.right_margin(i, 2)],              ...
                           [obj.rings{j}.left_margin(i, 3),                ...
                            obj.rings{j}.left_margin(i + sample_skip, 3),  ...
                            obj.rings{j}.right_margin(i + sample_skip, 3), ...
                            obj.rings{j}.right_margin(i, 3)],              ...
                           rgb( 'DarkSlateGray' ),                         ...
                           'EdgeColor', 'none',                            ...
                           'FaceLighting', 'none',                         ...
                           'FaceAlpha', 1.0 );

                    % Plot the right kerb
                    patch( ax,                                             ...
                           [obj.rings{j}.right_margin(i, 1),               ...
                            obj.rings{j}.right_margin(i + sample_skip, 1), ...
                            obj.rings{j}.right_kerb(i + sample_skip, 1),   ...
                            obj.rings{j}.right_kerb(i, 1)],                ...
                           [obj.rings{j}.right_margin(i, 2),               ...
                            obj.rings{j}.right_margin(i + sample_skip, 2), ...
                            obj.rings{j}.right_kerb(i + sample_skip, 2),   ...
                            obj.rings{j}.right_kerb(i, 2)],                ...
                           [obj.rings{j}.right_margin(i, 3),               ...
                            obj.rings{j}.right_margin(i + sample_skip, 3), ...
                            obj.rings{j}.right_kerb(i + sample_skip, 3),   ...
                            obj.rings{j}.right_kerb(i, 3)],                ...
                           right_colour,                                   ...
                           'EdgeColor', 'none',                            ...
                           'FaceLighting', 'none',                         ...
                           'FaceAlpha', 1.0 );
                end

            end

        end

    end

end

