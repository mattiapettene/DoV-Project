%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

classdef RDFTrack < BaseTrack
    %% RDFTrack
    % This class handles the creation of an RDFtrack object. It contains all the
    % information about the map of the RDFtrack. A RDFtrack is a track
    % which is defined by an RDF file.
    %
    % Disclaimer
    % ----------
    % For optimal results, please make sure the RDF file is structured as
    % follows:
    %
    %  [NODES]
    %  {id x_coord, y_coord, z_coord}
    %  0 0.0 0.0 0.0
    %  1 1.0 0.0 0.0
    %   ...
    %  [ELEMENTS]
    %  {node_1_id node_2_id node_3_id prop1 prop2 ...}
    %  0 1 2 0.0 0.0
    %  1 2 3 0.0 0.0
    %   ...
    %
    % The only requirements are that the sections [NODES] and [ELEMENTS] are
    % present, that the [NODES] section is structured as above, and that the
    % [ELEMENTS] section is structured as above up to the properties. The
    % properties are not used by this class, except for the colour of the
    % elements. The user can enable the use of the colour property by setting
    % the 'UseColour' to the name of the property when calling the constructor.
    % For example, if the colour property is called 'colour', then the user
    % should call the constructor as follows:
    %
    %  - obj = RDFTrack( rdf_filepath, 'UseColour', 'colour' ).
    %
    % We remind the user that the colour property must be a RGB triplet in the
    % range [0, 1]. We advise the user to not use blank spaces whenever they are
    % intended to be used as separators. For example, if the colour property is
    % defined as [0.0, 0.0, 0.0], then it should be written as [0.0,0.0,0.0].
    %
    % Properties
    % ----------
    %  - ax    -> axes handle.
    %  - rdf   -> the RDF file.
    %
    % Methods
    % -------
    %  - RDFTrack( rdf )  -> constructor.
    %  - parse_rdf( rdf ) -> parse the RDF file.
    %  - plot( ax )       -> plot the RDFtrack.
    %
    % Usage
    % -----
    %  - obj = RDFTrack( rdf )
    %  - obj.create_track( kerb_width )
    %  - obj.plot( ax )
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        ax  % axes handle
        rdf % the RDF file

    end

    %% Methods
    methods

        % Constructor
        function obj = RDFTrack( rdf_filepath, varargin )
            % RDFTrack constructor.
            %
            % Arguments
            % ---------
            %  - rdf_filepath -> path to the RDF file.
            %  - 'UseColour'  -> name of the property to use for the colour.
            %                    Default: ''. If set to '', then the colour
            %                    of the elements is set to the default colour,
            %                    which is [1.0, 1.0, 1.0].
            %
            % Outputs
            % -------
            %  - obj -> the RDFtrack object.
            %
            % Usage
            % -----
            %  - obj = RDFTrack( rdf_filepath, varargin )
            %

            % Parse the inputs
            p = inputParser;
            p.addRequired( 'rdf_filepath', @ischar );
            p.addParameter( 'UseColour', '', @ischar );
            p.parse( rdf_filepath, varargin{:} );

            % Parse the RDF file
            obj.parse_rdf( p.Results.rdf_filepath, p.Results.UseColour );

        end

        % Parse the RDF file
        function parse_rdf( obj, rdf_filepath, use_colour )
            % Parse the RDF file.
            %
            % Arguments
            % ---------
            %  - rdf_filepath -> path to the RDF file.
            %  - use_colour   -> name of the property to use for the colour.
            %
            % Outputs
            % -------
            %  - rdf -> the parsed RDF file.
            %
            % Usage
            % -----
            %  - obj.parse_rdf( rdf_filepath )
            %

            % Initialise the RDF file
            obj.rdf = struct( 'nodes', [], 'elements', [], 'colours', [] );

            % Open the RDF file
            fid = fopen( rdf_filepath, 'r' );

            % Read the file line by line
            line = fgetl( fid );
            while ischar( line )

                % Remove the trailing spaces
                line = strtrim( line );

                % If the line is a comment, skip it
                if line(1) == '#'
                    line = fgetl( fid );
                    continue;
                end

                % If the line is a section, read it
                if line(1) == '['

                    % Read the section name
                    section_name = line(2:end-1);

                    % Read the section
                    switch section_name

                        % Read the nodes
                        case 'NODES'
                            % Read the next line and skip it
                            line = fgetl( fid );

                            % Read the next line
                            line = fgetl( fid );

                            while ischar( line ) && line(1) ~= '['

                                % Remove the trailing spaces
                                line = strtrim( line );

                                % Split the line
                                line_split = strsplit( line, ' ' );

                                % Get the node coordinates
                                node_coords = str2double( line_split(2:4) );

                                % Store the node
                                obj.rdf.nodes(end + 1, :) = node_coords;

                                % Read the next line
                                line = fgetl( fid );

                            end

                        % Read the elements
                        case 'ELEMENTS'
                            % Read the next line
                            line = fgetl( fid );

                            % Check if the colour property is used
                            if ~isempty( use_colour ) && ischar( line )

                                % Remove the trailing spaces
                                line = strtrim( line );

                                % Split the line
                                line_split = strsplit( line, ' ' );

                                % Check if the colour property is present
                                compare_str = strcmp( line_split, ...
                                                      use_colour );
                                if ~any( compare_str )
                                    error( 'The colour property is not present.' );
                                end

                            end

                            % Read the next line
                            line = fgetl( fid );

                            while ischar( line ) && line(1) ~= '['

                                % Remove the trailing spaces
                                line = strtrim( line );

                                % Split the line
                                line_split = strsplit( line, ' ' );

                                % Get the element nodes
                                element_nodes = str2double( ...
                                                line_split(1:3) );

                                % Get the element colour if present
                                if ~isempty( use_colour ) && ...
                                    any( compare_str )
                                    element_colour = str2num( ...
                                                     line_split{compare_str} );
                                else
                                    element_colour = [0.0, 0.0, 0.0];
                                end

                                % Store the elements
                                obj.rdf.elements(end + 1, :) = element_nodes;

                                % Store the colours
                                obj.rdf.colours(end + 1, :) = element_colour;

                                % Read the next line
                                line = fgetl( fid );

                            end

                        % Skip the other sections
                        otherwise
                            line = fgetl( fid );
                            while ischar( line ) && line(1) ~= '['
                                line = fgetl( fid );
                            end

                    end

                end

            end

            % Close the file
            fclose( fid );

        end

        % Plot the RDFtrack
        function plot( obj, ax )
            % Plot the RDFtrack.
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

            % Plot the elements
            patch( 'Parent', obj.ax, ...
                   'Faces', obj.rdf.elements, ...
                   'Vertices', obj.rdf.nodes, ...
                   'FaceColor', 'flat', ...
                   'FaceVertexCData', obj.rdf.colours, ...
                   'EdgeColor', 'none', ...
                   'FaceAlpha', 1.0 );

        end

    end

end

