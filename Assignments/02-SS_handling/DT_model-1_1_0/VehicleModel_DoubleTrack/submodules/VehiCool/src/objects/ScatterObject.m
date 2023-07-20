%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

classdef ScatterObject < BaseObject
    %% ScatterObject
    % This class creates a scatter object in 3D space with 6 Degrees of Freedom
    % (i.e., [x, y, z, roll, pitch, yaw]). It is a subclass of BaseObject,
    % please refer to the documentation of BaseObject for more information.
    %
    % Properties
    % ----------
    %  - data         -> the object's data. It is a N x 3 matrix, where N is the
    %                    number of samples and each row contains the [x, y, z]
    %                    coordinates of the points for the scatter plot.
    %  - colour       -> the colour for the scatter plot, 'White' by default. It
    %                    can be a triplet of RGB values in the range [0, 1] or a
    %                    matrix of N x 3 values where each row contains the RGB
    %                    values for the corresponding point. If instead an N x 1
    %                    array is provided, the colour will be inherited from
    %                    that array using MATLAB's default colour palette. To
    %                    change palette, set the parameter 'ColourPalette' in
    %                    the VehiCool constructor (since this setting
    %                    unfortunately applies to the entire animation).
    %  - size         -> the size for the markers in the scatter plot. Default
    %                    is 10.
    %  - bounding_box -> the bounding box of the object. It is a 2 x 3 matrix
    %                    composed of the two corners of the bounding box. The
    %                    first row contains the minimum values for each
    %                    coordinate, while the second row contains the maximum
    %                    values.
    %
    % Methods
    % -------
    %  - ScatterObject( state, varargin ) -> constructor.
    %  - create_descriptor()              -> create the object's descriptor.
    %
    % Usage
    % -----
    %  - obj = ScatterObject( state, varargin )
    %  - out = obj.create_descriptor()
    %

    %% Properties
    properties (SetAccess = private, Hidden = true)

        data         % the object's data
        colour       % the colour of the object
        size         % the size of the markers
        bounding_box % the bounding box of the object

    end

    %% Methods
    methods

        % Constructor
        function obj = ScatterObject( data, varargin )
            % ScatterObject constructor
            %
            % Arguments
            % ---------
            %  - data             -> the object's data. It is a N x 3 matrix,
            %                        where N is the number of samples and each
            %                        row contains the [x, y, z] coordinates of
            %                        the points for the scatter plot. If the
            %                        data is an N x 3 x M matrix, where M is the
            %                        number of sets of data, then the object
            %                        will be updated with each set of data based
            %                        on the 'DeferredUpdate' parameter defined
            %                        in the constructor.
            %  - 'State'          -> the state of the object. It is a N x 3
            %                        matrix, where N is the number of samples
            %                        and each row contains the object state
            %                        [x, y, z, a1, a2, a3]. a1, a2, a3 are
            %                        generic angles that can be used to
            %                        represent the object orientation. The
            %                        rotation sequence is defined by the
            %                        'RotSeq' parameter defined in the
            %                        constructor. Default is [0, 0, 0, 0, 0, 0]
            %                        for all samples.
            %  - 'DeferredUpdate' -> multiple of the simulation sampling time to
            %                        update the object. This enables
            %                        asynchronous updates of objects. It must be
            %                        a positive integer. Default is 1.
            %  - 'InitTrans'      -> the initial transformation for the STL. It
            %                        is a 4 x 4 transformation matrix. Default
            %                        is eye(4).
            %  - 'RotSeq'         -> rotation sequence for the transformations.
            %                        It is a string of 3 characters, where each
            %                        character can be 'x', 'y' or 'z'. The
            %                        rotation sequence is defined by the
            %                        'RotSeq' parameter defined in the
            %                        constructor. Default is 'zxy'.
            %  - 'Colour'         -> the colour for the scatter plot. It can be
            %                        a triplet of RGB values in the range [0, 1]
            %                        or a matrix of N x 3 values where each row
            %                        contains the RGB values for the
            %                        corresponding point. Default is 'White'. If
            %                        instead an N x 1 array is provided, the
            %                        colour will be inherited from that array
            %                        using MATLAB's default colour palette. To
            %                        change palette, set the parameter
            %                        'ColourPalette' in the VehiCool constructor
            %                        (since this setting unfortunately applies
            %                        to the entire animation). If the data is an
            %                        N x 3 x M matrix, where M is the number of
            %                        sets of data, then the object will be
            %                        updated with each set of data based on the
            %                        'DeferredUpdate' parameter defined in the
            %                        constructor. It could also be an N x 1 x M
            %                        matrix.
            %  - 'Size'           -> the size for the markers in the scatter
            %                        plot. Default is 10.
            %  - 'Parent'         -> handle of the parent object. Default is [].
            %
            % Outputs
            % -------
            %  - obj  -> the object.
            %
            % Usage
            % -----
            %  - obj = ScatterObject( state, varargin )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'data', @isnumeric );
            addParameter( p, 'State', zeros( 1, 6 ), @isnumeric );
            addParameter( p, 'DeferredUpdate', 1, @isnumeric );
            addParameter( p, 'InitTrans', eye( 4 ), @isnumeric );
            addParameter( p, 'RotSeq', 'zxy', @ischar );
            addParameter( p, 'Colour', rgb( 'White' ), @isnumeric );
            addParameter( p, 'Size', 10, @isnumeric );
            addParameter( p, 'Parent', [], @(x) or( isa( x, 'BaseObject' ), isempty( x ) ) );
            parse( p, data, varargin{:} );

            % Call the superclass constructor
            obj@BaseObject( p.Results.State,                            ...
                            'DeferredUpdate', p.Results.DeferredUpdate, ...
                            'InitTrans', p.Results.InitTrans,           ...
                            'RotSeq', p.Results.RotSeq,                 ...
                            'Parent', p.Results.Parent );

            % Set the object properties
            obj.data   = p.Results.data;
            obj.colour = p.Results.Colour;
            obj.size   = p.Results.Size;

            % Compute the bounding box
            obj.bounding_box = obj.compute_bounding_box();

        end

        % Get the current set of data from a dataset
        function out = get_set( obj, dataset )
            % Get the current set of data from a dataset
            %
            % Usage
            % -----
            %  - out = obj.get_set( dataset )
            %

            % Check if the object's data is a 3D matrix
            if ndims( dataset ) == 3

                out = dataset(:, :, obj.idx);

            else

                out = dataset;

            end

        end

        % Create the object's descriptor
        function out = create_descriptor( obj )
            % Create the object's descriptor
            %
            % Usage
            % -----
            %  - out = obj.create_descriptor()
            %

            % Get the object's data, colour, and size
            obj_data   = obj.get_set( obj.data );
            obj_colour = obj.get_set( obj.colour );
            obj_size   = obj.get_set( obj.size );

            % Create the patch object from the STL file
            out = scatter3( obj_data(:, 1), obj_data(:, 2), obj_data(:, 3), ...
                            obj_size, obj_colour );

        end

        % Compute the object's bounding box
        function out = compute_bounding_box( obj )
            % Compute the object's bounding box
            %
            % Usage
            % -----
            %  - out = obj.compute_bounding_box()
            %

            % Get the object's data
            obj_data = obj.get_set( obj.data );

            % Compute the bounding box
            out = [ min( obj_data, [], 1 ); ...
                    max( obj_data, [], 1 ) ];

        end

        % Update the bounding box
        function update_bounding_box( obj )
            % Update the bounding box
            %
            % Usage
            % -----
            %  - obj.update_bounding_box()
            %

            % Update the bounding box
            obj.bounding_box = obj.compute_bounding_box();

        end

        % Update the object
        function update( obj )
            % Plot the object
            %
            % Arguments
            % ---------
            %  - ax -> the axes handle.
            %
            % Usage
            % -----
            %  - obj.plot( ax )
            %

            % Get the object's data, colour, and size
            obj_data   = obj.get_set( obj.data );
            obj_colour = obj.get_set( obj.colour );
            obj_size   = obj.get_set( obj.size );

            % Update the object's descriptor only if it needs to be updated
            if ndims( obj.data ) == 3
                obj.descriptor.XData    = obj_data(:, 1);
                obj.descriptor.YData    = obj_data(:, 2);
                obj.descriptor.ZData    = obj_data(:, 3);
            end
            if ndims( obj.colour ) == 3
                obj.descriptor.CData    = obj_colour;
            end
            if ndims( obj.size ) == 3
                obj.descriptor.SizeData = obj_size;
            end

            % Call the superclass method at the end since it updates the index
            update@BaseObject( obj );

        end

    end

end