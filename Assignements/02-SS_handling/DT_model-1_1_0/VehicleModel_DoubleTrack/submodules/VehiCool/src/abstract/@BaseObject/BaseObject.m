%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

classdef BaseObject < handle
    %% BaseObject
    % This class creates an object in 3D space with 6 Degrees of Freedom
    % (i.e., [x, y, z, roll, pitch, yaw]).
    %
    % This is an abstract class. It is meant to be inherited by other classes
    % that represent specific objects. The inherited classes must implement the
    % 'create_descriptor' method.
    %
    % Properties - abstract
    % ---------------------
    %  - bounding_box -> the bounding box of the object. It is a 2 x 3 matrix
    %                    composed of the two corners of the bounding box. The
    %                    first row contains the minimum values for each
    %                    coordinate, while the second row contains the maximum
    %                    values.
    %
    % Properties
    % ----------
    %  - ax                -> axes handle.
    %  - parent            -> handle of the parent object.
    %  - children          -> cell array of the handles to child objects.
    %  - state             -> the state of the object. It is a N x 3 matrix,
    %                         where N is the number of samples and each row
    %                         contains the object state [x, y, z, a1, a2, a3].
    %                         a1, a2, a3 are generic angles that can be used
    %                         to represent the object orientation. The rotation
    %                         sequence is defined by the 'RotSeq' parameter
    %                         defined in the constructor.
    %  - idx               -> the index of the current state.
    %  - deferred_update   -> multiple of the simulation sampling time to
    %                         update the object. This enables asynchronous
    %                         updates of objects. It must be a positive integer.
    %                         Default is 1.
    %  - rot_seq           -> the rotation sequence. It is a string of 3
    %                         characters, where each character can be 'x', 'y'
    %                         or 'z'. The rotation sequence is defined by the
    %                         'RotSeq' parameter defined in the constructor.
    %  - descriptor        -> the descriptor representing the object.
    %  - patch_transform   -> the transform object for the patch.
    %  - default_transform -> initial transform. It is a 4 x 4 transformation
    %                         matrix.
    %
    % Methods - abstract
    % ------------------
    %  - create_descriptor( obj )   -> create the object's descriptor.
    %  - update_bounding_box( obj ) -> update the bounding box.
    %
    % Methods
    % -------
    %  - BaseObject( state, varargin )              -> constructor.
    %  - set_state( state )                         -> set the object's state.
    %  - get_state()                                -> get the current state of
    %                                                  the object.
    %  - set_default_transform( default_transform ) -> set the object's default
    %                                                  transform.
    %  - get_transform()                            -> get the object's
    %                                                  transform.
    %  - define_tree( parent )                      -> define the object tree.
    %  - transform( T )                             -> transform the object.
    %  - skip( idx )                                -> skip to the given index.
    %  - plot( ax, varargin )                       -> plot the object.
    %  - plot_children( ax, varargin )              -> plot the children of the
    %                                                  object.
    %  - update( varargin )                         -> update the object pose.
    %  - update_children( varargin )                -> update the children of
    %                                                  the object.
    %
    % Usage
    % -----
    %  - obj = BaseObject( state, varargin )
    %  - obj.set_state( state )
    %  - out = obj.get_state()
    %  - obj.set_default_transform( default_transform )
    %  - T = obj.get_transform()
    %  - obj.define_tree( parent )
    %  - obj.transform( T )
    %  - obj.skip( idx )
    %  - obj.plot( ax, varargin )
    %  - obj.plot_children( ax, varargin )
    %  - obj.update( varargin )
    %  - obj.update_children( varargin )
    %

    %% Properties - abstract
    properties (Abstract, SetAccess = private, Hidden = true)

        bounding_box % the bounding box of the object

    end

    %% Properties
    properties (SetAccess = private, Hidden = true)

        ax                % axes handle
        parent            % handle of the parent object
        children          % list of the handles to child objects
        state             % the state of the object
        idx               % the index of the current state
        deferred_update   % multiple of the simulation sampling time to update
                          % the object
        rot_seq           % the rotation sequence
        descriptor        % the descriptor representing the object
        patch_transform   % the transform object for the patch
        default_transform % initial transform

    end

    %% Methods - abstract
    methods (Abstract)

        out = create_descriptor( obj ) % create the object's descriptor
        update_bounding_box( obj )     % update the bounding box

    end

    %% Methods
    methods

        % Constructor
        function obj = BaseObject( state, varargin )
            % BaseObject constructor
            %
            % Arguments
            % ---------
            %  - state            -> the state of the object. It is a N x 3
            %                        matrix, where N is the number of samples
            %                        and each row contains the object state
            %                        [x, y, z, a1, a2, a3]. a1, a2, a3 are
            %                        generic angles that can be used to
            %                        represent the object orientation. The
            %                        rotation sequence is defined by the
            %                        'RotSeq' parameter defined in the
            %                        constructor.
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
            %                        'RotSeq' parameter defined in the %
            %                        constructor. Default is 'zxy'.
            %  - 'Parent'         -> handle of the parent object. Default is [].
            %
            % Outputs
            % -------
            %  - obj  -> the object.
            %
            % Usage
            % -----
            %  - obj = BaseObject( state, varargin )
            %

            % Parse the inputs
            p = inputParser;
            addRequired( p, 'state', @isnumeric );
            addParameter( p, 'DeferredUpdate', 1, @(x) isnumeric(x) && x > 0 );
            addParameter( p, 'InitTrans', eye(4), @isnumeric );
            addParameter( p, 'RotSeq', 'zxy', @ischar );
            addParameter( p, 'Parent', [], @(x) or( isa( x, 'BaseObject' ), isempty( x ) ) );
            parse( p, state, varargin{:} );

            % Set the object properties
            obj.state             = p.Results.state;
            obj.idx               = 1;
            obj.deferred_update   = p.Results.DeferredUpdate;
            obj.default_transform = p.Results.InitTrans;
            rot_seq               = lower( p.Results.RotSeq );
            obj.rot_seq           = {[rot_seq(1), 'rotate'], ...
                                     [rot_seq(2), 'rotate'], ...
                                     [rot_seq(3), 'rotate']};

            % Set the parent
            obj.define_tree( p.Results.Parent );

        end

        % Set the object's state
        function set_state( obj, state )
            % Set the object's state
            %
            % Arguments
            % ---------
            %  - state -> the state of the object. It is a N x 3 matrix, where N
            %             is the number of samples and each row contains the
            %             object state [x, y, z, a1, a2, a3]. a1, a2, a3 are
            %             generic angles that can be used to represent the
            %             object orientation. The rotation sequence is defined
            %             by the 'RotSeq' parameter.
            %
            % Usage
            % -----
            %  - obj.set_state( state )
            %

            % Set the state
            obj.state = state;

        end

        % Get the current state
        out = get_state( obj );

        % Set the object's default transform
        function set_default_transform( obj, default_transform )
            % Set the object's default transformation matrix
            %
            % Arguments
            % ---------
            %  - default_transform -> the default transform. It is a 4 x 4
            %                         transformation matrix.
            %
            % Usage
            % -----
            %  - obj.set_default_transform( default_transform )
            %

            % Set the default transform
            obj.default_transform = default_transform;

        end

        % Get the object's transform
        function out = get_transform( obj )
            % Get the object's current transformation matrix
            %
            % Outputs
            % -------
            %  - out -> the transform. It is a 4 x 4 transformation matrix.
            %
            % Usage
            % -----
            %  - transform = obj.get_transform()
            %

            % Get the transform
            out = get( obj.patch_transform, 'Matrix' );

        end

        % Define the object tree
        define_tree( obj, parent );

        % Transform the object
        transform( obj, T );

        % Skip to a specific state
        skip( obj, idx );

        % Plot the object
        plot( obj, ax );

        % Update the object's position and orientation
        update( obj );

    end

end