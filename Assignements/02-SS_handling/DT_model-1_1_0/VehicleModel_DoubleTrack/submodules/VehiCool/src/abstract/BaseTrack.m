%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

classdef (Abstract) BaseTrack < handle
    %% BaseTrack
    % Base class for all tracks in the library.
    %
    % This class is abstract and cannot be instantiated. It serves as a
    % template for all tracks in the library.
    %
    % Usage
    % -----
    % Whenever you want to create a new track, you should inherit from this
    % class. You should implement the following properties and methods:
    %  - properties (SetAccess = private, Hidden = true)
    %    - ax -> axes handle.
    %  - methods
    %    - plot( obj, ax ) -> plot the object.
    %

    %% Properties - abstract
    properties (Abstract, SetAccess = private, Hidden = true)

        ax % axes handle

    end

    %% Methods - abstract
    methods (Abstract)

        plot( obj, ax ) % plot the object

    end

end
