%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function [fig, ax] = render( obj, varargin )
    % Render the scenario.
    %
    % Arguments
    % ---------
    %  - 'FigSize'     -> 1 x 2 matrix containing the size of the figure (i.e.,
    %                     [width, height]). Default is [960, 540].
    %  - 'ShowFigure' -> flag to show the figure or not. Default is
    %                    'on'. Options are 'on' and 'off'.
    %
    % Outputs
    % -------
    %  - fig -> figure handle.
    %  - ax  -> axes handle.
    %
    % Usage
    % -----
    %  - obj.render( varargin )
    %

    % Parse the inputs
    p = inputParser;
    addParameter( p, 'FigSize', [960, 540], @(x) isnumeric( x ) && isequal( size( x ), [1, 2] ) );
    addParameter( p, 'ShowFigure', 'on', @(x) ischar( x ) && ismember( x, {'on', 'off'} ) );
    parse( p, varargin{:} );

    % Create the figure and axes
    fig = figure( 'Name', 'VehiCool', 'NumberTitle', 'off', ...
                  'Position', [0, 0, p.Results.FigSize],    ...
                  'Visible', p.Results.ShowFigure );
    ax  = axes( 'Parent', fig, 'Visible', 'off' );

    % Set the figure properties
    hold( ax, 'on' );
    pbaspect( ax, [1, 1. 1] );   % axis aspect ratio
    daspect( ax, [1, 1, 1] );    % axis ticks aspect ratio
    colormap( ax, obj.palette ); % set the colormap

    % Set the scene lighting
    light( ax, 'Style', 'infinite' );
    lighting( ax, 'gouraud' );

    % Plot the track
    obj.track.plot( ax );

    % Plot the objects
    VehiCool.plot_objects( obj.objects, ax );

    % Plot the camera
    obj.camera.plot( ax );

end