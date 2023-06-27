%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function animate( obj, tf, varargin )
    % Animate the scenario.
    %
    % Arguments
    % ---------
    %  - tf            -> total animation time.
    %  - 'FrameRate'   -> frame rate of the animation. Default is 30.
    %  - 'SampleTime'  -> sample time of the data. Default is 0.001.
    %  - 'FigSize'     -> 1 x 2 matrix containing the size of the figure (i.e.,
    %                     [width, height]). Default is [960, 540].
    %  - 'ShowProgress'-> flag to show the progress bar or not. Default
    %                     is false.
    %  - 'ShowFigure'  -> flag to show the figure or not. Default is
    %                     'on'. Options are 'on' and 'off'.
    %  - 'SaveVideo'   -> flag to save the video or not. Default is
    %                     false.
    %  - 'FileName'    -> name of the video file. Default is 'VehiCool'.
    %  - 'FileFormat'  -> format of the video file. Default is 'MPEG-4'. Any of
    %                     the MATLAB's default formats can be used.
    %  - 'FileQuality' -> quality of the video file expressed as a value in the
    %                     range [0, 100]. Default is 100.
    %
    % Usage
    % -----
    %  - obj.animate( tf, varargin )
    %

    % Parse the inputs
    p = inputParser;
    addRequired( p, 'tf', @isnumeric );
    addParameter( p, 'FrameRate', 30, @(x) isnumeric( x ) && x > 0 );
    addParameter( p, 'SampleTime', 0.001, @(x) isnumeric( x ) && x > 0 );
    addParameter( p, 'FigSize', [960, 540], @(x) isnumeric( x ) && isequal( size( x ), [1, 2] ) );
    addParameter( p, 'ShowProgress', false, @islogical );
    addParameter( p, 'ShowFigure', 'on', @(x) ischar( x ) && ismember( x, {'on', 'off'} ) );
    addParameter( p, 'SaveVideo', false, @islogical );
    addParameter( p, 'FileName', 'VehiCool', @ischar );
    addParameter( p, 'FileFormat', 'MPEG-4', @ischar );
    addParameter( p, 'FileQuality', 100, @(x) isnumeric( x ) && x >= 0 && x <= 100 );
    parse( p, tf, varargin{:} );

    % Check that SampleTime <= 1 / FrameRate
    if p.Results.SampleTime > 1 / p.Results.FrameRate
        error( 'ERROR: SampleTime > 1 / FrameRate.' );
    end

    % Render the scenario
    [fig, ax] = obj.render( 'FigSize', p.Results.FigSize, ...
                            'ShowFigure', p.Results.ShowFigure );

    % Create the video file if needed
    if p.Results.SaveVideo
        vidfile           = VideoWriter( p.Results.FileName, ...
                                         p.Results.FileFormat );
        vidfile.FrameRate = p.Results.FrameRate;
        vidfile.Quality   = p.Results.FileQuality;
        open( vidfile );
    end

    % Round the inverse of the frame rate to the same order of the
    % sample time
    frame_time = round( 1 / p.Results.FrameRate, ...
                        -floor( log10( p.Results.SampleTime ) ) );

    % Skip to the current index
    idx = 1;
    VehiCool.skip_to( obj.objects, idx );

    % Simulate the scenario
    for t = 0:p.Results.SampleTime:tf

        % Advance the scenario
        s_adv = tic;
        obj.advance( idx );
        e_adv = toc( s_adv );

        % Simulate the scenario only at the specified frame rate
        if mod( t, frame_time ) == 0

            % Render the scenario
            s_ren = tic;
            drawnow nocallbacks; % nocallbacks disables the callbacks of the
                                 % figure, which would slow down the
                                 % animation if the user is interacting with it
            e_ren = toc( s_ren ) + e_adv;

            % Check if the user wants to save the animation
            if p.Results.SaveVideo
                A = getframe( fig );
                writeVideo( vidfile, A );
            else
                % Try to pause for the remaining time
                if e_ren < 1 / p.Results.FrameRate
                    pause( 1 / p.Results.FrameRate - e_ren );
                end
            end
        end

        % Advance the progress bar
        if p.Results.ShowProgress
            progress_bar( t, tf, 0, p.Results.SampleTime, ...
                          'Type', 'BarPercentage' );
        end

        % Increase the index
        idx = idx + 1;

    end

    % Stop holding onto the figure
    hold( ax, 'off' );

    % Close the video if needed
    if p.Results.SaveVideo
        close( vidfile );
    end

end