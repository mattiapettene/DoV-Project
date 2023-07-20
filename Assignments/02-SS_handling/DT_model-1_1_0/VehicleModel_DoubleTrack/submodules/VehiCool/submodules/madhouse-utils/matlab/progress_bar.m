%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%------------------------------------------------------------------------------%

function progress_bar( input, varargin )
    % Simple utility developed to easily incorporate a progress bar in any
    % MATLAB loop.
    %
    % Arguments
    % ---------
    %  - input       -> the input to be processed. If stop index is not
    %                   specified, it is assumed that input is the percentage
    %                   to be displayed. If stop index is specified, it is
    %                   assumed that input is the current index.
    %  - varargin{1} -> stop index.
    %  - varargin{2} -> start index. Default is 1.
    %  - varargin{3} -> step size of the loop. Default is 1.
    %  - 'Type'      -> type of progress bar:
    %                    - 'Percentage'    -> XXX%;
    %                    - 'Bar'           -> [===    ];
    %                    - 'BarPercentage' -> [===    ] XXX%.
    %                   Default is 'percentage'.
    %  - 'StepSize'  -> step size of the progress bar. Must be an integer
    %                   between 1 and 100. For optimal results, it should be a
    %                   divisor of 100 (otherwise the bar may not reach 100%).
    %                   Default is 2.
    %
    % Usage
    % -----
    %  - progress_bar( percentage )
    %  - progress_bar( index, stop_index )
    %  - progress_bar( index, stop_index, start_index )
    %  - progress_bar( index, stop_index, start_index, step_size )
    %  - progress_bar( ..., 'Type', type )
    %  - progress_bar( ..., 'StepSize', step_size )
    %

    % Process the inputs
    p = inputParser;
    addRequired(p, 'input', @isnumeric);
    addOptional(p, 'stop', [], @isnumeric);
    addOptional(p, 'start', 1, @isnumeric);
    addOptional(p, 'step', 1, @isnumeric);
    addParameter(p, 'Type', 'Percentage', @ischar);
    addParameter(p, 'StepSize', 2, @isnumeric);
    parse(p, input, varargin{:});

    % Check if the stop index is specified
    if isempty(p.Results.stop)

        % We assume that the input is the percentage to be displayed
        perc = p.Results.input;

    else

        % We assume that the input is the current index
        idx = p.Results.input;

        % Extract the loop parameters
        start = p.Results.start; % start index of the loop
        stop  = p.Results.stop;  % stop index of the loop
        step  = p.Results.step;  % step size of the loop

        % Compute the index of the last element of the loop. When defining a
        % range to loop over, the last element may not be included in the
        % range. For example, if we want to loop over the range [1, 10] with
        % a step size of 2, the last element of the range is 9, not 10.
        last = floor( (stop - start) / step ) * step + start;

        % Compute the percentage
        perc = ((idx - start) / (last - start)) * 100;

    end

    % Extract the bar properties
    type     = p.Results.Type;           % type of progress bar
    bar_step = p.Results.StepSize;       % step size of the progress bar
    bar_len  = floor( 100 / bar_step );  % length of the progress bar
    bar_adv  = floor( perc / bar_step ); % advancement of the bar

    % Check the type of progress bar
    switch type

        case 'Percentage'

            % If the index is not the first one, we need to delete the previous
            % percentage
            if idx ~= start
                fprintf( '\b\b\b\b\b%3d%%\n', bar_adv * bar_step );
            else
                fprintf( '%3d%%\n', bar_adv * bar_step );
            end

        case 'Bar'

            % If the index is not the first one, we need to delete the previous
            % percentage
            if idx ~= start
                fprintf([repmat( '\b', 1, 3 + bar_len ),      ...
                         '[',                                 ...
                         repmat( '=', 1, bar_adv ),           ...
                         repmat( ' ', 1, bar_len - bar_adv ), ...
                         ']\n']);
            else
                fprintf(['[',                               ...
                         repmat('=', 1, bar_adv),           ...
                         repmat(' ', 1, bar_len - bar_adv), ...
                         ']\n']);
            end

        case 'BarPercentage'

            % If the index is not the first one, we need to delete the previous
            % percentage
            if idx ~= start
                fprintf([repmat( '\b', 1, 3 + bar_len + 5 ),  ...
                         '[',                                 ...
                         repmat( '=', 1, bar_adv ),           ...
                         repmat( ' ', 1, bar_len - bar_adv ), ...
                         '] ', ...
                         '%3d%%\n'], bar_adv * bar_step);
            else
                fprintf(['[',                               ...
                         repmat('=', 1, bar_adv),           ...
                         repmat(' ', 1, bar_len - bar_adv), ...
                         '] ', ...
                         '%3d%%\n'], bar_adv * bar_step);
            end

        otherwise

            error('Invalid progress bar type.');

        % Note
        % ----
        % One could avoid to duplicate some of the code in the switch, but it
        % unfortuntely creates a strange jitter in the progress bar.

    end

end
