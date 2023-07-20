%------------------------------------------------------------------------------%
%% Authors
%  - Sebastiano Taddei.
%  - Mattia Piazza.
%------------------------------------------------------------------------------%

function skip_to( objects, idx )
    % Skip the scenario to the given index.
    %
    % Arguments
    % ---------
    %  - objects -> cell array of objects to skip for.
    %  - idx     -> index to skip to. It must be an integer greater than or
    %               equal to 1.
    %
    % Usage
    % -----
    %  - VehiCool.skip( objects, idx )
    %

    % Parse the inputs
    p = inputParser;
    p.addRequired( 'objects', @iscell );
    p.addRequired( 'idx',     @(x) floor(x) == x && x >= 1 );
    p.parse( objects, idx );

    % Update the objects
    for i = 1:length( p.Results.objects )
        % Update the root object
        p.Results.objects{i}.skip( p.Results.idx );

        % Update its children
        if ~isempty( p.Results.objects{i}.children )
            VehiCool.skip_to( p.Results.objects{i}.children, p.Results.idx );
        end
    end

end