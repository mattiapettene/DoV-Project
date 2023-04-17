function previewPointParams = previewPointControllerParams()

    % ----------------------------------------------------------------
    %% Function purpose: define the preview point lateral controller parameters
    % ----------------------------------------------------------------
    
    % ke error to position
    previewPointParams.ke = 0.5; 
    % kp error to angle
    previewPointParams.kp = 1.2; 
    % lookhead
    previewPointParams.lookAhead = 10;

end

