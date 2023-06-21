function effect_suspensions(vehicle_data,Ts)
    
    % Do a for cycle with different ratios for the suspensions:
    %
    epsilon_phi_original = (vehicle_data.front_suspension.Ks_f)/(vehicle_data.front_suspension.Ks_f + vehicle_data.rear_suspension.Ks_r);
    % ca = 0.4468
    epsilon_phi_vec = [0.3 0.4 epsilon_phi_original 0.5 0.6]
    Ks_f_vec = ((vehicle_data.rear_suspension.Ks_r).*epsilon_phi_vec)./(ones(5,1) +epsilon_phi_vec);
    
    % Same amplitude but different frequency
    p_sine_f_vec = [4000 f_n_num 2000 200 20];
    vo_sine_sim = zeros(round((Tsim/dt),TieBreaker="tozero") + 1,5);
    po_sin_vec = zeros(round((Tsim/dt),TieBreaker="tozero") + 1,5);
    for i = 1:5
        p_sine_freq = p_sine_f_vec(1,i);
        model = sim('Vehicle_Model_2Track');
        vo_sine_sim(:,i) = vo_sim.Data;
        po_sin_vec(:,i) = p_sim.Data;
    end



    model = sim('Vehicle_Model_2Track');

    % Put here the plots:

    
end
