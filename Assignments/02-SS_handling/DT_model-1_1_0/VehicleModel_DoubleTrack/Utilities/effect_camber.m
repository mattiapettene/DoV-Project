function effect_camber(vehicle_data,Ts,Tf)

    set(0,'defaulttextinterpreter','latex')
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
        
    % Do a for cycle with different camber angles:
    camber_vec = [ -1 -0.5 0 0.5 1 ];
    
    % Dovranno essere importati!
    time_sim = Tf;
    dt = Ts;

    alphaR_dt = zeros(round((time_sim/dt),TieBreaker="tozero") + 1,5);
    alphaF_dt = zeros(round((time_sim/dt),TieBreaker="tozero") + 1,5);
    Fyr_dt_norm = zeros(round((time_sim/dt),TieBreaker="tozero") + 1,5);
    Fyf_dt_norm = zeros(round((time_sim/dt),TieBreaker="tozero") + 1,5);
    delta_alpha_dt = zeros(round((time_sim/dt),TieBreaker="tozero") + 1,5);
    Ay = zeros(round((time_sim/dt),TieBreaker="tozero") + 1 -1 ,5);

    for i = 1:5
        vehicle_data.front_wheel.static_camber = camber_vec(1,i);
        assignin('base', "vehicle_data", vehicle_data);
        model_sim = sim('Vehicle_Model_2Track');
        fprintf('Starting simulation with $\\gamma_{front} = %.3f$\n', camber_vec(1,i));
        
        m = vehicle_data.vehicle.m;
        Lf = vehicle_data.vehicle.Lf;
        Lr = vehicle_data.vehicle.Lr;
        L = Lr + Lf;
        alpha_rr   = model_sim.states.alpha_rr.data;
        alpha_rl   = model_sim.states.alpha_rl.data;
        alpha_fr   = model_sim.states.alpha_fr.data;
        alpha_fl   = model_sim.states.alpha_fl.data;
        u          = model_sim.states.u.data;
        v          = model_sim.states.v.data;
        Omega      = model_sim.states.Omega.data;

        Fx_fr      = model_sim.extra_params.Fx_fr.data;
        Fx_fl      = model_sim.extra_params.Fx_fl.data;
        Fy_rr      = model_sim.extra_params.Fy_rr.data;
        Fy_rl      = model_sim.extra_params.Fy_rl.data;
        Fy_fr      = model_sim.extra_params.Fy_fr.data;
        Fy_fl      = model_sim.extra_params.Fy_fl.data;
        delta_fr   = model_sim.extra_params.delta_fr.data;
        delta_fl   = model_sim.extra_params.delta_fl.data;
        g = vehicle_data.vehicle.g;

        alphaR_dt(:,i) = (alpha_rr + alpha_rl)/2;
        alphaF_dt(:,i) = (alpha_fr + alpha_fl)/2;

        Fzr_0_static = m * (Lf/L) * g;
        Fzf_0_static = m * (Lr/L) * g;

        Fyr_dt_norm(:,i) = (Fy_rl + Fy_rr)./(Fzr_0_static);
        Fyf_dt_norm(:,i) = (sin(delta_fl).*Fx_fl + Fy_fl + sin(delta_fr).*Fx_fr + Fy_fr)./(Fzf_0_static);

        delta_alpha_dt(:,i) = alphaR_dt(:,i) - alphaF_dt(:,i);
        Ay(:,i) = diff(v(1:end))/dt + Omega(2:end).*u(2:end);
        %fprintf('Size of Ay: %4.2f - size of v = %4.2f - size of omega = %4.2f - size of u = %4.2f \n',size(Ay),size(diff(v(1:end))/dt),size(Omega),size(u));
        fprintf('End simulation #%d /5\n',i);
    end
    
    colorMap = colormap(lines(5));
    figure('Name','Normalized lateral forces for diff front camber angle','NumberTitle','off'), clf
    hold on
    for i = 1:5
        plot(alphaR_dt(:,i),Fyr_dt_norm(:,i),'LineWidth',2, 'DisplayName',("Yr(ay) with $\gamma_{front}$ = " +(-camber_vec(1,i))), 'Color', colorMap(i,:));
        plot(alphaF_dt(:,i),Fyf_dt_norm(:,i),'LineWidth',2, 'DisplayName',("Yf(ay) with $\gamma_{front}$ = " +(-camber_vec(1,i))), 'Color', colorMap(i,:));
        title({'Normalized lateral forces as function of $\gamma_{front}$', ' '})
        grid on
        ylabel('$Fyr/Fz0$, $Fyf/Fz0$ [-]')
        xlabel('$\alpha_{R}$, $\alpha_{F}$ [deg]')
        legend('Location','southeastoutside')
    end
    hold off

    figure('Name','Handling diagram for front camber angle','NumberTitle','off'), clf
    hold on
    for i = 1:5
        plot(Ay(:,i)./g, -delta_alpha_dt(2:end,i), 'LineWidth',2, 'DisplayName',("$-\Delta\alpha(ay)$ with $\gamma_{front}$ = " +(-camber_vec(1,i))), 'Color', colorMap(i,:));
        title({'Handling diagram as function of $\gamma_{front}$', ' '})
        grid on
        ylabel('$-\Delta\alpha$ [deg]')
        xlabel('$ay/g$ [-]')
        legend('Location','southeastoutside')
    end
    hold off
    
end
