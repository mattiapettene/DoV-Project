

% Vehicle data

new_vehicle.Lf       =  1.0   ; % [m]     distance of vehicle c.g from front axle
new_vehicle.Lr       =  0.7   ; % [m]     distance of vehicle c.g rear front axle
new_vehicle.h_G      =  0.295 ; % [m]     vehicle centre of gravity height
new_vehicle.m_s      =  225 ; % [kg]     sprung mass 
new_vehicle.Jx_s     =  65  ; % [kgm^2]  sprung mass roll inertia
new_vehicle.Jy_s     =  90  ; % [kgm^2]  sprung mass pitch inertia
new_vehicle.Jz_s     =  130 ; % [kgm^2]  sprung mass yaw inertia
new_vehicle.Jxz_s    =  20  ; % [kgm^2]   sprung mass ox-oz product of inertia
new_vehicle.m_uf     =  14     ; % [kg]   front single unsprung mass
new_vehicle.Jxu_f    =  0.2    ; % [kgm^2] front unsprung mass cambering moment of inertia
new_vehicle.Jyu_f    =  0.2    ; % [kgm^2] front unsprung mass pitching moment of inertia 
new_vehicle.Jzu_f    =  0.6    ; % [kgm^2] front unsprung mass yawing moment of inertia 
new_vehicle.m_u_r    =  18.5   ; % [kg]    rear axle unsprung mass
new_vehicle.Jxu_r    =  0.3    ; % [kgm^2] rear unsprung mass cambering moment of inertia 
new_vehicle.Jyu_r    =  0.3    ; % [kgm^2] rear unsprung mass pitching moment of inertia 
new_vehicle.Jzu_r    =  0.5    ; % [kgm^2] rear unsprung mass yawing moment of inertia 
new_vehicle.Jw_f     =  0.25   ; % [kgm^2] front wheel spin moment of inertia
new_vehicle.Jw_f     =  0.25   ; % [kgm^2] rear wheel spin moment of inertia 
new_vehicle.W_f      =  1.260  ; %  [m]      front track width
new_vehicle.W_r      =  1.210  ; %  [m]     rear track width
new_vehicle.delta_f0 =  0      ; % [deg]  static front toe angle
new_vehicle.delta_r0 =  0      ; % [deg] static rear toe angl
new_vehicle.gamma_f  = 0      ; % [deg] static front camber angle
new_vehicle.gamma_r  = 0      ; % [deg] static rear camber angle

new_vehicle.h_rf     = 0.025 ; %[m] front roll centre height
new_vehicle.h_rr     = 0.045 ; %[m] rear roll centre height
new_vehicle.k_phi_f  = 205;  %[Nm/deg] front roll stiffness
new_vehicle.k_phi_r  = 175;  %[Nm/deg]rear roll stiffness
new_vehicle.c_phi_f  = 225;  %[Nm s/deg]front roll damping
new_vehicle.c_phi_r  = 270;  %[Nm s/deg]rear roll damping
new_vehicle.k_phi_c  = 380; %[Nm/deg]chassis torsional stiffness


new_vehicle.tau_H          = 1/12; %[-] steering ratio
new_vehicle.Jr_f           = 30; % [kgm^2] front half roll inertia
new_vehicle.Jr_r           = 30; % [kgm^2] rear half roll inertia
new_vehicle.J_z            = 165; % [kgm^2]total yaw inertia
new_vehicle.epsilon_roll_f = 0.05; %[-]front roll steer coefficient 
new_vehicle.epsilon_roll_f = 0.1; %[-] rear roll steer coefficient


%  Front left suspension hard points, rest configuration [X,Y,Z]=[mm,mm,mm]

new_vehicle.P_A_f =[ -280.0  , 249.7   , 145.0 ]; % lower wishbone rear pivot
new_vehicle.P_B_f =[   80.0  , 249.7   , 145.0 ]; % lower wishbone front pivot
new_vehicle.P_C_f =[    5.0  , 585.6   , 140.0 ]; % lower wishbone outer ball joint
new_vehicle.P_D_f =[    5.0  , 550.5   , 165.0 ]; % push rod wishbone end
new_vehicle.P_E_f =[ -280.0  , 259.7   , 345.0 ]; % upper wishbone rear pivot
new_vehicle.P_F_f =[    80.0  , 259.8  , 345.0 ]; % upper wishbone front pivot
new_vehicle.P_H_f =[    -5.0  , 560.9  , 380.3 ]; % upper wishbone outer ball joint
new_vehicle.P_L_f =[   -61.0  , 248.8  , 315.0 ]; % inner track rod ball joint
new_vehicle.P_M_f =[   -61.0  , 540.0  , 344.0 ]; % outer track rod ball joint
new_vehicle.P_P_f =[   100.0  , 210.0  , 432.7 ]; % rocker axis 1st point
new_vehicle.P_Q_f =[   100.0  , 310.0  , 595.8 ]; % rocker axis 2nd point
new_vehicle.P_R_f =[    158.0  , 250.0 , 475.0 ]; % push rod rocker end
new_vehicle.P_S_f =[    100.0  , 160.0 , 530.0 ]; % damper to rocker point
new_vehicle.P_T_f =[   -170.0  , 160.0 , 530.0 ]; % damper to body point

% Rear left suspension hard points, rest configuration. [X,Y,Z]=[mm,mm,mm]
new_vehicle.P_A_r = [-280.0   , 249.7, 145.0]; % lower wishbone rear pivot
new_vehicle.P_B_r = [  80.0  , 249.7,  145.0]; % lower wishbone front pivot
new_vehicle.P_C_r = [   5.0  , 585.6,  140.0]; % lower wishbone outer ball joint
new_vehicle.P_D_r = [   5.0  , 550.5,  165.0]; % push rod wishbone end
new_vehicle.P_E_r = [-280.0   , 259.7, 345.0]; % upper wishbone rear pivot
new_vehicle.P_F_r = [  80.0  , 259.8,  345.0]; % upper wishbone front pivot
new_vehicle.P_H_r = [  -5.0  , 560.9,  380.3]; % upper wishbone outer ball joint
new_vehicle.P_L_r = [ -61.0  , 248.8,  315.0]; % inner track rod ball joint
new_vehicle.P_M_r = [ -61.0  , 540.0,  344.0]; % outer track rod ball joint
new_vehicle.P_P_r = [ 100.0  , 210.0,  432.7]; % rocker axis 1st point
new_vehicle.P_Q_r = [ 100.0  , 310.0,  595.8]; % rocker axis 2nd point
new_vehicle.P_R_r = [ 158.0  , 250.0,  475.0]; % push rod rocker end
new_vehicle.P_S_r = [ 100.0  , 160.0,  530.0]; % damper to rocker point
new_vehicle.P_T_r = [-170.0   , 160.0, 530.0]; % damper to body point







































