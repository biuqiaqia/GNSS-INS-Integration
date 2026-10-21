function [TCKF_result,GNSS_llh] =...
    TCKF_INS_GNSS_GH(CPT_data,GNSS_data,IMU_data,TC_KF_config)

% ------------------------------------------------------------------------
% Main function for TCKF
% modified by Paul's example
%
% by GH.Zhang 2019/04/29
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------

D2R = pi/180;
R2D = 180/pi;
%% Begins
% Loosely coupled ECEF Inertial navigation and GNSS integrated navigation
% simulation
no_epochs_GT = size(CPT_data,1);
no_epochs_IMU = size(IMU_data,1);
no_epochs_GNSS = size(GNSS_data,1);
user_pos_xyz0 = llh2xyz(CPT_data(1,4:6).*[D2R,D2R,1]);

% True Velocity & Doppler
temp.GT_t = CPT_data(:,1);
for id = 1:1:no_epochs_GNSS
    [~,temp.GT_id] = min(abs(temp.GT_t-GNSS_data{id,12}));
    temp.V_enu = [CPT_data(temp.GT_id,8),CPT_data(temp.GT_id,7),-CPT_data(temp.GT_id,9)];
    temp.VR_xyz = enu2xyz(temp.V_enu,user_pos_xyz0);
    temp.V_xyz = temp.VR_xyz - user_pos_xyz0';
    GNSS_vel_error(id,1) = norm(temp.V_xyz-GNSS_data{id,13});
%     GNSS_data{id,13} = temp.V_xyz;
    temp.GNSS_doppler_mea = [GNSS_data{id,10}(:,1),GNSS_data{id,11}(:,1),...
                             GNSS_data{id,10}(:,3:5),GNSS_data{id,11}(:,2:5)];
    temp.Tdoppler = True_doppler_by_velocity(temp.GNSS_doppler_mea,temp.V_xyz,user_pos_xyz0);
%     GNSS_data{id,11}(:,1) = temp.Tdoppler;
end
% True Doppler

%% Initialize true navigation solution
% old_time = IMU_data(1,2)+IMU_data(1,3)*1e-9;
old_time = IMU_data(1,1);
true_L_b = CPT_data(1,4)*D2R;
true_lambda_b = CPT_data(1,5)*D2R;
true_h_b = CPT_data(1,6);
true_v_eb_n = CPT_data(1,7:9)';
true_eul_nb = CPT_data(1,10:12)';
true_C_b_n = Euler_to_CTM(true_eul_nb)';
[~,~,c_e_n,old_true_C_b_e] =...
    NED_to_ECEFdebug(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n); %coordinate transformation [what is NED?]
% old_est_C_b_e =old_true_C_b_e;

% Determine Least-squares GNSS position solution and use to initialize INS
[GNSS_r_eb_e,~,~,temp.clock] = Least_square_positioning(GNSS_data{1,10}(:,1:5),user_pos_xyz0);
temp.GNSS_doppler_mea = [GNSS_data{1,10}(:,1),GNSS_data{1,11}(:,1),...
                         GNSS_data{1,10}(:,3:5),GNSS_data{1,11}(:,2:5)];
[temp.VR,temp.dclock] = Least_square_velocity(temp.GNSS_doppler_mea,user_pos_xyz0);
temp.clock_info = [GNSS_data{1,10}(:,1),temp.clock,temp.dclock];
est_clock(1,:) = mean(temp.clock_info(temp.clock_info(:,1)<=32,2:3));
est_clock(2,:) = mean(temp.clock_info(temp.clock_info(:,1)>32,2:3));
GNSS_r_eb_e = GNSS_r_eb_e';
GNSS_v_eb_e = GNSS_data{1,13};
old_est_r_eb_e = GNSS_r_eb_e;
old_est_v_eb_e = GNSS_v_eb_e;
[old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n] =...
    pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e);
est_L_b = old_est_L_b;

% Initialize estimated attitude solution
old_est_C_b_n = Euler_to_CTM(CPT_data(1,10:12)'.*[D2R,D2R,D2R])';% insert imu 1st data (now CPT!!!)
[~,~,old_est_C_b_e] = NED_to_ECEF(old_est_L_b,...
    old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);

% Initialize output profile record and errors record
TCKF_result = zeros(no_epochs_IMU,10);

% Generate output profile record
TCKF_result(1,1) = old_time;
TCKF_result(1,2) = old_est_L_b*R2D;
TCKF_result(1,3) = old_est_lambda_b*R2D;
TCKF_result(1,4) = old_est_h_b;
TCKF_result(1,5:7) = old_est_v_eb_n';
TCKF_result(1,8:10) = CTM_to_Euler(old_est_C_b_n')'.*[R2D,R2D,R2D];

% Determine errors and generate output record
% [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
%     old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n,...
%     true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
% out_errors(1,1) = old_time;
% out_errors(1,2:4) = delta_r_eb_n';
% out_errors(1,5:7) = delta_v_eb_n';
% out_errors(1,8:10) = delta_eul_nb_n';

% Initialize Kalman filter P matrix and IMU bias states
P_matrix = Initialize_TC_P_matrix(TC_KF_config);
est_IMU_bias = zeros(6,1);

% Generate IMU bias and clock output records
out_IMU_bias_est(1,1) = old_time;
out_IMU_bias_est(1,2:7) = est_IMU_bias';

% Generate KF uncertainty record
out_KF_SD(1,1) = old_time;
for i =1:15
    out_KF_SD(1,i+1) = sqrt(P_matrix(i,i));
end % for i

% Initialize GNSS model timing
time_last_GNSS = old_time;
GNSS_epoch = 1;
GNSS_xyz(1,:) = GNSS_r_eb_e';
GNSS_llh(1,:) = xyz2llh(GNSS_r_eb_e);
GNSS_llh(1,:) = GNSS_llh(1,:).*[R2D,R2D,1];

%% Main loop
for epoch = 2:no_epochs_IMU
    
%     if epoch >= 22389
%         epoch = epoch+1;
%     end
    
    disp(['GNSS Epoch ---> ',num2str(GNSS_epoch),'/',num2str(no_epochs_GNSS)]);
    % Input data from motion profile
%     time = IMU_data(epoch,2)+IMU_data(epoch,3)*1e-9;
    time = IMU_data(epoch,1);

    % Time interval
    tor_i = time - old_time;
    
    % Correct IMU errors
    meas_f_ib_b = [IMU_data(epoch,5);...
                   IMU_data(epoch,4);...
                   -IMU_data(epoch,6)];
    meas_omega_ib_b = [IMU_data(epoch,8);...
                       IMU_data(epoch,7);...
                       -IMU_data(epoch,9)];
    meas_f_ib_b = meas_f_ib_b - est_IMU_bias(1:3);
    meas_omega_ib_b = meas_omega_ib_b - est_IMU_bias(4:6);
    
    % Update estimated navigation solution
%     [est_r_eb_e,est_v_eb_e,est_C_b_e] = Nav_equations_ECEF(tor_i,...
%         old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
%         meas_omega_ib_b);
    [est_r_eb_e,est_v_eb_e,est_C_b_e,f_ib_e_debug] = Nav_equations_ECEFdelete1(tor_i,...
        old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
        meas_omega_ib_b);
    
%     est_C_b_n = c_e_n * est_C_b_e;
    
%     est_v_n = c_e_n *est_v_eb_e;
%     f_ib_n_debug = c_e_n *f_ib_e_debug;
%     f_ib_b_debug = est_C_b_n' *f_ib_n_debug;

    % Determine whether to update GNSS simulation and run Kalman filter
    if  GNSS_epoch<size(GNSS_data,1) && time >= GNSS_data{GNSS_epoch+1,12}  && 1
        GNSS_epoch = GNSS_epoch + 1;
        temp.prn = GNSS_data{GNSS_epoch,10}(:,1);
        temp.nsv0 = 3+any(temp.prn>0&temp.prn<=32)+...
                    any(temp.prn>32&temp.prn<=56)+...
                    any(temp.prn>56&temp.prn<=86)+...
                    any(temp.prn>86&temp.prn<=123);
        
        if size(GNSS_data{GNSS_epoch,10},1)>=temp.nsv0||1
            tor_s = time - time_last_GNSS;  % KF time interval
            time_last_GNSS = time;

            % Determine GNSS position solution
            [GNSS_r_eb_e,~,~,~] = WLS_positioning(GNSS_data{GNSS_epoch,10}(:,1:5),...
                                                  GNSS_data{GNSS_epoch,9}(:,2),...
                                                  GNSS_data{GNSS_epoch,9}(:,1),...
                                                  user_pos_xyz0);
            GNSS_xyz(GNSS_epoch,:) = GNSS_r_eb_e';
            GNSS_llh(GNSS_epoch,:) = xyz2llh(GNSS_r_eb_e);
            GNSS_llh(GNSS_epoch,:) = GNSS_llh(GNSS_epoch,:).*[R2D,R2D,1];
            GNSS_r_eb_e = GNSS_r_eb_e';
            GNSS_v_eb_e = GNSS_data{GNSS_epoch,13};
            GNSS_measurements = [GNSS_data{GNSS_epoch,10}(:,2),...          %pseudorange
                                 -GNSS_data{GNSS_epoch,11}(:,1).*GNSS_data{GNSS_epoch,11}(:,5),...%pseudorange rate
                                 GNSS_data{GNSS_epoch,10}(:,3:5),...        %SV position
                                 GNSS_data{GNSS_epoch,11}(:,2:4),...        %SV velocity
                                 GNSS_data{GNSS_epoch,10}(:,1)];            %PRN
%             GNSS_measurements(:,1) = GNSS_measurements(:,1) - (GNSS_measurements(:,9)>32).*17;

            % Run Integration Kalman filter
            [est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix] =...
                TC_KF_Epoch(GNSS_measurements,tor_s,est_C_b_e,...
                est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                meas_f_ib_b,est_L_b,TC_KF_config);

            % Generate IMU bias and clock output records
            out_IMU_bias_est(GNSS_epoch,1) = time;
            out_IMU_bias_est(GNSS_epoch,2:7) = est_IMU_bias';
    %         out_clock(GNSS_epoch,1) = time;
    %         out_clock(GNSS_epoch,2:3) = est_clock;

            % Generate KF uncertainty output record
            out_KF_SD(GNSS_epoch,1) = time;
            for i =1:17
                out_KF_SD(GNSS_epoch,i+1) = sqrt(P_matrix(i,i));
            end % for i
        else
            GNSS_llh(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_xyz(GNSS_epoch,:) = [nan,nan,nan];
        end
    end % if time    
    
    % Convert navigation solution to NED
    [est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n] =...
        ECEF_to_NED(est_r_eb_e,est_v_eb_e,est_C_b_e);
    
    est_C_b_n = c_e_n * est_C_b_e;
    est_v_b = (est_C_b_n') * (c_e_n*est_v_eb_e);

    % Generate output profile record
    TCKF_result(epoch,1) = time;
    TCKF_result(epoch,2) = est_L_b*R2D;
    TCKF_result(epoch,3) = est_lambda_b*R2D;
    TCKF_result(epoch,4) = est_h_b;
    TCKF_result(epoch,5:7) = est_v_eb_n';
    TCKF_result(epoch,8:10) = CTM_to_Euler(est_C_b_n')'.*[R2D,R2D,R2D];
    if TCKF_result(epoch,10) < 0
        TCKF_result(epoch,10) = TCKF_result(epoch,10)+360;
    end
    % Determine errors and generate output record
%     true_L_b = in_profile(epoch,2);
%     true_lambda_b = in_profile(epoch,3);
%     true_h_b = in_profile(epoch,4);
%     true_v_eb_n = in_profile(epoch,5:7)';
%     true_eul_nb = in_profile(epoch,8:10)';
%     true_C_b_n = Euler_to_CTM(true_eul_nb)';
%     [true_r_eb_e,true_v_eb_e,true_C_b_e] =...
%         NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
%     [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
%         est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
%         true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
%     out_errors(epoch,1) = time;
%     out_errors(epoch,2:4) = delta_r_eb_n';
%     out_errors(epoch,5:7) = delta_v_eb_n';
%     out_errors(epoch,8:10) = delta_eul_nb_n';
    
    % Reset old values
    old_time = time;
    old_est_r_eb_e = est_r_eb_e;
    old_est_v_eb_e = est_v_eb_e;
    old_est_C_b_e = est_C_b_e;
end