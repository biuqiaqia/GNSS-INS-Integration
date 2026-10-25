function [TCKF_result,TCKF_result_type,TCKF_result_uncert,GNSS_llh,GNSS_llh_1,GNSS_llh_2] =...
    TCKF_INS_GNSS_GH_1(CPT_data,GNSS_data,IMU_data,TC_KF_config)

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
% simulation
no_epochs_GT = size(CPT_data,1);
no_epochs_IMU = size(IMU_data,1);
no_epochs_GNSS = size(GNSS_data,1);
user_pos_xyz0 = llh2xyz(CPT_data(1,4:6).*[D2R,D2R,1]);

% True Velocity & Doppler 【看一下这里是不是可以不要啊】
temp.GT_t = CPT_data(:,1);
for id = 1:1:no_epochs_GNSS
    % [~,temp.GT_id] = min(abs(temp.GT_t-GNSS_data{id,12}));
    % temp.V_enu = [CPT_data(temp.GT_id,8),CPT_data(temp.GT_id,7),-CPT_data(temp.GT_id,9)];
    % temp.VR_xyz = enu2xyz(temp.V_enu,user_pos_xyz0);
    % temp.V_xyz = temp.VR_xyz - user_pos_xyz0';
    % GNSS_vel_error(id,1) = norm(temp.V_xyz-GNSS_data{id,13});
%     GNSS_data{id,13} = temp.V_xyz;
    temp.GNSS_doppler_mea = [GNSS_data{id,10}(:,1),GNSS_data{id,11}(:,1),...
                             GNSS_data{id,10}(:,3:5),GNSS_data{id,11}(:,2:5)];
    % temp.Tdoppler = True_doppler_by_velocity(temp.GNSS_doppler_mea,temp.V_xyz,user_pos_xyz0);
%     GNSS_data{id,11}(:,1) = temp.Tdoppler;
end
% True Doppler

%% Initialize true navigation solution 
% old_time = IMU_data(1,2)+IMU_data(1,3)*1e-9;%什么意思

old_time = IMU_data(1,1);%时间
true_L_b = CPT_data(1,4)*D2R;%纬度
true_lambda_b = CPT_data(1,5)*D2R;%经度
true_h_b = CPT_data(1,6);%高度
true_v_eb_n = CPT_data(1,7:9)';%从地球指向载体（body）的速度，在 n 系（NED）中表示【都为0】
true_eul_nb = CPT_data(1,10:12)'*D2R;%n 系到 b 系的欧拉角（姿态角）【改成了rad】
true_C_b_n = Euler_to_CTM(true_eul_nb)';%真实 roll pitch yaw (rad)
[~,~,c_e_n,old_true_C_b_e] =...
    NED_to_ECEFdebug(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n); %coordinate transformation 此处主要是获取了C_e_n 其他的矩阵没用
% old_est_C_b_e =old_true_C_b_e;

% Determine Least-squares GNSS position solution and use to initialize INS
% 都只有一行
[GNSS_r_eb_e,~,~,temp.clock] = Least_square_positioning(GNSS_data{1,10}(:,1:5),user_pos_xyz0);

temp.GNSS_doppler_mea = [GNSS_data{1,10}(:,1),GNSS_data{1,11}(:,1),...
                         GNSS_data{1,10}(:,3:5),GNSS_data{1,11}(:,2:5)];%输入是对的
[temp.VR,temp.dclock] = Least_square_velocity(temp.GNSS_doppler_mea,user_pos_xyz0);

%钟差  由于目前只有GPS，所以est clock只有一行
temp.clock_info = [GNSS_data{1,10}(:,1),temp.clock,temp.dclock];
est_clock(1,:) = mean(temp.clock_info(temp.clock_info(:,1)<=32,2:3));
est_clock(2,:) = mean(temp.clock_info(temp.clock_info(:,1)>32,2:3));

%Type’s
est_clock_1(1,:) = mean(temp.clock_info(temp.clock_info(:,1)<=32,2:3));
est_clock_1(2,:) = mean(temp.clock_info(temp.clock_info(:,1)>32,2:3));
%Uncertainty’s
est_clock_2(1,:) = mean(temp.clock_info(temp.clock_info(:,1)<=32,2:3));
est_clock_2(2,:) = mean(temp.clock_info(temp.clock_info(:,1)>32,2:3));

% 这里在做什么？为什么est pos要放两个？为什么vr要放两个不一样来源的？
%好像一个是用来初始化的，一个是后来main loop循环更新用的？
GNSS_r_eb_e = GNSS_r_eb_e';%ECEF position est
GNSS_v_eb_e = GNSS_data{1,13};%VRx VRy VRz
old_est_r_eb_e = GNSS_r_eb_e;%ECEF position est
old_est_v_eb_e = GNSS_v_eb_e;%VRx VRy VRz

%Type’s
old_est_r_eb_e_1 = GNSS_r_eb_e;%ECEF position est
old_est_v_eb_e_1 = GNSS_v_eb_e;%VRx VRy VRz
%Uncertainty’s
old_est_r_eb_e_2 = GNSS_r_eb_e;%ECEF position est
old_est_v_eb_e_2 = GNSS_v_eb_e;%VRx VRy VRz

[old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n] =...
    pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e);
est_L_b = old_est_L_b;%latitude 这里为什么要用两个值？%好像一个是用来初始化的，一个是后来main loop循环更新用的？
est_L_b_1 = old_est_L_b;
est_L_b_2 = old_est_L_b;
% Initialize estimated attitude solution
old_est_C_b_n = Euler_to_CTM(CPT_data(1,10:12)'.*[D2R,D2R,D2R])';% insert imu 1st data (now CPT!!!)  Body to NED matrix
[~,~,old_est_C_b_e] = NED_to_ECEF(old_est_L_b,...
    old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);% NED to ECEF

% get body to ECEF  old_est_C_b_e

%Type’s
[~,~,old_est_C_b_e_1] = NED_to_ECEF(old_est_L_b,...
    old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);

%Uncertainty’s
[~,~,old_est_C_b_e_2] = NED_to_ECEF(old_est_L_b,...
    old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);


% Initialize output profile record and errors record
TCKF_result = zeros(no_epochs_IMU,10);

TCKF_result_type = zeros(no_epochs_IMU,10);

TCKF_result_uncert = zeros(no_epochs_IMU,10);


% Generate output profile record 这里都是初始值，所以三个一样。这里的经纬度是GNSS data第一行最小二乘的，不是CPT。速度同样
TCKF_result(1,1) = old_time;
TCKF_result(1,2) = old_est_L_b*R2D; %latitude deg
TCKF_result(1,3) = old_est_lambda_b*R2D; %longitude deg
TCKF_result(1,4) = old_est_h_b; %height
TCKF_result(1,5:7) = old_est_v_eb_n'; %V
TCKF_result(1,8:10) = CTM_to_Euler(old_est_C_b_n')'.*[R2D,R2D,R2D]; %Euler

TCKF_result_type(1,1) = old_time;
TCKF_result_type(1,2) = old_est_L_b*R2D; %latitude deg
TCKF_result_type(1,3) = old_est_lambda_b*R2D; %longitude deg
TCKF_result_type(1,4) = old_est_h_b;  %height
TCKF_result_type(1,5:7) = old_est_v_eb_n';  %V
TCKF_result_type(1,8:10) = CTM_to_Euler(old_est_C_b_n')'.*[R2D,R2D,R2D]; %Euler

TCKF_result_uncert(1,1) = old_time;
TCKF_result_uncert(1,2) = old_est_L_b*R2D; %latitude deg
TCKF_result_uncert(1,3) = old_est_lambda_b*R2D; %longitude deg
TCKF_result_uncert(1,4) = old_est_h_b;  %height
TCKF_result_uncert(1,5:7) = old_est_v_eb_n';  %V
TCKF_result_uncert(1,8:10) = CTM_to_Euler(old_est_C_b_n')'.*[R2D,R2D,R2D]; %Euler

% % Determine errors and generate output record 这里没有Calculate_errors_NED的函数代码
% [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
%     old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n,...
%     true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
% out_errors(1,1) = old_time;
% out_errors(1,2:4) = delta_r_eb_n';
% out_errors(1,5:7) = delta_v_eb_n';
% out_errors(1,8:10) = delta_eul_nb_n';

% Initialize Kalman filter P matrix and IMU bias states
% Origin
P_matrix = Initialize_TC_P_matrix(TC_KF_config);
% Type's
P_matrix_1 = Initialize_TC_P_matrix(TC_KF_config);
% Uncertainty's
P_matrix_2 = Initialize_TC_P_matrix(TC_KF_config);

est_IMU_bias = zeros(6,1);

est_IMU_bias_1 = zeros(6,1);

est_IMU_bias_2 = zeros(6,1);

% Generate IMU bias and clock output records
out_IMU_bias_est(1,1) = old_time;
out_IMU_bias_est(1,2:7) = est_IMU_bias';

out_IMU_bias_est_1(1,1) = old_time;
out_IMU_bias_est_1(1,2:7) = est_IMU_bias_1';

out_IMU_bias_est_2(1,1) = old_time;
out_IMU_bias_est_2(1,2:7) = est_IMU_bias_2';

% Generate KF uncertainty record 
out_KF_SD(1,1) = old_time;
out_KF_SD_1(1,1) = old_time;
out_KF_SD_2(1,1) = old_time;
for i =1:15
    out_KF_SD(1,i+1) = sqrt(P_matrix(i,i));
    out_KF_SD_1(1,i+1) = sqrt(P_matrix_1(i,i));
    out_KF_SD_2(1,i+1) = sqrt(P_matrix_2(i,i));
end % for i

% Initialize GNSS model timing 初始的，三个都一样
time_last_GNSS = old_time;
GNSS_epoch = 1;
GNSS_xyz(1,:) = GNSS_r_eb_e'; %est pos ecef
GNSS_llh(1,:) = xyz2llh(GNSS_r_eb_e); %est pos llh rad
GNSS_llh(1,:) = GNSS_llh(1,:).*[R2D,R2D,1]; %est pos llh deg

GNSS_xyz_1(1,:) = GNSS_r_eb_e';
GNSS_llh_1(1,:) = xyz2llh(GNSS_r_eb_e);
GNSS_llh_1(1,:) = GNSS_llh_1(1,:).*[R2D,R2D,1];

GNSS_xyz_2(1,:) = GNSS_r_eb_e';
GNSS_llh_2(1,:) = xyz2llh(GNSS_r_eb_e);
GNSS_llh_2(1,:) = GNSS_llh_2(1,:).*[R2D,R2D,1];

%% Main loop
for epoch = 2:no_epochs_IMU

    disp(['GNSS Epoch ---> ',num2str(GNSS_epoch),'/',num2str(no_epochs_GNSS)]);
    % Input data from motion profile 这是什么？
%     time = IMU_data(epoch,2)+IMU_data(epoch,3)*1e-9;
    time = IMU_data(epoch,1); %IMU现在的时间
    % if GNSS_epoch >= 143
    %     keyboard;
    % end
    % Time interval
    tor_i = time - old_time;%IMU现在的时间减去 IMU上一个时间（后面有对old time 更新）
    
    % Correct IMU errors 这里为什么是54-6，87-9【可能错了】 Mti-320是 前左上  【改过了】
    % origin
    meas_f_ib_b = [IMU_data(epoch,4);...
                   -IMU_data(epoch,5);...
                   -IMU_data(epoch,6)];
    meas_f_ib_b = meas_f_ib_b - est_IMU_bias(1:3); % 比力

    meas_omega_ib_b = [IMU_data(epoch,7);...
                       -IMU_data(epoch,8);...
                       -IMU_data(epoch,9)];
    meas_omega_ib_b = meas_omega_ib_b - est_IMU_bias(4:6); % 陀螺

    % type
    meas_f_ib_b_1 = [IMU_data(epoch,4);...
                   -IMU_data(epoch,5);...
                   -IMU_data(epoch,6)];
    meas_f_ib_b_1 = meas_f_ib_b_1 - est_IMU_bias_1(1:3);
    meas_omega_ib_b_1 = [IMU_data(epoch,7);...
                       -IMU_data(epoch,8);...
                       -IMU_data(epoch,9)];
    meas_omega_ib_b_1 = meas_omega_ib_b_1 - est_IMU_bias_1(4:6);

    % uncertainty
    meas_f_ib_b_2 = [IMU_data(epoch,4);...
                   -IMU_data(epoch,5);...
                   -IMU_data(epoch,6)];
    meas_f_ib_b_2 = meas_f_ib_b_2 - est_IMU_bias(1:3);
    meas_omega_ib_b_2 = [IMU_data(epoch,7);...
                       -IMU_data(epoch,8);...
                       -IMU_data(epoch,9)];
    meas_omega_ib_b_2 = meas_omega_ib_b_2 - est_IMU_bias_2(4:6);
    
    % Update estimated navigation solution
%     [est_r_eb_e,est_v_eb_e,est_C_b_e] = Nav_equations_ECEF(tor_i,...
%         old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
%         meas_omega_ib_b);

    % 这里xyz究竟是什么顺序，非常重要。【破案了，目前感觉应该就是我数据收的问题，导致yaw错误比较严重，有些姿态角收敛不了】
    % ECEF大转换
    [est_r_eb_e,est_v_eb_e,est_C_b_e,f_ib_e_debug] = Nav_equations_ECEFdelete1(tor_i,...
        old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
        meas_omega_ib_b);

    [est_r_eb_e_1,est_v_eb_e_1,est_C_b_e_1,~] = Nav_equations_ECEFdelete1(tor_i,...
        old_est_r_eb_e_1,old_est_v_eb_e_1,old_est_C_b_e_1,meas_f_ib_b_1,...
        meas_omega_ib_b_1); 
    
    [est_r_eb_e_2,est_v_eb_e_2,est_C_b_e_2,~] = Nav_equations_ECEFdelete1(tor_i,...
        old_est_r_eb_e_2,old_est_v_eb_e_2,old_est_C_b_e_2,meas_f_ib_b_2,...
        meas_omega_ib_b_2);

%     est_C_b_n = c_e_n * est_C_b_e;
    
%     est_v_n = c_e_n *est_v_eb_e;
%     f_ib_n_debug = c_e_n *f_ib_e_debug;
%     f_ib_b_debug = est_C_b_n' *f_ib_n_debug;

    % Determine whether to update GNSS simulation and run Kalman filter
    % 这里要看哦，要理解。前面的还有有点模糊的地方，别忘了

    if  GNSS_epoch<size(GNSS_data,1) && time >= GNSS_data{GNSS_epoch+1,12}  %&& 1 %这里的time是IMU的time所以要判断time刚大于或等于下一个gnss epoch的时候，这个时候要做更新。&&1（和1），该去掉的
        GNSS_epoch = GNSS_epoch + 1;
        temp.prn = GNSS_data{GNSS_epoch,10}(:,1);
        temp.nsv0 = 3+any(temp.prn>0&temp.prn<=32)+...
                    any(temp.prn>32&temp.prn<=56)+...
                    any(temp.prn>56&temp.prn<=86)+...
                    any(temp.prn>86&temp.prn<=123);%这里是有检测到卫星就加一为4
        
        if size(GNSS_data{GNSS_epoch,10},1)>=temp.nsv0 %||1 %这里判断epoch的卫星数有咩有大于等于4，||1（或1）恒为True，该去掉的
            tor_s = time - time_last_GNSS;  % KF time interval
            time_last_GNSS = time;

            %% Origin
            % Determine GNSS position solution
            [GNSS_r_eb_e,~,~,~] = WLS_positioning(GNSS_data{GNSS_epoch,10}(:,1:5),...
                                                  GNSS_data{GNSS_epoch,9}(:,2),...
                                                  GNSS_data{GNSS_epoch,9}(:,1),...
                                                  user_pos_xyz0); %输入正确
            GNSS_xyz(GNSS_epoch,:) = GNSS_r_eb_e';
            GNSS_llh(GNSS_epoch,:) = xyz2llh(GNSS_r_eb_e);
            GNSS_llh(GNSS_epoch,:) = GNSS_llh(GNSS_epoch,:).*[R2D,R2D,1];
            GNSS_r_eb_e = GNSS_r_eb_e';
            %% Type

            % Determine GNSS position solution
            [GNSS_r_eb_e_1,~,~,~] = WLS_positioning_1(GNSS_data{GNSS_epoch,10}(:,1:5),...
                                                     GNSS_data{GNSS_epoch,9}(:,2),...
                                                  GNSS_data{GNSS_epoch,9}(:,1),...
                                                  GNSS_data{GNSS_epoch,15},...
                                                  user_pos_xyz0);
            GNSS_xyz_1(GNSS_epoch,:) = GNSS_r_eb_e_1';
            GNSS_llh_1(GNSS_epoch,:) = xyz2llh(GNSS_r_eb_e_1);
            GNSS_llh_1(GNSS_epoch,:) = GNSS_llh_1(GNSS_epoch,:).*[R2D,R2D,1];
            GNSS_r_eb_e_1 = GNSS_r_eb_e_1';
            
            %% Uncertainty

            % Determine GNSS position solution
            [GNSS_r_eb_e_2,~,~,~] = WLS_positioning_2(GNSS_data{GNSS_epoch,10}(:,1:5),...
                                                    GNSS_data{GNSS_epoch,9}(:,2),...
                                                  GNSS_data{GNSS_epoch,9}(:,1),...
                                                  GNSS_data{GNSS_epoch,16},...
                                                  user_pos_xyz0);
            GNSS_xyz_2(GNSS_epoch,:) = GNSS_r_eb_e_2';
            GNSS_llh_2(GNSS_epoch,:) = xyz2llh(GNSS_r_eb_e_2);
            GNSS_llh_2(GNSS_epoch,:) = GNSS_llh_2(GNSS_epoch,:).*[R2D,R2D,1];
            GNSS_r_eb_e_2 = GNSS_r_eb_e_2';
            

            GNSS_v_eb_e = GNSS_data{GNSS_epoch,13};%【意义是啥？VRxyz】
            GNSS_measurements = [GNSS_data{GNSS_epoch,10}(:,2),...          %pseudorange
                                 -GNSS_data{GNSS_epoch,11}(:,1).*GNSS_data{GNSS_epoch,11}(:,5),...%pseudorange rate
                                 GNSS_data{GNSS_epoch,10}(:,3:5),...        %SV position
                                 GNSS_data{GNSS_epoch,11}(:,2:4),...        %SV velocity
                                 GNSS_data{GNSS_epoch,10}(:,1),...          %PRN
                                 GNSS_data{GNSS_epoch,15},...               %Type
                                 GNSS_data{GNSS_epoch,16},...               %Uncertainty  
                                 GNSS_data{GNSS_epoch,9}];                  %CNR, Elev 
%             GNSS_measurements(:,1) = GNSS_measurements(:,1) - (GNSS_measurements(:,9)>32).*17;

            
            %% Origin   
            % Run Integration Kalman filter
            [est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix] =...
                TC_KF_Epoch(GNSS_measurements,tor_s,est_C_b_e,...
                est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                meas_f_ib_b,est_L_b,TC_KF_config);
            %% Type
            % Run Integration Kalman filter
            [est_C_b_e_1,est_v_eb_e_1,est_r_eb_e_1,est_IMU_bias_1,est_clock_1,P_matrix_1] =...
                TC_KF_Epoch_1(GNSS_measurements,tor_s,est_C_b_e_1,...
                est_v_eb_e_1,est_r_eb_e_1,est_IMU_bias_1,est_clock_1,P_matrix_1,...
                meas_f_ib_b_1,est_L_b_1,TC_KF_config);
            %% Uncertainty
            % Run Integration Kalman filter
            [est_C_b_e_2,est_v_eb_e_2,est_r_eb_e_2,est_IMU_bias_2,est_clock_2,P_matrix_2] =...
                TC_KF_Epoch_2(GNSS_measurements,tor_s,est_C_b_e_2,...
                est_v_eb_e_2,est_r_eb_e_2,est_IMU_bias_2,est_clock_2,P_matrix_2,...
                meas_f_ib_b_2,est_L_b_2,TC_KF_config);


            %%
            % Generate IMU bias and clock output records
            out_IMU_bias_est(GNSS_epoch,1) = time;
            out_IMU_bias_est(GNSS_epoch,2:7) = est_IMU_bias';

            out_IMU_bias_est_1(GNSS_epoch,1) = time;
            out_IMU_bias_est_1(GNSS_epoch,2:7) = est_IMU_bias_1';

            out_IMU_bias_est_2(GNSS_epoch,1) = time;
            out_IMU_bias_est_2(GNSS_epoch,2:7) = est_IMU_bias_2';
    %         out_clock(GNSS_epoch,1) = time;
    %         out_clock(GNSS_epoch,2:3) = est_clock;

            % Generate KF uncertainty output record
            out_KF_SD(GNSS_epoch,1) = time;
            out_KF_SD_1(GNSS_epoch,1) = time;
            out_KF_SD_2(GNSS_epoch,1) = time;
            for i =1:17
                out_KF_SD(GNSS_epoch,i+1) = sqrt(P_matrix(i,i));
                out_KF_SD_1(GNSS_epoch,i+1) = sqrt(P_matrix_1(i,i));
                out_KF_SD_2(GNSS_epoch,i+1) = sqrt(P_matrix_2(i,i));
            end % for i
        else
            GNSS_llh(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_xyz(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_llh_1(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_xyz_1(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_llh_2(GNSS_epoch,:) = [nan,nan,nan];
            GNSS_xyz_2(GNSS_epoch,:) = [nan,nan,nan];
        end
    end % if time    
    
    % Convert navigation solution to NED
    [est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n] =...
        ECEF_to_NED(est_r_eb_e,est_v_eb_e,est_C_b_e);

    [est_L_b_1,est_lambda_b_1,est_h_b_1,est_v_eb_n_1,est_C_b_n_1] =...
        ECEF_to_NED(est_r_eb_e_1,est_v_eb_e_1,est_C_b_e_1);

    [est_L_b_2,est_lambda_b_2,est_h_b_2,est_v_eb_n_2,est_C_b_n_2] =...
        ECEF_to_NED(est_r_eb_e_2,est_v_eb_e_2,est_C_b_e_2);
    
    est_C_b_n = c_e_n * est_C_b_e;
    est_v_b = (est_C_b_n') * (c_e_n*est_v_eb_e);

    est_C_b_n_1 = c_e_n * est_C_b_e_1;
    est_v_b_1 = (est_C_b_n_1') * (c_e_n*est_v_eb_e_1);

    est_C_b_n_2 = c_e_n * est_C_b_e_2;
    est_v_b_2 = (est_C_b_n_2') * (c_e_n*est_v_eb_e_2);

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

     % Generate output profile record
    TCKF_result_type(epoch,1) = time;
    TCKF_result_type(epoch,2) = est_L_b_1*R2D;
    TCKF_result_type(epoch,3) = est_lambda_b_1*R2D;
    TCKF_result_type(epoch,4) = est_h_b_1;
    TCKF_result_type(epoch,5:7) = est_v_eb_n_1';
    TCKF_result_type(epoch,8:10) = CTM_to_Euler(est_C_b_n_1')'.*[R2D,R2D,R2D];
    if TCKF_result_type(epoch,10) < 0
        TCKF_result_type(epoch,10) = TCKF_result_type(epoch,10)+360;
    end

     % Generate output profile record
    TCKF_result_uncert(epoch,1) = time;
    TCKF_result_uncert(epoch,2) = est_L_b_2*R2D;
    TCKF_result_uncert(epoch,3) = est_lambda_b_2*R2D;
    TCKF_result_uncert(epoch,4) = est_h_b_2;
    TCKF_result_uncert(epoch,5:7) = est_v_eb_n_2';
    TCKF_result_uncert(epoch,8:10) = CTM_to_Euler(est_C_b_n_2')'.*[R2D,R2D,R2D];
    if TCKF_result_uncert(epoch,10) < 0
        TCKF_result_uncert(epoch,10) = TCKF_result_uncert(epoch,10)+360;
    end

    %% 目前没有这么详细的真值，算不了
    % Determine errors and generate output record
    % true_L_b = in_profile(epoch,2);
    % true_lambda_b = in_profile(epoch,3);
    % true_h_b = in_profile(epoch,4);
    % true_v_eb_n = in_profile(epoch,5:7)';
    % true_eul_nb = in_profile(epoch,8:10)';
    % true_C_b_n = Euler_to_CTM(true_eul_nb)';
    % [true_r_eb_e,true_v_eb_e,true_C_b_e] =...
    %     NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
    % [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
    %     est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
    %     true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
    % out_errors(epoch,1) = time;
    % out_errors(epoch,2:4) = delta_r_eb_n';
    % out_errors(epoch,5:7) = delta_v_eb_n';
    % out_errors(epoch,8:10) = delta_eul_nb_n';
    %%
    % Reset old values
    old_time = time;
    old_est_r_eb_e = est_r_eb_e;
    old_est_v_eb_e = est_v_eb_e;
    old_est_C_b_e = est_C_b_e;

    % Reset old values
    old_est_r_eb_e_1 = est_r_eb_e_1;
    old_est_v_eb_e_1 = est_v_eb_e_1;
    old_est_C_b_e_1 = est_C_b_e_1;

    % Reset old values
    old_est_r_eb_e_2 = est_r_eb_e_2;
    old_est_v_eb_e_2 = est_v_eb_e_2;
    old_est_C_b_e_2 = est_C_b_e_2;
end

%% Visualization


Bias = 0;
if Bias == 1
    % 可视化bias
    time = out_IMU_bias_est(:,1);
    bias_ax = out_IMU_bias_est(:,2);
    bias_ay = out_IMU_bias_est(:,3);
    bias_az = out_IMU_bias_est(:,4);
    bias_gx = out_IMU_bias_est(:,5);
    bias_gy = out_IMU_bias_est(:,6);
    bias_gz = out_IMU_bias_est(:,7);
    
    bias_gx_dps = bias_gx * 180/pi;
    bias_gy_dps = bias_gy * 180/pi;
    bias_gz_dps = bias_gz * 180/pi;
    
    % 创建子图
    figure('Position', [100, 100, 1000, 800]);
    
    % 子图1：加速度计偏置
    subplot(2,1,1);
    plot(time, bias_ax, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_ay, 'g', 'LineWidth', 1.5);
    plot(time, bias_az, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_ax (x)', 'bias\_ay (y)', 'bias\_az (z)', 'Location', 'best');
    title('IMU Accelerometer Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (m/s^2)');
    xlim([min(time) max(time)]);
    
    % 子图2：陀螺仪偏置
    subplot(2,1,2);
    plot(time, bias_gx_dps, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_gy_dps, 'g', 'LineWidth', 1.5);
    plot(time, bias_gz_dps, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_gx (x)', 'bias\_gy (y)', 'bias\_gz (z)', 'Location', 'best');
    title('IMU Gyroscope Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (°/s)');
    xlim([min(time) max(time)]);
    
    
    sgtitle('TCKF IMU Bias Estimation Over Time', 'FontSize', 14, 'FontWeight', 'bold');
    
    
    % Type
    time = out_IMU_bias_est_1(:,1);
    bias_ax = out_IMU_bias_est_1(:,2);
    bias_ay = out_IMU_bias_est_1(:,3);
    bias_az = out_IMU_bias_est_1(:,4);
    bias_gx = out_IMU_bias_est_1(:,5);
    bias_gy = out_IMU_bias_est_1(:,6);
    bias_gz = out_IMU_bias_est_1(:,7);
    
    bias_gx_dps = bias_gx * 180/pi;
    bias_gy_dps = bias_gy * 180/pi;
    bias_gz_dps = bias_gz * 180/pi;
    
    % 创建子图
    figure('Position', [100, 100, 1000, 800]);
    
    % 子图1：加速度计偏置
    subplot(2,1,1);
    plot(time, bias_ax, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_ay, 'g', 'LineWidth', 1.5);
    plot(time, bias_az, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_ax (x)', 'bias\_ay (y)', 'bias\_az (z)', 'Location', 'best');
    title('IMU Accelerometer Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (m/s^2)');
    xlim([min(time) max(time)]);
    
    % 子图2：陀螺仪偏置
    subplot(2,1,2);
    plot(time, bias_gx_dps, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_gy_dps, 'g', 'LineWidth', 1.5);
    plot(time, bias_gz_dps, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_gx (x)', 'bias\_gy (y)', 'bias\_gz (z)', 'Location', 'best');
    title('IMU Gyroscope Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (°/s)');
    xlim([min(time) max(time)]);
    
    
    sgtitle('TCKF IMU Bias Estimation (Type) Over Time', 'FontSize', 14, 'FontWeight', 'bold');
    
    
    % Uncert
    time = out_IMU_bias_est_2(:,1);
    bias_ax = out_IMU_bias_est_2(:,2);
    bias_ay = out_IMU_bias_est_2(:,3);
    bias_az = out_IMU_bias_est_2(:,4);
    bias_gx = out_IMU_bias_est_2(:,5);
    bias_gy = out_IMU_bias_est_2(:,6);
    bias_gz = out_IMU_bias_est_2(:,7);
    
    bias_gx_dps = bias_gx * 180/pi;
    bias_gy_dps = bias_gy * 180/pi;
    bias_gz_dps = bias_gz * 180/pi;
    
    % 创建子图
    figure('Position', [100, 100, 1000, 800]);
    
    % 子图1：加速度计偏置
    subplot(2,1,1);
    plot(time, bias_ax, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_ay, 'g', 'LineWidth', 1.5);
    plot(time, bias_az, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_ax (x)', 'bias\_ay (y)', 'bias\_az (z)', 'Location', 'best');
    title('IMU Accelerometer Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (m/s^2)');
    xlim([min(time) max(time)]);
    
    % 子图2：陀螺仪偏置
    subplot(2,1,2);
    plot(time, bias_gx_dps, 'r', 'LineWidth', 1.5); hold on;
    plot(time, bias_gy_dps, 'g', 'LineWidth', 1.5);
    plot(time, bias_gz_dps, 'b', 'LineWidth', 1.5);
    grid on; hold off;
    legend('bias\_gx (x)', 'bias\_gy (y)', 'bias\_gz (z)', 'Location', 'best');
    title('IMU Gyroscope Bias Estimation');
    xlabel('Time (s)');
    ylabel('Bias (°/s)');
    xlim([min(time) max(time)]);
    
    
    sgtitle('TCKF IMU Bias Estimation (Uncert) Over Time', 'FontSize', 14, 'FontWeight', 'bold');
end

% % 姿态误差标准差（Roll/Pitch/Yaw）
% figure('Position', [100, 100, 1000, 600]);
% time = out_KF_SD(:,1);
% 
% att_std = out_KF_SD(:,2:4) * 180/pi;        % °
% att_std_1 = out_KF_SD_1(:,2:4) * 180/pi;
% att_std_2 = out_KF_SD_2(:,2:4) * 180/pi;
% 
% plot(time, att_std, 'LineWidth',1.2); hold on;
% plot(time, att_std_1, '--', 'LineWidth',1.2);
% plot(time, att_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% 
% title('Attitude Error STD Comparison');
% xlabel('Time (s)'); ylabel('STD (°)');
% legend('Roll (Baseline)', 'Pitch (Baseline)', 'Yaw (Baseline)', ...
%        'Roll (Type)', 'Pitch (Type)', 'Yaw (Type)', ...
%        'Roll (Uncert)', 'Pitch (Uncert)', 'Yaw (Uncert)', ...
%        'Location', 'best');
% 
% % 速度误差标准差（ECEF）
% figure('Position', [100, 100, 1000, 600]);
% vel_std = out_KF_SD(:,5:7);
% vel_std_1 = out_KF_SD_1(:,5:7);
% vel_std_2 = out_KF_SD_2(:,5:7);
% 
% plot(time, vel_std, 'LineWidth',1.2); hold on;
% plot(time, vel_std_1, '--', 'LineWidth',1.2);
% plot(time, vel_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% 
% title('Velocity Error STD Comparison');
% xlabel('Time (s)'); ylabel('STD (m/s)');
% legend('East (Baseline)', 'North (Baseline)', 'Up (Baseline)', ...
%        'East (Type)', 'North (Type)', 'Up (Type)', ...
%        'East (Uncert)', 'North (Uncert)', 'Up (Uncert)', ...
%        'Location', 'best');
% 
% % 位置误差标准差（ECEF）
% figure('Position', [100, 100, 1000, 600]);
% pos_std = out_KF_SD(:,8:10);
% pos_std_1 = out_KF_SD_1(:,8:10);
% pos_std_2 = out_KF_SD_2(:,8:10);
% 
% plot(time, pos_std, 'LineWidth',1.2); hold on;
% plot(time, pos_std_1, '--', 'LineWidth',1.2);
% plot(time, pos_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% 
% title('Position Error STD Comparison');
% xlabel('Time (s)'); ylabel('STD (m)');
% legend('X_e (Baseline)', 'Y_e (Baseline)', 'Z_e (Baseline)', ...
%        'X_e (Type)', 'Y_e (Type)', 'Z_e (Type)', ...
%        'X_e (Uncert)', 'Y_e (Uncert)', 'Z_e (Uncert)', ...
%        'Location', 'best');
% 
% % IMU 零偏不确定性（加计 + 陀螺）
% figure('Position', [100, 100, 1000, 800]);
% 
% % 加计零偏 STD (m/s²)
% subplot(2,1,1);
% acc_std = out_KF_SD(:,11:13);
% acc_std_1 = out_KF_SD_1(:,11:13);
% acc_std_2 = out_KF_SD_2(:,11:13);
% plot(time, acc_std, 'LineWidth',1.2); hold on;
% plot(time, acc_std_1, '--', 'LineWidth',1.2);
% plot(time, acc_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% title('Accelerometer Bias STD');
% ylabel('STD (m/s^2)');
% legend('a_x (Baseline)', 'a_y (Baseline)', 'a_z (Baseline)', ...
%        'a_x (Type)', 'a_y (Type)', 'a_z (Type)', ...
%        'a_x (Uncert)', 'a_y (Uncert)', 'a_z (Uncert)', ...
%        'Location', 'best');
% 
% % 陀螺零偏 STD (°/h)
% subplot(2,1,2);
% gyro_std = out_KF_SD(:,14:16) * 180/pi * 3600;  % 转为 °/h
% gyro_std_1 = out_KF_SD_1(:,14:16) * 180/pi * 3600;
% gyro_std_2 = out_KF_SD_2(:,14:16) * 180/pi * 3600;
% plot(time, gyro_std, 'LineWidth',1.2); hold on;
% plot(time, gyro_std_1, '--', 'LineWidth',1.2);
% plot(time, gyro_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% title('Gyroscope Bias STD');
% ylabel('STD (°/h)');
% legend('\omega_x (Baseline)', '\omega_y (Baseline)', '\omega_z (Baseline)', ...
%        '\omega_x (Type)', '\omega_y (Type)', '\omega_z (Type)', ...
%        '\omega_x (Uncert)', '\omega_y (Uncert)', '\omega_z (Uncert)', ...
%        'Location', 'best', 'Interpreter', 'latex');
% 
% % 时钟不确定性（钟差 + 钟漂）
% figure('Position', [100, 100, 1000, 600]);
% clock_std = out_KF_SD(:,17:18);
% clock_std_1 = out_KF_SD_1(:,17:18);
% clock_std_2 = out_KF_SD_2(:,17:18);
% 
% plot(time, clock_std, 'LineWidth',1.2); hold on;
% plot(time, clock_std_1, '--', 'LineWidth',1.2);
% plot(time, clock_std_2, '-.', 'LineWidth',1.2);
% hold off; grid on;
% 
% title('Clock Error STD Comparison');
% xlabel('Time (s)'); ylabel('STD (m, m/s)');
% legend('Offset (Baseline)', 'Drift (Baseline)', ...
%        'Offset (Type)', 'Drift (Type)', ...
%        'Offset (Uncert)', 'Drift (Uncert)', 'Location', 'best');


% --- 确保 time 长度与 TCKF_result 一致 ---
n_result = size(TCKF_result, 1);

% 方法1：假设每历元1秒，从0开始
time = (0:n_result-1)';  % 列向量

% 方法2：如果你有外部时间向量（如 GNSS 时间），请对齐长度
% time = your_time_vector(1:n_result);  % 截取前n_result个

% --- 提取速度误差（假设第5~7列是速度误差）---
v_base   = TCKF_result(:, 8:10);
v_type   = TCKF_result_type(:, 8:10);
v_uncert = TCKF_result_uncert(:, 8:10);

% --- 检查尺寸 ---
if ~isequal(size(v_base,1), length(time))
    error('时间向量与数据行数不匹配：time(%d), v_base(%d)', length(time), size(v_base,1));
end

% --- 绘图 ---
figure('Position', [100, 100, 1000, 600]);
plot(time, v_base,     'LineWidth', 1.2);      hold on;
plot(time, v_type,     '--', 'LineWidth', 1.2);
plot(time, v_uncert,   '-.', 'LineWidth', 1.2);
hold off;
grid on;
title('Velocity Error Comparison (East, North, Up)');
xlabel('Time (s)');
ylabel('Velocity Error (m/s)');
legend('East (Baseline)', 'North (Baseline)', 'yaw (Baseline)', ...
       'East (Type)', 'North (Type)', 'yaw (Type)', ...
       'East (Uncert)', 'North (Uncert)', 'yaw (Uncert)', ...
       'Location', 'best');