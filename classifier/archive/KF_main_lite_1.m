% ------------------------------------------------------------------------
% Main code for LCKF
% modified by Paul's example
%
% by GH.Zhang 2019/04/26
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------
clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

% Constants
deg_to_rad = 0.01745329252;
rad_to_deg = 1/deg_to_rad;
micro_g_to_meters_per_second_squared = 9.80665E-6;

%% CONFIGURATION
% load('open3s.mat');
% 读取 GNSS 数据（CSV 格式）
GNSS_csv = readtable('GNSS_data.csv');
IMU_csv = readtable('IMU_data.csv');
IMU_data = table2array(IMU_csv);  % 转为普通 double 矩阵
CPT_csv = readtable('CPT_data.csv');
CPT_data = table2array(CPT_csv);  % 转为普通 double 矩阵
Stime = 0;
IMU_static_time = 30;
Etime = 150;

%% Starting Time
%比起更改代码，更改数据形式更方便，但是gnss data的读取只能更改

% 提取关键字段
GPS_Time = GNSS_csv.GPS_Time;          % 时间戳
% PRN = GNSS_csv.PRN;                    % 卫星编号
% Elevation = GNSS_csv.Elevation;        % 仰角
% Azimuth = GNSS_csv.Azimuth;            % 方位角
% CNR = GNSS_csv.CNR;                    % 载噪比
% Pr_Residual = GNSS_csv.Pr_Residual;    % 伪距残差
% Doppler = GNSS_csv.Doppler;            % 多普勒频移
% GT_Latitude = GNSS_csv.GT_Latitude;    % 真值纬度
% GT_Longitude = GNSS_csv.GT_Longitude;  % 真值经度
% GT_Height = GNSS_csv.GT_Height;        % 真值高度

uniqueTimes = unique(GPS_Time);
% 初始化元胞数组
GNSS_data = cell(length(uniqueTimes), 15); % 每行一个 epoch，最多 15 列

% 遍历每一行，构造元胞数组
for i = 1:size(uniqueTimes,1)
    % 构造卫星观测数据
    current_time = uniqueTimes(i);
    GNSS_data{i,1} = uniqueTimes(i);
    % 找出所有该时间点的行
    idx = GNSS_csv.GPS_Time == current_time;
    epoch_data = GNSS_csv(idx, :);  % 提取该 epoch 所有卫星数据
    
    GNSS_data{i,2} = [epoch_data.WLS_Latitude(1), epoch_data.WLS_Longitude(1), epoch_data.WLS_Height(1)];
    GNSS_data{i,3} = epoch_data.PRN; 
    GNSS_data{i,4} = epoch_data.Elevation; 
    GNSS_data{i,5} = epoch_data.Azimuth; 
    GNSS_data{i,9} = [epoch_data.CNR,epoch_data.Elevation]; 
    GNSS_data{i,10} = [epoch_data.PRN,epoch_data.corrected_Pr,epoch_data.SaX,epoch_data.SaY,epoch_data.SaZ];
    GNSS_data{i,11} = [epoch_data.Doppler,epoch_data.SaVX,epoch_data.SaVY,epoch_data.SaVZ,epoch_data.lambda];  
    GNSS_data{i,12} = [i];
    GNSS_data{i,13} = [epoch_data.VRx(1);epoch_data.VRy(1);epoch_data.VRz(1)];  
    GNSS_data{i,14} = [epoch_data.PRN,epoch_data.Pr_Residual];    
    GNSS_data{i,15} = [epoch_data.Type]; 
    GNSS_data{i,16} = [epoch_data.Uncertainty]; 

end
%对齐开始时间戳
% 在实际收数据中，这部分是人工完成的，尤其是cpt data和gnss data是人工标注的时间点（1，2，3.。。）
[~,temp.CPT_ids] = min(abs(CPT_data(:,1)-Stime));
CPT_data(1:temp.CPT_ids,:)=[];

[~,temp.GNSS_ids] = min(abs([GNSS_data{:,12}]-Stime));
GNSS_data(1:temp.GNSS_ids,:)=[];

[~,temp.IMU_ids] = min(abs(IMU_data(:,1)-Stime));
IMU_data(1:temp.IMU_ids,:)=[];
%对齐结束时间戳
CPT_data(CPT_data(:,1)>Etime,:)=[];
IMU_data(IMU_data(:,1)>Etime,:)=[];
GNSS_data([GNSS_data{:,12}]'>Etime,:)=[];
% 计算IMU零偏
temp.g_xyz = llh2xyz([GNSS_data{1,2}(1),GNSS_data{1,2}(2),GNSS_data{1,2}(3)].*[D2R,D2R,1]);
temp.g = Gravity_ECEF(temp.g_xyz');
[~,temp.CPT_id_st] = min(abs(CPT_data(:,1)-IMU_static_time));
[~,temp.IMU_id_st] = min(abs(IMU_data(:,1)-IMU_static_time));
IMU_bias = mean(IMU_data(1:temp.IMU_id_st,4:9));
% IMU_bias(1,1:3) = [0.2653-0.265,0.0742,0];
% IMU_bias(1,3) = 0;
% IMU_bias(1,1) = 0;
IMU_bias(1,3) = IMU_bias(1,3)-norm(temp.g);%z零偏观测值-g=z零偏

for idm = 1:1:size(IMU_data,1)
    IMU_data(idm,4:9) = IMU_data(idm,4:9)-IMU_bias;%清除零偏
%     IMU_data(idm,4) = 0;
end
% for idt = 1:1:size(GNSS_data,1)
%     temp.Lid = (GNSS_data{idt,6}'>32 & GNSS_data{idt,6}'<87);
%     GNSS_data{idt,6}(temp.Lid)=[];
%     GNSS_data{idt,7}(temp.Lid)=[];
%     GNSS_data{idt,8}(temp.Lid)=[];
% end

%% LCKF Setup
% Initial attitude uncertainty per axis (deg, converted to rad)
LC_KF_config.init_att_unc = deg2rad(2);
% Initial velocity uncertainty per axis (m/s)
LC_KF_config.init_vel_unc = 0.1;
% Initial position uncertainty per axis (m)
LC_KF_config.init_pos_unc = 10;
% Initial accelerometer bias uncertainty per instrument (micro-g, converted
% to m/s^2)
LC_KF_config.init_b_a_unc = 10000 * micro_g_to_meters_per_second_squared;
% Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
LC_KF_config.init_b_g_unc = 200 * deg_to_rad / 3600;

% Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
% LC_KF_config.gyro_noise_PSD = 0.01^2;
LC_KF_config.gyro_noise_PSD = 0.005^2;
% Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
% LC_KF_config.accel_noise_PSD = 0.2^2;
LC_KF_config.accel_noise_PSD = 0.008^2;
% Accelerometer bias random walk PSD (m^2 s^-5)
LC_KF_config.accel_bias_PSD = 1.0E-5;
% Gyro bias random walk PSD (rad^2 s^-3)
LC_KF_config.gyro_bias_PSD = 4.0E-11;

% % Position measurement noise SD per axis (m)
% LC_KF_config.pos_meas_SD = 25;
% % Velocity measurement noise SD per axis (m/s)
% LC_KF_config.vel_meas_SD = 10;

% Position measurement noise SD per axis (m)
LC_KF_config.pos_meas_SD = 25;
% Velocity measurement noise SD per axis (m/s)
LC_KF_config.vel_meas_SD = 10;

%% TCKF Setup
% Initial attitude uncertainty per axis (deg, converted to rad)
TC_KF_config.init_att_unc = deg2rad(2);
% Initial velocity uncertainty per axis (m/s)
TC_KF_config.init_vel_unc = 0.1*2;
% Initial position uncertainty per axis (m)
TC_KF_config.init_pos_unc = 10*2;
% Initial accelerometer bias uncertainty per instrument (micro-g, converted
% to m/s^2)
TC_KF_config.init_b_a_unc = 10000 * micro_g_to_meters_per_second_squared;
% Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
TC_KF_config.init_b_g_unc = 200 * deg_to_rad / 3600;
% Initial clock offset uncertainty per axis (m)
% TC_KF_config.init_clock_offset_unc = 20*10/10*10;%200
TC_KF_config.init_clock_offset_unc = 80;%200
% Initial clock drift uncertainty per axis (m/s)
% TC_KF_config.init_clock_drift_unc = 1*10*2/2;
TC_KF_config.init_clock_drift_unc = 40;

% Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
% TC_KF_config.gyro_noise_PSD = 0.01^2;
TC_KF_config.gyro_noise_PSD = 0.005^2;

% Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
% TC_KF_config.accel_noise_PSD = 0.2^2;
TC_KF_config.accel_noise_PSD = 0.008^2;

% Accelerometer bias random walk PSD (m^2 s^-5)
TC_KF_config.accel_bias_PSD = 1.0E-5;
% Gyro bias random walk PSD (rad^2 s^-3)
TC_KF_config.gyro_bias_PSD = 4.0E-11;
% Receiver clock frequency-drift PSD (m^2/s^3)
% TC_KF_config.clock_freq_PSD = 1*10/2;
TC_KF_config.clock_freq_PSD = 1;
% Receiver clock phase-drift PSD (m^2/s)
% TC_KF_config.clock_phase_PSD = 1*10/2;
TC_KF_config.clock_phase_PSD = 1;

% Pseudo-range measurement noise SD (m)
% TC_KF_config.pseudo_range_SD = 2.5*60/1.5;
TC_KF_config.pseudo_range_SD = 50;
% Pseudo-range rate measurement noise SD (m/s)
% TC_KF_config.range_rate_SD = 0.5*80/1.5;
TC_KF_config.range_rate_SD = 25;%20 before
TC_KF_config.bds_tag = 1;

%% Data Synchronization
temp.GT_t0 = CPT_data(2,1);
temp.GT_t = CPT_data(2:end,1);
temp.GNSS_t = [GNSS_data{:,12}]';
% temp.IMU_t = IMU_data(:,2)+IMU_data(:,3).*1e-9;
temp.IMU_t = IMU_data(:,1);
[~,temp.GNSS_ids] = min(abs(temp.GNSS_t-temp.GT_t0));
[~,temp.IMU_ids] = min(abs(temp.IMU_t-temp.GT_t0));
CPT_data(1,:) = [];
IMU_data(1:(temp.IMU_ids-1),:)=[];
GNSS_data(1:(temp.GNSS_ids-1),:)=[];

%check quaternion 四元数到欧拉角
eul = quat2eul(IMU_data(1,[10:13]));
eul = eul.*[R2D,R2D,R2D];


%% GNSS/INS Integration
int_type = 0;%0-origin//1-S-based//2-G-based
% Loosely coupled
[LCKF_result,~] = LCKF_INS_GNSS_GH(CPT_data,GNSS_data,IMU_data,LC_KF_config);
disp('LCKF Complete...');
% Tightly coupled
[TCKF_result,TCKF_result_type,TCKF_result_uncert,GNSS_llh,GNSS_llh_1,GNSS_llh_2] = TCKF_INS_GNSS_GH_1(CPT_data,GNSS_data,IMU_data,TC_KF_config);
disp('TCKF Complete...');

%% Error Evaluation
Error_2D = [];
Error_3D = [];
for idt = 1:1:size(CPT_data,1)
    CPT_xyz(idt,:) = llh2xyz(CPT_data(idt,4:6).*[D2R,D2R,1]);
    CPT_enu(idt,:) = xyz2enu(CPT_xyz(idt,:),CPT_xyz(1,:));
    id_GNSS_pos = find(round([GNSS_data{:,1}])'==CPT_data(idt,3));
    if ~isempty(id_GNSS_pos)
        id_GNSS_pos = id_GNSS_pos(1);
        temp.GNSS_llh = GNSS_llh(id_GNSS_pos,:).*[D2R,D2R,1];
    else
        temp.GNSS_llh = [nan,nan,nan];
    end
        temp.GNSS_xyz3 = llh2xyz([temp.GNSS_llh(1),temp.GNSS_llh(2),temp.GNSS_llh(3)]);
        temp.GNSS_xyz = llh2xyz([temp.GNSS_llh(1),temp.GNSS_llh(2),CPT_data(idt,6)]);
        temp.err2_GNSS = norm(CPT_xyz(idt,:)-temp.GNSS_xyz);
        temp.err3_GNSS = norm(CPT_xyz(idt,:)-temp.GNSS_xyz3);

        [~,id_LCKF_pos] = min(abs(LCKF_result(:,1)-CPT_data(idt,1)));
        temp.LCKF_xyz = llh2xyz(LCKF_result(id_LCKF_pos,2:4).*[D2R,D2R,1]);
        temp.LCKF_xyz2 = llh2xyz([LCKF_result(id_LCKF_pos,2:3),CPT_data(idt,6)].*[D2R,D2R,1]);
        temp.err2_LCKF = norm(CPT_xyz(idt,:)-temp.LCKF_xyz2);
        temp.err3_LCKF = norm(CPT_xyz(idt,:)-temp.LCKF_xyz);


        [~,id_TCKF_pos] = min(abs(TCKF_result(:,1)-CPT_data(idt,1)));
        temp.TCKF_xyz = llh2xyz(TCKF_result(id_TCKF_pos,2:4).*[D2R,D2R,1]);
        temp.TCKF_xyz2 = llh2xyz([TCKF_result(id_TCKF_pos,2:3),CPT_data(idt,6)].*[D2R,D2R,1]);
        temp.err2_TCKF = norm(CPT_xyz(idt,:)-temp.TCKF_xyz2);
        temp.err3_TCKF = norm(CPT_xyz(idt,:)-temp.TCKF_xyz);
       
        [~,id_TCKF_pos_1] = min(abs(TCKF_result_type(:,1)-CPT_data(idt,1)));
        temp.TCKF_xyz_1 = llh2xyz(TCKF_result_type(id_TCKF_pos_1,2:4).*[D2R,D2R,1]);
        temp.TCKF_xyz2_1 = llh2xyz([TCKF_result_type(id_TCKF_pos_1,2:3),CPT_data(idt,6)].*[D2R,D2R,1]);
        temp.err2_TCKF_1 = norm(CPT_xyz(idt,:)-temp.TCKF_xyz2_1);
        temp.err3_TCKF_1 = norm(CPT_xyz(idt,:)-temp.TCKF_xyz_1);

        [~,id_TCKF_pos_2] = min(abs(TCKF_result_uncert(:,1)-CPT_data(idt,1)));
        temp.TCKF_xyz_2 = llh2xyz(TCKF_result_uncert(id_TCKF_pos_2,2:4).*[D2R,D2R,1]);
        temp.TCKF_xyz2_2 = llh2xyz([TCKF_result_uncert(id_TCKF_pos_2,2:3),CPT_data(idt,6)].*[D2R,D2R,1]);
        temp.err2_TCKF_2 = norm(CPT_xyz(idt,:)-temp.TCKF_xyz2_2);
        temp.err3_TCKF_2 = norm(CPT_xyz(idt,:)-temp.TCKF_xyz_2);

        
        % Error_2D = [Error_2D;CPT_data(idt,1),temp.err2_GNSS,temp.err2_LCKF,temp.err2_TCKF];
        % Error_3D = [Error_3D;CPT_data(idt,1),temp.err3_GNSS,temp.err3_LCKF,temp.err3_TCKF];
        Error_2D = [Error_2D; ...
        CPT_data(idt,1), ...
        temp.err2_GNSS, ...
        temp.err2_LCKF, ...
        temp.err2_TCKF, ...
        temp.err2_TCKF_1, ...
        temp.err2_TCKF_2];
        
        Error_3D = [Error_3D; ...
        CPT_data(idt,1), ...
        temp.err3_GNSS, ...
        temp.err3_LCKF, ...
        temp.err3_TCKF, ...
        temp.err3_TCKF_1, ...
        temp.err3_TCKF_2];

end
% fprintf('GNSS RMSE: 2D->%5.1f(%5.1f) 3D->%5.1f(%5.1f)\n',...
%         rms(Error_2D(Error_2D(:,2)>0,2)),sum(Error_2D(:,2)>0)/size(Error_2D,1)*100,...
%         rms(Error_3D(Error_3D(:,2)>0,2)),sum(Error_3D(:,2)>0)/size(Error_3D,1)*100);
% fprintf('LCKF RMSE: 2D->%5.1f(%5.1fm) 3D->%5.1f(%5.1fm)\n',...
%         rms(Error_2D(Error_2D(:,3)>0,3)),rms(Error_2D(Error_2D(:,2)>0,3)),...
%         rms(Error_3D(Error_3D(:,3)>0,3)),rms(Error_3D(Error_3D(:,2)>0,3)));
% fprintf('TCKF RMSE: 2D->%5.1f(%5.1fm) 3D->%5.1f(%5.1fm)\n',...
%         rms(Error_2D(Error_2D(:,4)>0,4)),rms(Error_2D(Error_2D(:,2)>0,4)),...
%         rms(Error_3D(Error_3D(:,4)>0,4)),rms(Error_3D(Error_3D(:,2)>0,4)));
% 
% 
% 
% for idt = 1:1:size(GNSS_llh,1)
%     GNSS_xyz(idt,:) = llh2xyz(GNSS_llh(idt,:).*[D2R,D2R,1]);
%     temp.enu = xyz2enu(GNSS_xyz(idt,:),CPT_xyz(1,:));
%     GNSS_enu(idt,:) = [GNSS_data{idt,12},temp.enu'];
% end
% LCKF_result(LCKF_result(:,1)==0,:)=[];
% for idt = 1:1:size(LCKF_result,1)
%     LCKF_result_xyz(idt,:) = llh2xyz(LCKF_result(idt,2:4).*[D2R,D2R,1]);
%     LCKF_result_xyz_enu(idt,:) = xyz2enu(LCKF_result_xyz(idt,:),CPT_xyz(1,:));
% end
% TCKF_result(TCKF_result(:,1)==0,:)=[];
% for idt = 1:1:size(TCKF_result,1)
%     TCKF_result_xyz(idt,:) = llh2xyz(TCKF_result(idt,2:4).*[D2R,D2R,1]);
%     TCKF_result_xyz_enu(idt,:) = xyz2enu(TCKF_result_xyz(idt,:),CPT_xyz(1,:));
% end
% 
% 
% %% open plot
% figure(11);
% hold on;
% axis off;
% plot(CPT_data(:,5),CPT_data(:,4),'g-','LineWidth',5);
% plot(GNSS_llh(:,2),GNSS_llh(:,1),'o-','MarkerSize',5,...
%      'MarkerFaceColor',[0.95, 0.85, 0.1],'MarkerEdgeColor','k','LineWidth',0.5,'Color','k');
% plot(LCKF_result(LCKF_result(:,3)~=0,3),LCKF_result(LCKF_result(:,2)~=0,2),'o','MarkerSize',5,...
%      'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',0.5);
% plot(TCKF_result(:,3),TCKF_result(:,2),'o-','MarkerSize',5,...
%      'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',0.5,'Color','r');
% hBase = plot_openstreetmap('Alpha', 1, 'Scale', 50);
% plot(CPT_data(:,5),CPT_data(:,4),'g-','LineWidth',5);
% plot(GNSS_llh(:,2),GNSS_llh(:,1),'o-','MarkerSize',5,...
%      'MarkerFaceColor',[0.95, 0.85, 0.1],'MarkerEdgeColor','k','LineWidth',0.5,'Color','k');
% plot(LCKF_result(LCKF_result(:,3)~=0,3),LCKF_result(LCKF_result(:,2)~=0,2),'o','MarkerSize',5,...
%      'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',0.5);
% plot(TCKF_result(:,3),TCKF_result(:,2),'o-','MarkerSize',5,...
%      'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',0.5,'Color','r');
% legend('Truth','WLS','LCKF','TCKF');
% 
% 
% 
% 





%% 输出 RMSE 报告
fprintf('GNSS RMSE: 2D->%5.1f(m) 3D->%5.1f(m)\n',...
    rms(Error_2D(Error_2D(:,2)>0,2)), rms(Error_3D(Error_3D(:,2)>0,2)));

fprintf('LCKF RMSE: 2D->%5.1f(m) 3D->%5.1f(m)\n',...
    rms(Error_2D(Error_2D(:,3)>0,3)), rms(Error_3D(Error_3D(:,3)>0,3)));

fprintf('TCKF-Origin RMSE: 2D->%5.1f(m) 3D->%5.1f(m)\n',...
    rms(Error_2D(Error_2D(:,4)>0,4)), rms(Error_3D(Error_3D(:,4)>0,4)));

fprintf('TCKF-Type RMSE: 2D->%5.1f(m) 3D->%5.1f(m)\n',...
    rms(Error_2D(Error_2D(:,5)>0,5)), rms(Error_3D(Error_3D(:,5)>0,5)));

fprintf('TCKF-Uncert RMSE: 2D->%5.1f(m) 3D->%5.1f(m)\n',...
    rms(Error_2D(Error_2D(:,6)>0,6)), rms(Error_3D(Error_3D(:,6)>0,6)));

%% 转换所有轨迹到 ENU（用于绘图）
for idt = 1:size(GNSS_llh,1)
    GNSS_xyz(idt,:) = llh2xyz(GNSS_llh(idt,:).*[D2R,D2R,1]);
    temp.enu = xyz2enu(GNSS_xyz(idt,:), CPT_xyz(1,:));
    GNSS_enu(idt,:) = [GNSS_data{idt,12}, temp.enu'];
end

LCKF_result(LCKF_result(:,1)==0,:) = [];
for idt = 1:size(LCKF_result,1)
    LCKF_result_xyz(idt,:) = llh2xyz(LCKF_result(idt,2:4).*[D2R,D2R,1]);
    LCKF_result_xyz_enu(idt,:) = xyz2enu(LCKF_result_xyz(idt,:), CPT_xyz(1,:));
end

TCKF_result(TCKF_result(:,1)==0,:) = [];
TCKF_result_type(TCKF_result_type(:,1)==0,:) = [];
TCKF_result_uncert(TCKF_result_uncert(:,1)==0,:) = [];

for idt = 1:size(TCKF_result,1)
    TCKF_result_xyz(idt,:) = llh2xyz(TCKF_result(idt,2:4).*[D2R,D2R,1]);
    TCKF_result_xyz_enu(idt,:) = xyz2enu(TCKF_result_xyz(idt,:), CPT_xyz(1,:));
end

for idt = 1:size(TCKF_result_type,1)
    TCKF_result_xyz_type(idt,:) = llh2xyz(TCKF_result_type(idt,2:4).*[D2R,D2R,1]);
    TCKF_result_xyz_enu_type(idt,:) = xyz2enu(TCKF_result_xyz_type(idt,:), CPT_xyz(1,:));
end

for idt = 1:size(TCKF_result_uncert,1)
    TCKF_result_xyz_uncert(idt,:) = llh2xyz(TCKF_result_uncert(idt,2:4).*[D2R,D2R,1]);
    TCKF_result_xyz_enu_uncert(idt,:) = xyz2enu(TCKF_result_xyz_uncert(idt,:), CPT_xyz(1,:));
end

%% 绘图
figure(11); clf;
hold on; axis off;

% 真值
plot(CPT_data(:,5), CPT_data(:,4), 'g-', 'LineWidth', 5);

% GNSS WLS
plot(GNSS_llh(:,2), GNSS_llh(:,1), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor',[0.95, 0.85, 0.1], 'MarkerEdgeColor','k', 'LineWidth',0.5, 'Color','k');

% LCKF
plot(LCKF_result(LCKF_result(:,3)~=0,3), LCKF_result(LCKF_result(:,2)~=0,2), ...
    'o', 'MarkerSize', 5, 'MarkerFaceColor','b', 'MarkerEdgeColor','k', 'LineWidth',0.5);

% TCKF Origin
plot(TCKF_result(:,3), TCKF_result(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'LineWidth',0.5, 'Color','r');

% TCKF Type
plot(TCKF_result_type(:,3), TCKF_result_type(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','m', 'MarkerEdgeColor','m', 'LineWidth',0.5, 'Color','m');

% TCKF Uncertainty
plot(TCKF_result_uncert(:,3), TCKF_result_uncert(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','c', 'MarkerEdgeColor','c', 'LineWidth',0.5, 'Color','c');

% 底图
hBase = plot_openstreetmap('Alpha', 1, 'Scale', 50);

% 重画一次（避免被底图覆盖）
plot(CPT_data(:,5), CPT_data(:,4), 'g-', 'LineWidth', 5);
plot(GNSS_llh(:,2), GNSS_llh(:,1), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor',[0.95, 0.85, 0.1], 'MarkerEdgeColor','k', 'LineWidth',0.5, 'Color','k');
plot(LCKF_result(LCKF_result(:,3)~=0,3), LCKF_result(LCKF_result(:,2)~=0,2), ...
    'o', 'MarkerSize', 5, 'MarkerFaceColor','b', 'MarkerEdgeColor','k', 'LineWidth',0.5);
plot(TCKF_result(:,3), TCKF_result(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'LineWidth',0.5, 'Color','r');
plot(TCKF_result_type(:,3), TCKF_result_type(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','m', 'MarkerEdgeColor','m', 'LineWidth',0.5, 'Color','m');
plot(TCKF_result_uncert(:,3), TCKF_result_uncert(:,2), 'o-', 'MarkerSize', 5, ...
    'MarkerFaceColor','c', 'MarkerEdgeColor','c', 'LineWidth',0.5, 'Color','c');

legend('Truth', 'WLS', 'LCKF', 'TCKF-Origin', 'TCKF-Type', 'TCKF-Uncert');
title('Trajectory Comparison of Different Weighting Strategies');









