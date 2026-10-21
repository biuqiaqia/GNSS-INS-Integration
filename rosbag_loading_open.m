% ------------------------------------------------------------------------
% Main code for rosbag loading
%
% by GH.Zhang 2019/04/26
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------
clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

% file_name = 'middle_urban';
% rosbag_path = ['data\rosbag\',file_name,'.bag'];
% filename_nav = 'data\eph\hksc090j.19n';  
% GNSS_type = 2;

% file_name = 'middle_urban2';
% rosbag_path = ['data\rosbag\',file_name,'.bag'];
% filename_nav = 'data\eph\hksc118m.19n';  
% GNSS_type = 1;
% GT_type = 1;

% file_name = 'test1';
% rosbag_path = ['data\rosbag\',file_name,'.bag'];
% filename_nav = 'data\eph\hksc340q.19n';  
% GNSS_type = 1;
% GT_type = 2;
% gt_waypoint = [22.306655,114.187038,20,1;...
%                22.30665501,114.18703801,20.01,0];%LLH and epoch

% Whampoa open-sky 1
file_name = 'open1';
rosbag_path = ['data\rosbag\',file_name,'.bag'];
filename_nav = 'data\eph\hksc341g.19n';  
refsat_obs = 'data\eph\hksc341g.19o';  
GNSS_type = 1;
GT_type = 2;
gt_waypoint = csvread(['data\rosbag\',file_name,'_gt.csv']);
gt_waypoint(2:end,4) = gt_waypoint(2:end,4)+7-1;
gt_heading = 193.36;

% % Whampoa open-sky 2
% file_name = 'open2';
% rosbag_path = ['data\rosbag\',file_name,'.bag'];
% filename_nav = 'data\eph\hksc341g.19n';  
% refsat_obs = 'data\eph\hksc341g.19o';  
% GNSS_type = 1;
% GT_type = 2;
% gt_waypoint = csvread(['data\rosbag\',file_name,'_gt.csv']);
% gt_waypoint(2:end,4) = gt_waypoint(2:end,4)+2;
% gt_heading = 193.36;

% % Whampoa open-sky 3
% file_name = 'open3';
% rosbag_path = ['data\rosbag\',file_name,'.bag'];
% filename_nav = 'data\eph\hksc341g.19n';  
% refsat_obs = 'data\eph\hksc341g.19o';  
% GNSS_type = 1;
% GT_type = 2;
% gt_waypoint = csvread(['data\rosbag\',file_name,'_gt.csv']);
% gt_waypoint(2:end,4) = gt_waypoint(2:end,4)+2;
% gt_heading = 193.36;

bag = rosbag(rosbag_path);

%% load ground-truth
if GT_type == 1
    disp('GT loading from rosbag...');
    bagselect2 = select(bag,'Topic','/novatel_data/bestpos');
    msgStructs2 = readMessages(bagselect2,'DataFormat','struct');
    bagselect2b = select(bag,'Topic','/novatel_data/inspvax');
    msgStructs2b = readMessages(bagselect2b,'DataFormat','struct');
    bagselect2c = select(bag,'Topic','/novatel_data/corrimudata');
    msgStructs2c = readMessages(bagselect2c,'DataFormat','struct');
    for idt = 1:1:size(bagselect2.MessageList,1)
        CPT_data(idt,1) = bagselect2.MessageList.Time(idt);%ROS Time
        CPT_data(idt,2) = double(msgStructs2{idt,1}.Header.GpsWeek);
        CPT_data(idt,3) = double(msgStructs2{idt,1}.Header.GpsWeekSeconds)/1000;
        CPT_data(idt,4) = double(msgStructs2{idt,1}.Latitude);
        CPT_data(idt,5) = double(msgStructs2{idt,1}.Longitude);
        CPT_data(idt,6) = double(msgStructs2{idt,1}.Altitude);
    end

    dt = double(msgStructs2b{1,1}.Header.GpsWeekSeconds)/1000-double(msgStructs2{1,1}.Header.GpsWeekSeconds)/1000;
    for idt = 1:1:size(bagselect2b.MessageList,1)
        if idt+dt>0
        CPT_data(idt+dt,7) = double(msgStructs2b{idt,1}.NorthVelocity);
        CPT_data(idt+dt,8) = double(msgStructs2b{idt,1}.EastVelocity);
        CPT_data(idt+dt,9) = -double(msgStructs2b{idt,1}.UpVelocity);%down
        CPT_data(idt+dt,10) = double(msgStructs2b{idt,1}.Roll);
        CPT_data(idt+dt,11) = double(msgStructs2b{idt,1}.Pitch);
        CPT_data(idt+dt,12) = double(msgStructs2b{idt,1}.Azimuth);
        end
    end

    dt = double(msgStructs2c{1,1}.Header.GpsWeekSeconds)/1000-double(msgStructs2{1,1}.Header.GpsWeekSeconds)/1000;
    for idt = 1:1:size(bagselect2c.MessageList,1)
        CPT_data(idt+dt,13) = double(msgStructs2c{idt,1}.XAccel);
        CPT_data(idt+dt,14) = double(msgStructs2c{idt,1}.YAccel);
        CPT_data(idt+dt,15) = double(msgStructs2c{idt,1}.ZAccel);
        CPT_data(idt+dt,16) = double(msgStructs2c{idt,1}.PitchRate);
        CPT_data(idt+dt,17) = double(msgStructs2c{idt,1}.RollRate);
        CPT_data(idt+dt,18) = double(msgStructs2c{idt,1}.YawRate);
    end
    user_pos_xyz0 = llh2xyz(CPT_data(1,4:6).*[D2R,D2R,1]);

elseif GT_type == 2
    user_pos_xyz0 = llh2xyz(gt_waypoint(1,1:3).*[D2R,D2R,1]);

end

%% load GNSS
% read from eph
global weights 
weights = 0;
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';
%enabled constellations
GPS_flag = 1;  GLO_flag = 1;  GAL_flag = 1;  BDS_flag = 1;  QZS_flag = 1; SBS_flag = 1;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';
[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;
goGNSS;
[Eph, iono] = load_RINEX_nav(filename_nav);            
fprintf('--> RINEX files loading complete\n');

% read from rosbag
disp('GNSS loading from rosbag...');
if GNSS_type==1
    bagselect1 = select(bag,'Topic','/ublox_gps_m8t/rxmraw');
    msgStructs = readMessages(bagselect1,'DataFormat','struct');
elseif GNSS_type == 2
    bagselect1 = select(bag,'Topic','/GNSS_');
    msgStructs = readMessages(bagselect1,'DataFormat','struct');
end

idx_run_data = 0;
XR0 = [];
el_mask = 0;

for t = 1:1:size(bagselect1.MessageList,1)
    if GNSS_type == 1
        [time_GPS(t,1)] = weektow2time(double(msgStructs{t,1}.Week), double(msgStructs{t,1}.RcvTOW), 'G');
        [~, all_sow(t,1)] = time2weektow(time_GPS(t,1));
        temp.gnssid = double([msgStructs{t,1}.Meas.GnssId]);
        temp.svid = double([msgStructs{t,1}.Meas.SvId])+(temp.gnssid==3).*86;%3-Beidou
%         temp.svid = double([msgStructs{t,1}.Meas.SvId])+(temp.gnssid==6).*32;%6-GLONASS
        temp.prn = temp.svid';
        for idm = 1:1:size(temp.prn,1)
            pr1_R(temp.prn(idm),t) = double(msgStructs{t,1}.Meas(idm).PrMes);
            ph1_R(temp.prn(idm),t) = double(msgStructs{t,1}.Meas(idm).CpMes);
            dop1_R(temp.prn(idm),t) = double(msgStructs{t,1}.Meas(idm).DoMes);
            snr1_R(temp.prn(idm),t) = double(msgStructs{t,1}.Meas(idm).Cno);
        end
    elseif GNSS_type == 2
        [time_GPS(t,1)] = weektow2time(CPT_data(1,2), double(msgStructs{t,1}.GNSSRaws(1).GNSSTime)/1000, 'G');
        [~, all_sow(t,1)] = time2weektow(time_GPS(t,1));
        temp.prn = [msgStructs{t,1}.GNSSRaws.PrnSatellitesIndex]';
        for idm = 1:1:size(temp.prn,1)
            pr1_R(temp.prn(idm),t) = double(msgStructs{t,1}.GNSSRaws(idm).Pseudorange)-...
                                     double(msgStructs{t,1}.GNSSRaws(idm).SatClkErr)+...
                                     double(msgStructs{t,1}.GNSSRaws(idm).ErrTropo)+...
                                     double(msgStructs{t,1}.GNSSRaws(idm).ErrIono);
            ph1_R(temp.prn(idm),t) = double(msgStructs{t,1}.GNSSRaws(idm).Visable);
            dop1_R(temp.prn(idm),t) = double(msgStructs{t,1}.GNSSRaws(idm).Visable);
            snr1_R(temp.prn(idm),t) = double(msgStructs{t,1}.GNSSRaws(idm).Snr);
        end
    end
    % Gogps calculation
    time_rx = time_GPS(t);       
    time_rx_ = mod(time_rx,604800);
    idx_run_data = idx_run_data + 1;    
   
    pr1 = pr1_R(:,t);
    dop1 = dop1_R(:,t);  
    snr = snr1_R(:,t);
    nSatTot = size(pr1,1);
    Eph_t = rt_find_eph(Eph, time_GPS(t), nSatTot);      
    lambda = goGNSS.getGNSSWavelengths(Eph_t, nSatTot);
    eph_avail = Eph(30,:);
    dtR = 0;
    % Rough Satellite Position Estimation
        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, 20200*1e3*ones(1,length(eph_avail)), eph_avail, Eph, [], [], zeros(1,length(eph_avail)), zeros(1,length(eph_avail)), dtR);
    
    sat = find(pr1 ~= 0);
    eph_avail = Eph(30,:);
    sat = sat(ismember(sat, eph_avail));
    
    data_all{t,1} = time_rx_; % GPS TOW; 
    GPST(t,1) = time_rx_; % GPS TOW;

    for jdx = 1 : length(XS)
        target_pos_xyz = [XS(jdx,1) XS(jdx,2) XS(jdx,3)];    
        target_pos_enu = xyz2enu(target_pos_xyz,user_pos_xyz0);
        azimuth(jdx) = ((pi/2)-atan2(target_pos_enu(1),target_pos_enu(2)))*R2D;
        if (azimuth(jdx)<0)
            azimuth(jdx) = azimuth(jdx) + 360 ;
        end
        elevation(jdx) = ( R2D*atan(target_pos_enu(3)/norm(target_pos_enu(1:2))));
    end
    
    azimuth_ned = azimuth.*(-1)+90;
    azimuth_ned(azimuth_ned<0) = azimuth_ned(azimuth_ned<0)+360;
    data_all{t,6} = eph_avail; % ALL PRN;
    data_all{t,7} = elevation; % elevation angle for each grid location;
    data_all{t,8} = azimuth_ned; % azimuth angle for each grid location;            
    % ols to estimate rough receiver position
    sat = find(pr1 ~= 0);
    eph_avail = Eph(30,:);
    sat = sat(ismember(sat, eph_avail));

    flag_gps_enough = 1;
%     if (length(find(sat <= 32))<4)&&(length(find(sat >= 87))<4) 
%     if (length(find(sat <= 32))+length(find(sat >= 87))<5)
%         flag_gps_enough = 0;  
%     end
    min_nsat = 4;
    
    XR_LLH(idx_run_data,:) = zeros(1,3);
    VR_LLH(idx_run_data,:) = zeros(1,3);
    Time_WLS(idx_run_data) = all_sow(t);
    
    if (size(sat,1) >= min_nsat && flag_gps_enough)
        lambda = lambda(sat);
        pseudorange = pr1(sat);
        doppler = dop1(sat);
        snr = snr(sat);

        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR);

        num_sys  = length(unique(sys(sys ~= 0)));
        min_nsat = 3 + num_sys;

        if (isempty(XR0))
            XR0 = zeros(3,1);
        end

        index = find(no_eph == 0);   
        nsat_avail = length(index);

        if (nsat_avail < min_nsat) %if available observations are not enough, return empty variables
            fprintf('not enough sv\n');
%             XR = [nan;nan;nan];
            [az, el, dist_] = topocent(XR, XS(index,:));  
            outlier_prn{t,1} = [];
            err_tropo = tropo_error_correction(el, hR);
            err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, iono, []);
            [XS_tt, dtS_tt, XS_tx_tt, VS_tx_tt, time_tx_tt, no_eph, sys] = satellite_positions(time_rx+tt, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));        
            VS = (XS_tt - XS) / tt;
            VR = [nan;nan;nan];
        else
            %iterative least-squares from XR0,i.e. given coordinates or the center of the Earth (i.e. [0; 0; 0])
            n_iter_max = 5;
            n_iter = 0;
            var_SPP(1) = Inf;
            SPP_threshold=4; %meters 
            outlier_prn{t,1} = [];
            while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
                [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), snr(index), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS(index), zeros(nsat_avail,1), zeros(nsat_avail,1), sys(index), SPP_threshold);       
                outlier_prn{t,1} = bad_obs;
                XR0 = XR;
                n_iter = n_iter + 1;
            end        
            %satellite topocentric coordinates (azimuth, elevation, distance)
            [az, el, dist_] = topocent(XR, XS(index,:));        
            %cartesian to geodetic conversion of ROVER coordinates
            [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
            %radians to degrees
            phiR = phiR * 180 / pi;
            lamR = lamR * 180 / pi;
            %computation of tropospheric errors
            err_tropo = tropo_error_correction(el, hR);
            %computation of ionospheric errors
            err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, iono, []);  

            %accurate Clock estimation
            [dtR, var_dtR, bad_obs, bad_epoch, var_SPP, ~, is_bias] = LS_SA_code_clock_Qmo(pseudorange(index), snr(index), el, dist_, dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);
            %accurate SV Pos
            [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));

            %SV velocity estimation (learned from RTKLIB)
            tt = 1e-3;
            [XS_tt, dtS_tt, XS_tx_tt, VS_tx_tt, time_tx_tt, no_eph, sys] = satellite_positions(time_rx+tt, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));        
            VS = (XS_tt - XS) / tt;
            dtSV = (dtS_tt - dtS) / tt;

            %aacurate rover pos estimation
            n_iter = 0;
            while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
                [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), snr(index), el, zeros(nsat_avail,1), dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);       
                outlier_prn{t,1} = bad_obs;
                XR0 = XR;
                n_iter = n_iter + 1;
            end
            %satellite topocentric coordinates (azimuth, elevation, distance)
            [az, el, dist_] = topocent(XR, XS(index,:));   
            az_enu = az.*(-1)+90;
            az_enu(az_enu<0,1) = az_enu(az_enu<0,1)+360;
            
            %aacurate rover vel estimation
           [VR, dtRV, pr_rate_lt,A,LOS_vector_lt] = LS_SA_code_Vel_Qmo_debug(XR, XS(index,:), VS(index,:), doppler(index), snr(index), el, lambda(index), sys(index));       
        end 
    end
    rcr_loc_temp = xyz2llh(XR);
    data_all{t,2} = [rcr_loc_temp(2)*R2D rcr_loc_temp(1)*R2D rcr_loc_temp(3)]; % location that receiver estimated (LLA)
    data_all{t,3} = sat(index,1); % PRN list that receiver received 
    data_all{t,4} = el; % elevation angle that receiver received
    data_all{t,5} = az; % azimuth angle that receiver received
    data_all{t,9} = []; % Measurement Features
    data_all{t,10} = [sat(index,1),pseudorange(index,1)+dtS(index,1).*goGNSS.V_LIGHT-err_tropo-err_iono,XS(index,:)];% Measurement data
    data_all{t,11} = [doppler(index,1),VS(index,:),lambda(index,1)];% doppler shift
    data_all{t,13} = VR;
    % Available SV & LOS Ratio
    nSV = size(data_all{t,3},1);
    [~,ia,~] = unique(data_all{t,6}); % negelct the duplicated ephs
    idx_eph_visible_sat = find(data_all{t,7}(ia)>el_mask);
    temp_prn_eph = data_all{t,6}(ia(idx_eph_visible_sat));
    data_all{t,6} = data_all{t,6}(ia(idx_eph_visible_sat));
    data_all{t,7} = data_all{t,7}(ia(idx_eph_visible_sat));
    data_all{t,8} = data_all{t,8}(ia(idx_eph_visible_sat));
    % Normalized Pseudorange Residual
    [~,pr_residual{t,1},~] = Least_square_positioning(data_all{t,10}(:,1:5),user_pos_xyz0);
    pr_residual{t,2} = pr_residual{t,1}./mean(abs(pr_residual{t,1}));
    % Pseudorange Rate Consistency
    if t == 1
        Pr_rate{1,1} = [ones(size(data_all{1,11},1),1).*9999,zeros(size(data_all{1,11},1),1)];
        for idm = 1:1:size(data_all{t,11},1)
            Pr_rate{t,1}(idm,1) = data_all{t,11}(idm,1)*data_all{t,11}(idm,5);
            Pr_rate{t,1}(idm,2) = Pr_rate{t,1}(idm,1)-9999;
        end
    else
        for idm = 1:1:size(data_all{t,11},1)
            if ~isempty(data_all{t-1,3})&&any(data_all{t,3}(idm,1)==data_all{t-1,3}(:,1))
                Pr_rate{t,1}(idm,1) = data_all{t,11}(idm,1)*data_all{t,11}(idm,5);            
                id_pr_last = find(data_all{t,10}(idm,1)==data_all{t-1,10}(:,1));
                Pr_rate{t,1}(idm,2) = -(data_all{t,10}(idm,2) - data_all{t-1,10}(id_pr_last,2));
            else
                Pr_rate{t,1}(idm,1) = 9999;
                Pr_rate{t,1}(idm,2) = 0;
            end
        end
    end
    Pr_rate{t,1}(:,3) = Pr_rate{t,1}(:,1)-Pr_rate{t,1}(:,2);
    % Measurement Features
    data_all{t,9} = [data_all{t,3},...                                      % 1 - PRN
                     ones(nSV,1).*nSV,...                                   % 2 - nSV
                     ones(nSV,1).*(nSV/size(data_all{t,6},2)),...           % 3 - LOS ratio
                     snr(index),...                                         % 4 - C/N0
                     el,...                                                 % 5 - Elevation
                     pr_residual{t,2},...                                   % 6 - Normalized Pr_residual
                     Pr_rate{t,1}(:,3),...                                  % 7 - Pr_rate consistency
                     Pr_rate{t,1}(:,1)];                                    % 8 - Pr_rate
	% Outlier exclusion by gogps
    if ~isempty(outlier_prn{t,1})&&1
        data_all{t,3}(outlier_prn{t,1},:)=[];
        data_all{t,4}(outlier_prn{t,1},:)=[];
        data_all{t,5}(outlier_prn{t,1},:)=[];
        data_all{t,9}(outlier_prn{t,1},:)=[];
        data_all{t,10}(outlier_prn{t,1},:)=[];
        data_all{t,11}(outlier_prn{t,1},:)=[];
    end
    % Exclusion by elevation mask
    mask_el_id = data_all{t,4}<=el_mask;
    data_all{t,3}(mask_el_id,:)=[];
    data_all{t,4}(mask_el_id,:)=[];
    data_all{t,5}(mask_el_id,:)=[];
    data_all{t,9}(mask_el_id,:)=[];
    data_all{t,10}(mask_el_id,:)=[];
    data_all{t,11}(mask_el_id,:)=[];
    data_all{t,12} = bagselect1.MessageList.Time(t);%ROS Time
end


% GNSS_data = data_all;

%% load IMU
disp('IMU loading from rosbag...');
bagselect3 = select(bag,'Topic','/imu/data');
msgStructs3 = readMessages(bagselect3,'DataFormat','struct');
for idt = 1:1:size(bagselect3.MessageList,1)
    IMU_data(idt,1) = bagselect3.MessageList.Time(idt);%ROS Time
    IMU_data(idt,2) = double(msgStructs3{idt,1}.Header.Stamp.Sec);
    IMU_data(idt,3) = double(msgStructs3{idt,1}.Header.Stamp.Nsec);
    IMU_data(idt,4) = double(msgStructs3{idt,1}.LinearAcceleration.X);%m/s^2
    IMU_data(idt,5) = double(msgStructs3{idt,1}.LinearAcceleration.Y);
    IMU_data(idt,6) = double(msgStructs3{idt,1}.LinearAcceleration.Z);
    IMU_data(idt,7) = double(msgStructs3{idt,1}.AngularVelocity.X);%rad/sec
    IMU_data(idt,8) = double(msgStructs3{idt,1}.AngularVelocity.Y);
    IMU_data(idt,9) = double(msgStructs3{idt,1}.AngularVelocity.Z);
    IMU_data(idt,10) = double(msgStructs3{idt,1}.Orientation.X);%quaternion
    IMU_data(idt,11) = double(msgStructs3{idt,1}.Orientation.Y);
    IMU_data(idt,12) = double(msgStructs3{idt,1}.Orientation.Z);
    IMU_data(idt,13) = double(msgStructs3{idt,1}.Orientation.W);
end

%% GT_label
if GT_type == 2
    [num_data_epoch, ~] = size(data_all);
    % gt_waypoint
    gt_lon_epoch = [];
    gt_lat_epoch = [];
    gt_alt_epoch = [];
    gt_waypoint(end,4)=t;
    for idwp = 1:1:size(gt_waypoint,1)-1
        gt_lon_epoch = [gt_lon_epoch,gt_waypoint(idwp,2):(gt_waypoint(idwp+1,2)-gt_waypoint(idwp,2))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,2)];
        gt_lat_epoch = [gt_lat_epoch,gt_waypoint(idwp,1):(gt_waypoint(idwp+1,1)-gt_waypoint(idwp,1))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,1)];
        gt_alt_epoch = [gt_alt_epoch,gt_waypoint(idwp,3):(gt_waypoint(idwp+1,3)-gt_waypoint(idwp,3))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,3)];
        gt_lon_epoch(end) = [];
        gt_lat_epoch(end) = [];
        gt_alt_epoch(end) = [];
    end
    gt_lon_epoch = [gt_lon_epoch,gt_waypoint(end,2)];
    gt_lat_epoch = [gt_lat_epoch,gt_waypoint(end,1)];
    gt_alt_epoch = [gt_alt_epoch,gt_waypoint(end,3)];
    
    %pseudo-CPT
    GPST_min = round(GPST(1,1));
    GPST_max = round(GPST(end,1));
    for idt = 1:1:(GPST_max-GPST_min+1)
        [~,idd] = min(abs(GPST(:,1)-(GPST_min-1+idt)));
        CPT_data(idt,1) = data_all{idd,12};%ROS Time
        CPT_data(idt,2) = double(msgStructs{1,1}.Week);
        CPT_data(idt,3) = GPST_min-1+idt;
        CPT_data(idt,4) = gt_lat_epoch(idd);
        CPT_data(idt,5) = gt_lon_epoch(idd);
        CPT_data(idt,6) = gt_alt_epoch(idd);
        CPT_data(idt,12) = gt_heading;
    end
end

%% Shift ROS time bias
ros_time0 = IMU_data(1,1)-1;
for idt = 1:1:size(data_all,1)
    data_all{idt,12} = data_all{idt,12}-ros_time0;
end
IMU_data(:,1) = IMU_data(:,1)-ros_time0;
IMU_data(:,2) = IMU_data(:,2)-ros_time0;
CPT_data(:,1) = CPT_data(:,1)-ros_time0;

%down to 1Hz sampling GNSS
for idm = 1:1:size(CPT_data,1)
    temp.GNSS_t = [data_all{:,12}]';
    [~,id_GNSS_sample] = min(abs(temp.GNSS_t-CPT_data(idm,1)));
    GNSS_data(idm,:) = data_all(id_GNSS_sample,:);
    GNSS_data{idm,1} = round(GNSS_data{idm,1});
end

for idt = 1:1:size(CPT_data,1)
    CPT_xyz(idt,:) = llh2xyz(CPT_data(idt,4:6).*[D2R,D2R,1]);
    CPT_enu(idt,:) = xyz2enu(CPT_xyz(idt,:),CPT_xyz(1,:));
end

%% Pseudorange error
ref_xyz = [-2414266.9197,5386768.9868,2407460.0314];%hksc


[Pr_error] = DD_pr_error_from_Ref(filename_nav,refsat_obs,ref_xyz,GNSS_data(:,[1,10,4]),CPT_xyz);
for idm = 1:1:size(GNSS_data,1)
    if ~isempty(Pr_error{idm,1})
        GNSS_data{idm,14} = Pr_error{idm,1};%DD_Pr_Error
        GNSS_data{idm,15} = GNSS_data{idm,10};
        if size(GNSS_data{idm,15},1) == size(GNSS_data{idm,14},1)
            GNSS_data{idm,15}(:,2) = GNSS_data{idm,15}(:,2)-GNSS_data{idm,14}(:,2);% Corrected GNSS by DD
        else
            disp(['Refsat missing at Epoch: ',num2str(idm)]);
            temp.match_id = ismember(GNSS_data{idm,15}(:,1),GNSS_data{idm,14}(:,1));
            GNSS_data{idm,15}(temp.match_id,2) = GNSS_data{idm,15}(temp.match_id,2)-GNSS_data{idm,14}(:,2);% Corrected GNSS by DD
        end
    else
        GNSS_data{idm,15} = GNSS_data{idm,10};
    end
end

%% Check data correctness
for idt = 1:1:size(GNSS_data,1)
    [ls_xyz(idt,:),~,~] = Least_square_positioning(GNSS_data{idt,10}(:,1:5),user_pos_xyz0);
    ls_llh(idt,:) = xyz2llh(ls_xyz(idt,:));
    ls_llh(idt,:) = ls_llh(idt,:).*[R2D,R2D,1];
    [cls_xyz(idt,:),~,~] = Least_square_positioning(GNSS_data{idt,15}(:,1:5),user_pos_xyz0);
    cls_llh(idt,:) = xyz2llh(cls_xyz(idt,:));
    cls_llh(idt,:) = cls_llh(idt,:).*[R2D,R2D,1];
end
% CPT_data(35,4:6) = CPT_data(34,4:6);

figure(1);
hold on;
plot(ls_llh(:,2),ls_llh(:,1),'o','MarkerSize',10,...
     'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
plot(cls_llh(:,2),cls_llh(:,1),'o','MarkerSize',10,...
     'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',2);
plot(CPT_data(:,5),CPT_data(:,4),'g-','LineWidth',3);
hBase = plot_openstreetmap('Alpha', 0.3, 'Scale', 50);

hBase = plot_openstreetmap('Alpha', 1, 'Scale', 50);
plot(ls_llh(1:100,2),ls_llh(1:100,1),'o','MarkerSize',10,...
     'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
plot(cls_llh(1:100,2),cls_llh(1:100,1),'o','MarkerSize',10,...
     'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',2);
plot(CPT_data(1:100,5),CPT_data(1:100,4),'g-','LineWidth',3);

% hSea  = plot_openstreetmap('Alpha', 0.5, 'Scale', 2, 'BaseUrl', "http://tiles.openseamap.org/seamark");
% plot_google_map('maptype', 'roadmap');

save(['data\mat\',file_name,'m1.mat'],'CPT_data','GNSS_data','IMU_data');

% 
% figure(2)
% subplot(4,1,1);
% title('Position')
% hold on;
% plot(CPT_data(:,1),CPT_data(:,12),'k-','LineWidth',2);
% ylabel('Yaw (deg)');
% subplot(4,1,2);
% hold on;
% plot(CPT_data(:,1),CPT_enu(:,2),'k-','LineWidth',2);
% ylabel('North (m)');
% subplot(4,1,3);
% hold on;
% plot(CPT_data(:,1),CPT_enu(:,1),'k-','LineWidth',2);
% ylabel('East (m)');
% subplot(4,1,4);
% hold on;
% plot(CPT_data(:,1),CPT_enu(:,3),'k-','LineWidth',2);
% ylabel('Down (m)');
% 
figure(3);
subplot(4,1,1);
title('Acceleration')
hold on;
plot(CPT_data(:,1),CPT_data(:,12),'k-','LineWidth',2);
ylabel('Yaw (deg)');
subplot(4,1,2);
hold on;
plot(IMU_data(:,1),IMU_data(:,4),'b-','LineWidth',2);
ylabel('X (m/s^2)');
subplot(4,1,3);
hold on;
plot(IMU_data(:,1),IMU_data(:,5),'b-','LineWidth',2);
ylabel('Y (m/s^2)');
subplot(4,1,4);
hold on;
plot(IMU_data(:,1),-IMU_data(:,6),'b-','LineWidth',2);
ylabel('Z (m/s^2)');
% 
% figure(4);
% subplot(4,1,1);
% title('Angular Rate')
% hold on;
% plot(CPT_data(:,1),CPT_data(:,12),'k-','LineWidth',2);
% ylabel('Yaw (deg)');
% subplot(4,1,2);
% hold on;
% plot(IMU_data(:,1),IMU_data(:,7).*R2D,'b-','LineWidth',2);
% ylabel('X (deg/s)');
% subplot(4,1,3);
% hold on;
% plot(IMU_data(:,1),IMU_data(:,8).*R2D,'b-','LineWidth',2);
% ylabel('Y (deg/s)');
% subplot(4,1,4);
% hold on;
% plot(IMU_data(:,1),-IMU_data(:,9).*R2D,'b-','LineWidth',2);
% ylabel('Z (deg/s)');
