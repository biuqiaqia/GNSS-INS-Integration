function [pos_est,pr_resi,DOP,dtR_est] = WLS_positioning_2(sv_data,el,cnr,uncertainty,pos_ini)

% -------------------------------------------------------------------------
% Least square positioning
% Multi-constellation (GPS/GLONASS/GALILEO/BEIDOU)
% Pseudorange residual and DOP supported
% Input:    sv_data =   column 1: PRN code
%                       column 2: Pseudorange
%                       column 3-5: Satellite position 
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           pr_resi =   pseudorange residual
%           DOP     =   HDOP,VDOP
%
% by GH.Zhang 2019/03/07
% -------------------------------------------------------------------------

% RMSE weighting
snr_1=50;
snr_0=20;
snr_A=50;
snr_a=30;
q_R = 1 ./ (sin(el * pi/180).^2) .* (10.^(-(cnr-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(cnr-snr_1)+1));
q_R(cnr >= snr_1) = 1;
Q = diag(q_R.*(uncertainty));
invQ=diag((diag(Q).^-1));


% uncert_weight_mapping = [1, 0.8, 0.6, 0.4, 0.3];
% uncert_weight = uncert_weight_mapping(uncertainty + 2);  % MATLAB 索引从1
% q_R_combined = q_R ./ uncert_weight';  % 组合后的“伪方差”
% Q_combined = diag(q_R_combined);
% invQ = diag((diag(Q_combined).^-1));

%% Step 2: 权重平衡压缩（收敛聚集高权重卫星）


%% tanh
% % 原始权重 = 1 / q_R （因为 q_R 是方差因子）
% original_weight = 1 ./ q_R;  % +eps 防止除零
% 
% % 设置压缩参数
% threshold = 0.05;                % 权重阈值：大于此值才压缩
% 
% mask = original_weight > threshold;
% if any(mask)
%     w_high = original_weight(mask);
%     w_mean = mean(w_high);      % 计算高权重组的均值
%     scale = 0.05;                 % 可调参数：越大压缩越弱
% 
%     % 中心化
%     centered = w_high - w_mean;
% 
%     % 用 tanh 压缩偏离程度
%     compressed = tanh(centered / scale) * scale;
% 
%     % 加回均值
%     adjusted_high = w_mean + compressed;
% 
%     % 更新权重
%     adjusted_weight = original_weight;
%     adjusted_weight(mask) = adjusted_high + 1;
% else
%     adjusted_weight = original_weight;
% end
% 
% % 转回方差因子
% q_R_adjusted = 1 ./ adjusted_weight;
% Q = diag(q_R_adjusted);
% invQ = diag(adjusted_weight);

original_weight = 1 ./ q_R;  % +eps 防止除零

% 设置压缩参数
threshold = 0.01;                % 权重阈值：大于此值才压缩
compression_factor = 2.0;      % 可调：越大压缩越弱（建议 1.0~5.0）

% 对高权重进行对数压缩
adjusted_weight = original_weight;
mask = original_weight > threshold;

if any(mask)
    ratio = original_weight(mask) / threshold;
    adjusted_weight(mask) = threshold + compression_factor * log10(ratio);
    % 安全保护：确保权重为正
    adjusted_weight(mask) = max(adjusted_weight(mask), 1e-6);
end

% 可选：归一化保持总权重（取消注释启用）
% total_orig = sum(original_weight);
% total_adj = sum(adjusted_weight);
% if total_adj > 0
%     adjusted_weight = adjusted_weight * total_orig / total_adj;
% end

% 转回方差因子
q_R_adjusted = 1 ./ adjusted_weight;

% 构造新权重矩阵
Q = diag(q_R_adjusted);        % 方差矩阵
invQ = diag(adjusted_weight);  % 权重矩阵（直接使用压缩后权重）


constellation = zeros(size(sv_data,1),1);
for idm = 1:1:size(sv_data,1)
    if sv_data(idm,1) <= 32
        constellation(idm) = 1; 
        H_t(idm,:) = [1,0,0,0];
    % elseif sv_data(idm,1) > 32 && sv_data(idm,1) <= 56
    %     constellation(idm) = 2;
    %     H_t(idm,:) = [1,1,0,0];
    % elseif sv_data(idm,1) > 56 && sv_data(idm,1) <= 86
    %     constellation(idm) = 3;
    %     H_t(idm,:) = [1,0,1,0];
    % elseif sv_data(idm,1) > 86 && sv_data(idm,1) <= 123
    elseif sv_data(idm,1) > 32 %by ZY
        constellation(idm) = 4;
        H_t(idm,:) = [1,0,0,1];
    % elseif sv_data(idm,1) > 123 && sv_data(idm,1) <= 127
    %     constellation(idm) = 1;
    %     H_t(idm,:) = [1,0,0,0];
    end
end
H_t = H_t(:,ismember([1;2;3;4],constellation));

if size(sv_data,1)>=(size(H_t,2)+3)
    dt = zeros(size(H_t,2),1);
    delta_pos = [1e9; 1e9; 1e9; dt];

    while norm(delta_pos(1:3,1))>1e-3
        clear H;
        for id_sv = 1:1:size(sv_data,1)
            pos_pr0(id_sv) = norm(sv_data(id_sv,3:5)-pos_ini);
            % y(id_sv,1) = sv_data(id_sv,2) - pos_pr0(id_sv) - dt(1,1);

            y(id_sv,1) = sv_data(id_sv,2) - 40 .*uncertainty(id_sv) - pos_pr0(id_sv) - dt(1,1);
            
            H_v(id_sv,1:3) = (pos_ini - sv_data(id_sv,3:5))./pos_pr0(id_sv);
            %DOP
            sv_enu = xyz2enu(sv_data(id_sv,3:5),pos_ini);
            H_dop(id_sv,1:3) = sv_enu'./pos_pr0(id_sv);
        end
        H = [H_v,H_t];
        N = (H'*(invQ)*H);% Normal matrix
        delta_pos = (N^-1)*H'*(invQ)*y;
%         delta_pos = inv(H'*H)*H'*y;
        pos_ini = pos_ini + delta_pos(1:3,1)';
        dt = dt*1 + delta_pos(4:end,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;%residual
    H_dop = [H_dop,H_t];
    Q = (H_dop'*H_dop)^-1;
    DOP = [sqrt(Q(1,1)+Q(2,2)),sqrt(Q(3,3))];%HDOP,VDOP
    dtR_est = H_t*dt;
else
    pos_est = [nan,nan,nan];
    pr_resi = nan(size(sv_data,1),1);
    DOP = [nan,nan];
    dtR_est = nan;
end

















