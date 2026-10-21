function [Doppler] = True_doppler_by_velocity(sv_data,true_v,pos_ini)

% -------------------------------------------------------------------------
% Least square doppler reverse calculation
% Multi-constellation (GPS/GLONASS/GALILEO/BEIDOU)
%
% by GH.Zhang 2019/03/07
% -------------------------------------------------------------------------

constellation = zeros(size(sv_data,1),1);
for idm = 1:1:size(sv_data,1)
    if sv_data(idm,1) <= 32
        constellation(idm) = 1; 
        H_t(idm,:) = [1,0,0,0];
    elseif sv_data(idm,1) > 32 && sv_data(idm,1) <= 56
        constellation(idm) = 2;
        H_t(idm,:) = [1,1,0,0];
    elseif sv_data(idm,1) > 56 && sv_data(idm,1) <= 86
        constellation(idm) = 3;
        H_t(idm,:) = [1,0,1,0];
    elseif sv_data(idm,1) > 86 && sv_data(idm,1) <= 123
        constellation(idm) = 4;
        H_t(idm,:) = [1,0,0,1];
    elseif sv_data(idm,1) > 123 && sv_data(idm,1) <= 127
        constellation(idm) = 1;
        H_t(idm,:) = [1,0,0,0];
    end
end
H_t = H_t(:,ismember([1;2;3;4],constellation));

if size(sv_data,1)>=(size(H_t,2)+3)
    for id_sv = 1:1:size(sv_data,1)
        pos_pr0(id_sv) = norm(sv_data(id_sv,3:5)-pos_ini);
        H_v(id_sv,1:3) = (pos_ini - sv_data(id_sv,3:5))./pos_pr0(id_sv);
    end
    dop_R = sv_data(:,2);
    lambda = sv_data(:,9);
    VS = sv_data(:,6:8);
    H = [H_v,H_t];
    y = (dop_R .* lambda)-sum(H_v(:,1:3).*VS,2);
    x = inv(H'*H)*H'*y;
    VR  = x(1:3);
    dtRV = H_t*x(4:end,1);
    
    Tx = [true_v;x(4:end,1)];
    Ty = H*Tx;
    Doppler = (Ty+sum(H_v(:,1:3).*VS,2))./lambda;
else
    Doppler = ones(size(sv_data,1),1).*nan;
end
