function [C,v,r,bias,clock,P] = TC_KF_Epoch_classifier(measurements, ...
    dt,C0,v0,r0,bias0,clock0,P0,specific_force,latitude,config)
% Class-conditioned TC update using the existing 17-state GPS filter.
% Caller must supply calibrated noise standard deviations for labels 0:3:
%   config.class_pseudorange_SD (metres), config.class_range_rate_SD (m/s).
% No class meaning or uncertainty-to-confidence conversion is inferred.
assert(size(measurements,2) == 13 && ~isempty(measurements), ...
    'classifier:ObservationWidth', 'Use build_classifier_measurements to construct observations.');
assert(all(measurements(:,9) >= 1 & measurements(:,9) <= 32), ...
    'classifier:GPSOnly', 'This adapter supports the consistent 17-state GPS path only.');
assert(isequal(size(P0),[17,17]), 'classifier:StateSize', ...
    'Expected a 17-state covariance and bds_tag = 0 at initialization.');
assert(isfield(config,'class_pseudorange_SD') && isfield(config,'class_range_rate_SD'), ...
    'classifier:MissingCalibration', 'Provide explicit class-to-noise calibration for this label version.');
validateattributes(config.class_pseudorange_SD, {'numeric'}, {'numel',4,'finite','positive'});
validateattributes(config.class_range_rate_SD, {'numeric'}, {'numel',4,'finite','positive'});
labels = measurements(:,10);
validateattributes(labels, {'numeric'}, {'finite','integer','>=',0,'<=',3});
assert(all(all(isfinite(measurements(:,1:9)))), ...
    'classifier:InvalidObservation', 'Navigation observations must be finite.');
config.observation_pseudorange_SD = config.class_pseudorange_SD(labels+1);
config.observation_range_rate_SD = config.class_range_rate_SD(labels+1);
[C,v,r,bias,clock,P] = TC_KF_Epoch(measurements(:,1:9),dt, ...
    C0,v0,r0,bias0,clock0,P0,specific_force,latitude,config);
end
