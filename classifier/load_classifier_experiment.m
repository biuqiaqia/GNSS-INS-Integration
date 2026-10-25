function [CPT_data, GNSS_data, IMU_data, info] = load_classifier_experiment(experiment_dir, label_field)
% Read a curated triplet without changing the archived CSV files.
% Default input is a genuine Predicted_Type column. Passing 'Type' explicitly
% enables reference-label replay and does not constitute AI inference.
if nargin < 2
    label_field = 'Predicted_Type';
end
g = readtable(fullfile(experiment_dir, 'GNSS.csv'));
i = readtable(fullfile(experiment_dir, 'IMU.csv'));
c = readtable(fullfile(experiment_dir, 'CPT.csv'));
required = {'GPS_Time','PRN','Elevation','Azimuth','CNR','corrected_Pr', ...
    'SaX','SaY','SaZ','Doppler','SaVX','SaVY','SaVZ','lambda', ...
    'VRx','VRy','VRz','Pr_Residual','WLS_Longitude','WLS_Latitude','WLS_Height'};
assert(all(ismember(required, g.Properties.VariableNames)), ...
    'classifier:MissingColumns', 'Missing required GNSS columns.');
assert(ismember(label_field, g.Properties.VariableNames), ...
    'classifier:MissingLabels', ...
    'Requested label column is absent. Original Type is not automatically a model prediction.');
assert(all(ismember({'Time','GPS_Time'}, c.Properties.VariableNames)), ...
    'classifier:ReferenceTime', 'CPT must contain Time and GPS_Time.');
assert(all(isfinite(c.GPS_Time)) && all(isfinite(c.Time)), ...
    'classifier:MissingTime', 'CPT timestamps must be finite. Original missing samples are not repaired.');
assert(numel(unique(c.GPS_Time)) == height(c), ...
    'classifier:AmbiguousTime', 'CPT repeats GPS_Time. Resolve the reference stream before replay.');
assert(ismember('TimePoint', i.Properties.VariableNames), ...
    'classifier:IMUTime', 'IMU must contain an explicit TimePoint column.');
assert(all(isfinite(i.TimePoint)), 'classifier:MissingTime', ...
    'IMU TimePoint contains missing timestamps.');
assert(size(unique([g.GPS_Time,g.PRN], 'rows'),1) == height(g), ...
    'classifier:DuplicateObservation', 'GNSS repeats an epoch/PRN key.');
if strcmp(label_field, 'Type')
    warning('classifier:ReferenceReplay', ...
        'Using original Type labels. This is label-assisted replay, not measured AI performance.');
end
labels = g.(label_field);
validateattributes(labels, {'numeric'}, {'finite','integer','>=',0,'<=',3});
times = unique(g.GPS_Time);
[found, reference_rows] = ismember(times, c.GPS_Time);
assert(all(found), 'classifier:UnmatchedTime', ...
    'Some GNSS epochs have no exact CPT GPS_Time match. No time origin is guessed.');
GNSS_data = cell(numel(times), 16);
for epoch = 1:numel(times)
    rows = g.GPS_Time == times(epoch);
    e = g(rows,:);
    GNSS_data{epoch,1} = times(epoch);
    GNSS_data{epoch,2} = [e.WLS_Longitude(1),e.WLS_Latitude(1),e.WLS_Height(1)];
    GNSS_data{epoch,3} = e.PRN;
    GNSS_data{epoch,4} = e.Elevation;
    GNSS_data{epoch,5} = e.Azimuth;
    GNSS_data{epoch,9} = [e.CNR,e.Elevation];
    GNSS_data{epoch,10} = [e.PRN,e.corrected_Pr,e.SaX,e.SaY,e.SaZ];
    GNSS_data{epoch,11} = [e.Doppler,e.SaVX,e.SaVY,e.SaVZ,e.lambda];
    GNSS_data{epoch,12} = c.Time(reference_rows(epoch));
    GNSS_data{epoch,13} = [e.VRx(1);e.VRy(1);e.VRz(1)];
    GNSS_data{epoch,14} = [e.PRN,e.Pr_Residual];
    GNSS_data{epoch,15} = labels(rows);
    if ismember('Uncertainty', e.Properties.VariableNames)
        GNSS_data{epoch,16} = e.Uncertainty;
    else
        GNSS_data{epoch,16} = nan(height(e),1);
    end
end
CPT_data = table2array(c);
IMU_data = table2array(i);
% Normalize timestamps only in the returned in-memory copy so the two
% historical drivers see the same physical time rather than packet counts.
IMU_data(:,1) = i.TimePoint;
IMU_data(:,2) = i.TimePoint;
assert(all(diff(i.TimePoint) > 0), 'classifier:NonmonotonicIMU', ...
    'IMU TimePoint must be strictly increasing. Original samples are not removed.');
info.label_field = label_field;
info.is_model_prediction = strcmp(label_field, 'Predicted_Type');
info.gnss_time = c.Time(reference_rows);
info.imu_time_range = [min(i.TimePoint),max(i.TimePoint)];
info.reference_time_range = [min(c.Time),max(c.Time)];
info.overlap = [max([min(info.gnss_time),min(i.TimePoint),min(c.Time)]), ...
    min([max(info.gnss_time),max(i.TimePoint),max(c.Time)])];
assert(info.overlap(2) > info.overlap(1), 'classifier:NoOverlap', ...
    'The three streams have no overlapping interval.');
end
