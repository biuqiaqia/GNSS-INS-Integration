function verify_classifier_interfaces()
% Contract and numerical checks. Does not claim end-to-end AI performance.
folder = fileparts(mfilename('fullpath'));
root = fileparts(folder);
addpath(root, folder);
warning('off','classifier:ReferenceReplay');
cleanup = onCleanup(@() warning('on','classifier:ReferenceReplay')); %#ok<NASGU>
for name = {'aaa_test','new','data1','data5'}
    [c,g,i,info] = load_classifier_experiment(fullfile(root,'experiments',name{1}), 'Type');
    assert(~isempty(c) && ~isempty(g) && ~isempty(i));
    assert(isequal(i(:,1),i(:,2)) && ~info.is_model_prediction);
    observations = build_classifier_measurements(g(1,:));
    assert(size(observations,2) == 13);
    assert(isequal(observations(:,10),g{1,15}));
    fprintf('LOAD OK: %s (%d epochs)\n',name{1},size(g,1));
end
must_fail(@() load_classifier_experiment(fullfile(root,'experiments','new')), ...
    'classifier:MissingLabels');
must_fail(@() load_classifier_experiment(fullfile(root,'experiments','data6'),'Type'), ...
    'classifier:MissingTime');
must_fail(@() load_classifier_experiment(fullfile(root,'experiments','data3'),'Type'), ...
    'classifier:MissingTime');

% Verify that equal class noise reproduces the pre-existing TC calculation,
% and a class-dependent calibration actually changes the covariance.
config.gyro_noise_PSD = 0.005^2;
config.accel_noise_PSD = 0.008^2;
config.accel_bias_PSD = 1e-5;
config.gyro_bias_PSD = 4e-11;
config.clock_phase_PSD = 1;
config.clock_freq_PSD = 1;
config.pseudo_range_SD = 25;
config.range_rate_SD = 2;
sat = [20200000,14000000,17000000; -15000000,21000000,12000000; ...
    18000000,-14000000,19000000; 16000000,17000000,-18000000];
r0 = [6378137;0;0];
rho = sqrt(sum((sat-r0').^2,2))+10;
m = [rho,zeros(4,1),sat,zeros(4,3),(1:4)',(0:3)',zeros(4,1),40*ones(4,1),45*ones(4,1)];
args = {1,eye(3),zeros(3,1),r0,zeros(6,1),[0,0],eye(17),zeros(3,1),0};
[~,~,~,~,~,baseline] = TC_KF_Epoch(m(:,1:9),args{:},config);
config.class_pseudorange_SD = [25,25,25,25];
config.class_range_rate_SD = [2,2,2,2];
[~,~,~,~,~,equal_noise] = TC_KF_Epoch_classifier(m,args{:},config);
assert(norm(baseline-equal_noise,'fro') < 1e-12);
config.class_pseudorange_SD = [5,10,20,40];
[~,~,~,~,~,class_noise] = TC_KF_Epoch_classifier(m,args{:},config);
assert(all(isfinite(class_noise(:))) && norm(class_noise-equal_noise,'fro') > 1e-8);
bad = m; bad(1,10) = 4;
must_fail(@() TC_KF_Epoch_classifier(bad,args{:},config),'');
bad = m; bad(1,9) = 87;
must_fail(@() TC_KF_Epoch_classifier(bad,args{:},config),'classifier:GPSOnly');
fprintf('CLASSIFIER INTERFACE CHECKS PASSED\n');
end

function must_fail(call, identifier)
failed = false;
try
    call();
catch exception
    failed = true;
    assert(isempty(identifier) || strcmp(exception.identifier,identifier), ...
        'Unexpected error: %s',exception.identifier);
end
assert(failed,'Expected the input to be rejected.');
end
