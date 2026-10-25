# Experimental Input Triplets

Each folder contains GNSS.csv, IMU.csv, and CPT.csv. These retain original source bytes; filename normalization is the only import change. Sources, SHA-256 hashes and checks are recorded in [audit.json](audit.json).

| Folder | Original workspace folder | GNSS / IMU / CPT rows | Status |
|---|---|---|---|
| aaa_test | System_ongoing/input/AAAtest/ | 1,996 / 23,463 / 326 | Structurally loadable for original-label replay |
| new | System_ongoing/input/NEW/ | 1,187 / 14,300 / 211 | Structurally loadable for original-label replay |
| data1 | AAAPNT/The Past/data1/ | 1,341 / 27,552 / 209 | Structurally loadable for original-label replay |
| data3 | AAAPNT/The Past/data3/ | 671 / 27,552 / 250 | Missing CPT/IMU timestamps; strict loader rejects it |
| data5 | AAAPNT/The Past/data5/ | 2,509 / 23,271 / 264 | Structurally loadable for original-label replay |
| data6 | AAAPNT/The Past/TCKF/ | 2,465 / 25,015 / 2,465 | Missing CPT time/position and repeated GPS times; strict loader rejects it |

The audit records the original filenames. Duplicate workspace copies were not imported. The alternative GNSS_data_scr.csv score variant is excluded. File grouping does not establish sensor synchronization or independence between experiments.

## File roles and schemas

- GNSS.csv: 34 columns, one satellite observation per row. Includes GPS_Time, PRN, geometry, C/N0, corrected pseudorange, Doppler, wavelength, satellite position/velocity, receiver velocity, WLS/reference position, ray-tracing fields, Type, and Uncertainty. Epoch/PRN keys are unique within each selected file.
- IMU.csv: 13 columns, including named acceleration, angular rate, quaternion, and TimePoint. The first two experiments place TimePoint first; numbered datasets place a packet counter first and TimePoint second. The maintained loader selects TimePoint by name.
- CPT.csv: 12 columns with reference Time, GPS week/time, position, velocity, and roll/pitch/yaw. Velocity fields use VR-style names; their frame must be checked before interpreting them as NED. Data3 has 131 missing CPT GPS_Time values and 3,237 missing IMU TimePoint values. Data6 has 2,192 rows missing CPT Time/position, in addition to repeated numeric GPS_Time values. Missing times are distinct from duplicate valid times; neither is repaired implicitly.

No selected GNSS file contains Predicted_Type or another explicit prediction column. Original Type labels are not automatically AI outputs. aaa_test/new have discrete uncertainty values; numbered datasets have continuous scores, including negative values. Their generation history is not established. Do not interpret them as probabilities or clamp/remap them silently.

## Loading and checks

```matlab
addpath(pwd, fullfile(pwd, 'classifier'));
[CPT_data, GNSS_data, IMU_data, info] = ...
    load_classifier_experiment(fullfile('experiments', 'new'), 'Type');
```

Passing Type explicitly selects original-label replay. Omitting it requires Predicted_Type and rejects the original inputs. Numeric labels are preserved; no equivalence to other label versions is assumed.

GNSS GPS_Time is matched exactly to CPT GPS_Time to obtain CPT Time. Both returned IMU timing columns use TimePoint. Raw CSVs are unchanged. The loader does not calibrate axes, choose static initialization, trim evaluation intervals, or repair ambiguous reference records.

Rebuild metadata with Python, NumPy, and pandas:

```bash
python tools/audit_experiments.py
```

The audit reads inputs and writes metadata only. Full navigation evaluation still requires sensor-frame and timing confirmation, initialization, verified label definitions, and class-noise calibration.
