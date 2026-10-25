# GNSS/INS Integration

MATLAB research code for loosely coupled and tightly coupled GNSS/inertial navigation, recovered classifier-fusion experiments, and a maintained class-conditioned filter interface. Six GNSS/IMU/reference experiment triplets accompany the original MAT example. The baseline propagates ECEF position, velocity, and attitude from IMU measurements, applies GNSS Kalman-filter corrections, and compares standalone GNSS, LCKF, and TCKF trajectories.

> This repository preserves experimental implementations. The classifier interface and four input triplets passed MATLAB loading/numerical checks; full trajectory performance and online inference remain unvalidated. Additional GNSS utilities are needed by the historical drivers. See Sections 13–14 for the maintained interface, original source, and data readiness.

## Table of Contents

1. [Research Objectives](#1-research-objectives)
2. [Integration Methods and State Definitions](#2-integration-methods-and-state-definitions)
3. [Input Data and Coordinate Conventions](#3-input-data-and-coordinate-conventions)
4. [Experiment and Data Preparation](#4-experiment-and-data-preparation)
5. [Integration and Kalman Filter Files](#5-integration-and-kalman-filter-files)
6. [GNSS Position, Velocity, and Doppler Files](#6-gnss-position-velocity-and-doppler-files)
7. [Inertial Navigation and Coordinate Utilities](#7-inertial-navigation-and-coordinate-utilities)
8. [Outputs and Evaluation](#8-outputs-and-evaluation)
9. [Recommended Reading and Execution Order](#9-recommended-reading-and-execution-order)
10. [Dependencies](#10-dependencies)
11. [Known Issues and Reproducibility](#11-known-issues-and-reproducibility)
12. [Third-party Code and Related Projects](#12-third-party-code-and-related-projects)
13. [Classifier-assisted Fusion](#13-classifier-assisted-fusion)
14. [Curated Experimental Data](#14-curated-experimental-data)

## 1. Research Objectives

The project investigates how inertial propagation and GNSS corrections can be combined to maintain a navigation solution between GNSS epochs. Its main objectives are to:

1. propagate an ECEF navigation solution at the IMU sampling rate;
2. estimate GNSS position and velocity from pseudorange and Doppler observations;
3. compare position/velocity-level fusion with observation-level fusion;
4. estimate accelerometer and gyroscope biases through closed-loop filtering; and
5. evaluate integrated trajectories against CPT reference position, velocity, and attitude data.
6. (GNSS signal classifier/GNSS error map could be integrated)

The implementation uses reference information for initialization and some diagnostic calculations. It should therefore be understood as an experimental comparison framework, not an autonomous initialization system.

## 2. Integration Methods and State Definitions

| Method | GNSS information used in the update | Main driver | Update function |
|---|---|---|---|
| LCKF | ECEF receiver position and velocity from weighted least squares | `LCKF_INS_GNSS_GH.m` | `LC_KF_Epoch.m` |
| TCKF | Per-satellite corrected pseudorange and Doppler-derived range rate, together with satellite position and velocity | `TCKF_INS_GNSS_GH.m` | `TC_KF_Epoch.m` |

Both methods call `Nav_equations_ECEFdelete1.m` for inertial propagation. Their update functions propagate the error covariance, form measurement innovations, calculate the Kalman gain, and feed corrections back into attitude, velocity, position, and IMU biases.

### Error-state ordering

| Indices | Quantity | Units |
|---|---|---|
| 1–3 | Attitude error | rad |
| 4–6 | ECEF velocity error | m/s |
| 7–9 | ECEF position error | m |
| 10–12 | Accelerometer bias error | m/s² |
| 13–15 | Gyroscope bias error | rad/s |
| 16 | Receiver clock offset, TCKF only | m |
| 17 | Receiver clock drift, TCKF only | m/s |

LCKF uses 15 states. The active transition and measurement matrices in TCKF use 17 states. Parts of the code attempt to introduce additional constellation clock states at indices 18–19, but this extension is inconsistent with the active 17-state update; see Section 11.

## 3. Input Data and Coordinate Conventions

### `open3s.mat`

The bundled prepared experiment contains the following variables, verified by inspecting the MAT-file metadata:

| Variable | Stored shape | Role |
|---|---|---|
| `CPT_data` | 131 × 12 double | Reference timing, geodetic position, NED velocity, and attitude |
| `GNSS_data` | 131 × 15 cell | Epoch-wise GNSS observations and derived quantities |
| `IMU_data` | 13,056 × 13 double | IMU timing, acceleration, angular rate, and orientation fields |

These are stored dimensions before the main script trims the experiment. Raw ROS bags, receiver logs, navigation files, and separate reference trajectories are not bundled. The MAT-file itself contains experimental navigation data; it is not a synthetic fixture.

### `CPT_data`

| Columns | Interpretation in the code |
|---|---|
| 1 | Reference timestamp used for synchronization |
| 2–3 | GPS week and GPS time-of-week in the loader's reference format |
| 4–6 | Latitude and longitude in degrees; height in metres |
| 7–9 | North, east, and down velocity in m/s |
| 10–12 | Roll, pitch, and heading in degrees, converted to radians where used for initialization |

The ROS loader can append reference inertial fields beyond column 12, but those extra columns are absent from the included example.

### `GNSS_data`

Each row represents one GNSS epoch. The following cells are relevant to the current processing chain:

| Cell column | Contents |
|---|---|
| 1 | GPS time-of-week |
| 2 | Receiver position ordered as longitude, latitude, height |
| 3 | Observed PRN list |
| 4–5 | Observed satellite elevation and azimuth |
| 6–8 | Auxiliary satellite PRN/elevation/azimuth arrays referenced by preprocessing |
| 9 | Per-satellite feature matrix; generated by the loader as PRN, satellite count, availability ratio, C/N0, elevation, normalized pseudorange residual, range-rate consistency, and range-rate proxy |
| 10 | `[PRN, corrected pseudorange, satellite ECEF X, Y, Z]` |
| 11 | `[Doppler, satellite ECEF VX, VY, VZ, wavelength]` |
| 12 | Timestamp aligned with the IMU/reference clock |
| 13 | Receiver ECEF velocity |
| 14 | Optional reference-derived double-difference pseudorange errors |
| 15 | Workflow-dependent intermediate matrix: the loader stores corrected pseudorange data here, while LCKF rebuilds it as a nine-column Doppler/geometry matrix |

Rows of the feature, pseudorange, and Doppler matrices must refer to the same satellites in the same order. The integration drivers obtain the weights through columns 1–2 of cell 9; this does not match the loader's documented-by-assignment C/N0/elevation positions, and must be reconciled before running a new dataset.

### `IMU_data`

The integration functions read acceleration from columns 4–6 and angular rate from columns 7–9. They map these to the body frame as `[column 5, column 4, -column 6]` and `[column 8, column 7, -column 9]`. The main script also reads columns 10–13 with `quat2eul` for a diagnostic attitude conversion.

Time handling is currently inconsistent: the main script and TCKF use column 1, while LCKF uses column 2. Both timestamps must be checked against GNSS cell 12 before comparison. Acceleration and angular rate supplied to mechanization must be in m/s² and rad/s respectively.

Internal position and velocity propagation uses ECEF metres and m/s. Local velocity outputs use north-east-down (NED). Geodetic conversion functions use radians internally; the result arrays store latitude, longitude, and Euler angles in degrees.

## 4. Experiment and Data Preparation

### `LCKF_main_lite.m`

Main entry point for both integration methods despite its LCKF-only filename. It loads `open3s.mat`, trims all streams using `Stime` and `Etime`, estimates a static IMU bias from the initial interval, compensates acceleration for gravity during that calibration, defines LC/TC filter parameters, aligns stream starts, calls both drivers, calculates trajectory errors, and plots all solutions over a map.

The stored experiment settings are `Stime = 1`, `Etime = 103`, and `IMU_static_time = 30`. The script removes selected PRNs from auxiliary cells 6–8; it does not apply that same operation to the observation matrices in cells 10–11. Configuration values are embedded in the script rather than loaded from a separate configuration file.

### `rosbag_loading_open.m`

Historical preprocessing script for reference, GNSS, and IMU ROS topics. It reads broadcast navigation, computes satellite states and observation corrections through goGPS helpers, builds epoch-wise GNSS cells, downsamples GNSS to reference epochs, computes optional reference-station corrections, and saves `CPT_data`, `GNSS_data`, and `IMU_data` under `data/mat/`.

The script is an unfinished preparation template: experiment path selections are commented out, several required variables are not initialized in the active code, and the active IMU-loading loop fills only the timestamp columns. It does not currently populate the complete 13-column IMU input required by the main experiment. Its final MAT filename is based on `file_name`; it does not automatically recreate `open3s.mat`.

### `.gitignore`

Excludes observation exports, raw receiver files, ROS bags, archives, autosaves, and generated data directories. Explicit exceptions track the curated `experiments/*/GNSS.csv`, `IMU.csv`, and `CPT.csv` inputs. MAT-files are not globally ignored, so `open3s.mat` remains tracked.

### `README.md`

Project handover document describing the research purpose, input contracts, file responsibilities, reading order, dependencies, and current implementation limits.

## 5. Integration and Kalman Filter Files

### `LCKF_INS_GNSS_GH.m`

Runs the loosely coupled IMU loop. It initializes position and velocity from GNSS weighted least squares and attitude from the CPT reference, corrects IMU measurements with estimated biases, propagates the ECEF solution, and applies a position/velocity update at available GNSS epochs. It returns the ten-column integrated trajectory and standalone GNSS latitude/longitude/height estimates.

The function also calculates reference-derived velocity diagnostics and rebuilds GNSS cell 15 for Doppler processing. These intermediate modifications are local to the function.

### `TCKF_INS_GNSS_GH.m`

Runs the tightly coupled IMU loop. It initializes GNSS position and receiver clock quantities using unweighted least squares, initializes attitude from CPT, propagates the ECEF solution, and assembles each update matrix as `[pseudorange, range rate, satellite XYZ, satellite velocity XYZ, PRN]`. Range rate is calculated as negative Doppler multiplied by wavelength.

It also computes standalone WLS positions for comparison. Reference-derived Doppler is calculated for diagnostics, but the assignments that would replace measured Doppler are commented out. The satellite-count gate is currently bypassed by `|| 1`.

### `LC_KF_Epoch.m`

One 15-state loosely coupled filter epoch. It builds the first-order error transition, propagates covariance with gyro/accelerometer and bias process noise, forms a six-component ECEF position/velocity innovation, and performs closed-loop correction. Measurement covariance uses the fixed `pos_meas_SD` and `vel_meas_SD` values from the configuration. The active attitude correction includes an additional factor of one-half relative to the corresponding TCKF correction.

### `TC_KF_Epoch.m`

One tightly coupled filter epoch. It predicts per-satellite range and range rate, adds receiver clock terms, and updates navigation, biases, and clock states using 17-dimensional matrices. Variances default to `pseudo_range_SD` and `range_rate_SD`. Optional `observation_pseudorange_SD` and `observation_range_rate_SD` configuration vectors supply one positive standard deviation per satellite. The classifier adapter uses this interface; existing calls retain fixed-noise behavior.

### `Initialize_LC_P_matrix.m`

Creates the 15 × 15 diagonal initial covariance from attitude, velocity, position, accelerometer-bias, and gyro-bias standard deviations in `LC_KF_config`.

### `Initialize_TC_P_matrix.m`

Creates the initial navigation and clock covariance. Its default allocation is 17 × 17; setting `bds_tag = 1` writes indices up to 19 and expands the matrix. The main script currently enables this flag, which conflicts with the active 17-state propagation in `TC_KF_Epoch.m`.

## 6. GNSS Position, Velocity, and Doppler Files

### `WLS_positioning.m`

Iterative weighted least-squares receiver positioning from `[PRN, pseudorange, satellite XYZ]`. It constructs observation weights using the supplied elevation and C/N0, estimates position and constellation-dependent clock terms, and returns position, residuals, horizontal/vertical DOP, and per-satellite clock estimates. Iteration stops when the position increment falls below 1 mm; no maximum iteration count is imposed.

### `WLS_velocity.m`

Weighted Doppler velocity solution. Its nine-column input is `[PRN, Doppler, satellite XYZ, satellite velocity XYZ, wavelength]`. It uses the same elevation/C/N0 weight model as WLS positioning and returns ECEF receiver velocity plus clock-rate terms.

### `Least_square_positioning.m`

Unweighted counterpart of the positioning solver, used for TCKF initialization and preprocessing diagnostics. It also estimates clock terms and returns residuals and DOP. Its insufficient-satellite branch does not assign the fourth output, so callers requesting clock estimates need to handle that case.

### `Least_square_velocity.m`

Unweighted Doppler velocity and clock-rate solver used during TCKF initialization. Its final velocity sign differs from `WLS_velocity.m`; the observation sign convention needs to be verified before interchanging the two functions.

### `True_doppler_by_velocity.m`

Diagnostic Doppler reconstruction using a supplied reference ECEF velocity, satellite geometry, wavelength, and estimated clock-rate terms. It is called by the TCKF driver to calculate a reference-based Doppler sequence; that sequence is not substituted into the active measurement update.

## 7. Inertial Navigation and Coordinate Utilities

| File | Responsibility |
|---|---|
| `Nav_equations_ECEFdelete1.m` | Propagates attitude using Earth rotation and body angular increments, transforms specific force to ECEF, integrates velocity with gravity and Coriolis terms, and integrates position using average velocity. Also returns the acceleration term for diagnostics. |
| `Gravity_ECEF.m` | Computes ECEF gravity including the WGS84 gravitational parameter, J2 term, and centrifugal acceleration. |
| `Skew_symmetric.m` | Converts a three-component vector into its cross-product matrix. |
| `Euler_to_CTM.m` | Converts roll, pitch, and yaw in radians into a coordinate transformation matrix. Callers transpose it where a body-to-navigation matrix is needed. |
| `CTM_to_Euler.m` | Extracts roll, pitch, and yaw in radians from a transformation matrix. |
| `NED_to_ECEF.m` | Converts geodetic position, NED velocity, and body-to-NED attitude to ECEF position, velocity, and attitude. |
| `NED_to_ECEFdebug.m` | Same conversion with the additional ECEF-to-NED transformation matrix returned for diagnostics. |
| `ECEF_to_NED.m` | Converts ECEF position, velocity, and attitude to geodetic position and the local NED frame. |
| `pv_ECEF_to_NED.m` | Position/velocity-only counterpart used during initialization. |
| `plot_openstreetmap.m` | Downloads map tiles for the current longitude/latitude figure axes and draws a configurable map background. Requires network access for uncached tiles. |

These coordinate utilities do not replace the separately called `llh2xyz`, `xyz2llh`, `xyz2enu`, and `enu2xyz` functions, which are absent from this repository.

## 8. Outputs and Evaluation

Both integration drivers return a trajectory with this layout:

| Columns | Output |
|---|---|
| 1 | IMU-loop timestamp |
| 2–4 | Latitude in degrees, longitude in degrees, height in metres |
| 5–7 | North, east, down velocity in m/s |
| 8–10 | Roll, pitch, yaw in degrees; negative yaw is shifted into the positive heading range |

`GNSS_llh` contains epoch-wise standalone GNSS latitude, longitude, and height. Bias and covariance-standard-deviation histories are built inside the drivers but are not exposed in their return signatures.

The main script creates `Error_2D` and `Error_3D` in the workspace, each ordered as `[reference time, GNSS error, LCKF error, TCKF error]`. Three-dimensional error is an ECEF position distance. The quantity called 2D error is calculated by replacing estimated height with reference height before taking that distance; it is not an explicit east/north projection.

The script prints RMS summaries and plots reference, standalone GNSS, LCKF, and TCKF trajectories. Integrated estimates are matched to reference epochs by nearest timestamp, without a maximum time-gap check. It does not automatically export a result table or save the generated figures.

No numerical navigation accuracy is claimed here. The classifier addition passed MATLAB R2024b checks for real-data loading and synthetic one-epoch filter equivalence. Full historical trajectory comparisons and online model inference have not been validated.

## 9. Recommended Reading and Execution Order

1. Read `LCKF_main_lite.m` for experiment settings, stream trimming, calibration, and evaluation.
2. Inspect the three variables in `open3s.mat` and compare them with Section 3.
3. Read `LCKF_INS_GNSS_GH.m` alongside `WLS_positioning.m`, `WLS_velocity.m`, and `LC_KF_Epoch.m`.
4. Read `TCKF_INS_GNSS_GH.m` alongside `TC_KF_Epoch.m` and `Initialize_TC_P_matrix.m`.
5. Read `Nav_equations_ECEFdelete1.m` and the coordinate utilities to follow inertial propagation.
6. Consult `rosbag_loading_open.m` only when preparing a new experiment, after completing its missing configuration and data extraction.

For a local reproduction, add this repository and the required external GNSS helpers to the MATLAB path. Inspect the input before launching the main script:

```matlab
whos('-file', 'open3s.mat')
which llh2xyz -all
which xyz2llh -all
which xyz2enu -all
which enu2xyz -all
```

Resolve the clock-state dimensions, input-column conventions, and timing issues listed below, then run:

```matlab
LCKF_main_lite
```

The script begins with `clear`, so configuration changes must be made in the script or in a deliberate refactoring rather than entered only in the workspace.

## 10. Dependencies

- MATLAB for the experiment, matrix solvers, plotting, and MAT-file inputs. The original MATLAB release has not been established.
- External GNSS coordinate helpers: `llh2xyz`, `xyz2llh`, `xyz2enu`, and `enu2xyz`.
- An installation providing `quat2eul` for the main script's quaternion diagnostic, and `rms` for its printed evaluation.
- ROS Toolbox or a compatible MATLAB installation providing `rosbag`, `select`, and `readMessages`, when using the ROS preprocessing script.
- goGPS and related RINEX, satellite-state, atmospheric-correction, and least-squares helpers for preprocessing. They are not bundled in this repository.
- `finv` support for the preprocessing script's statistics table, plus the separately referenced `DD_pr_error_from_Ref` helper and reference-station data if that branch is used.
- Network access for the optional map-tile background.

## 11. Known Issues and Reproducibility

- **Clock-state dimensions:** the main script sets `bds_tag = 1`, causing covariance initialization to expand to 19 states, while the TC transition and observation matrices remain 17-state. Branches for non-GPS PRNs reference clock states 18–19 even though their propagation is commented out.
- **Time columns:** LCKF reads IMU column 2, while TCKF and the main script use column 1. Equal meaning and alignment cannot be assumed for another dataset.
- **GNSS weight columns:** the drivers pass feature-matrix column 2 as elevation and column 1 as C/N0. The ROS loader actually assigns satellite count and PRN to those columns. The prepared input and the loader require an explicit schema reconciliation.
- **Satellite availability:** TCKF bypasses its calculated satellite-count requirement with `|| 1`, so insufficient or singular geometry can reach the solvers.
- **Reference-assisted initialization:** both drivers use CPT attitude and reference position information; reference-derived velocity/Doppler calculations also remain in the experimental code.
- **Constellation handling:** several solvers group every PRN greater than 32 into one category, while other routines distinguish multiple PRN ranges. Their clock models are not consistently implemented across initialization and filtering.
- **Numerical and sign conventions:** solvers use explicit matrix inverses, positioning iterations have no hard limit, weighted and unweighted velocity signs differ, and the LC attitude correction contains an extra one-half factor. These choices need validation before quantitative conclusions.
- **Preprocessing completeness:** experiment paths and variables remain unresolved in `rosbag_loading_open.m`; full IMU measurements are not extracted by its active loading loop.
- **Evaluation:** nearest-time matching has no tolerance, the reported 2D distance uses a height substitution, and no independent execution log or validated result artifact is included.

## 12. Third-party Code and Related Projects

Several filter, mechanization, gravity, and coordinate functions retain Paul Groves' copyright notices and refer to the BSD license accompanying the software for *Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems*, Second Edition. Driver and loader headers also identify modifications by G. H. Zhang. Preserve those notices; a repository-level license and the referenced `license.txt` are not currently included.

The related [GNSS Signal Classifier](https://github.com/biuqiaqia/GNSS-Signal-Classifier) project studies per-satellite signal conditions. Recovered classifier work and its current integration interfaces are documented below. The later error-map integration is outside this repository update.

## 13. Classifier-assisted Fusion

The classifier-related work has two paths: historical fusion experiments reading `Type`/`Uncertainty` from CSV, and a MATLAB-to-Python inference prototype. Neither establishes that the bundled CSV labels were generated by the archived neural model. This addition preserves that distinction and provides a consistent class-conditioned filtering interface.

### Original Source Code

The files in `classifier/archive/` are the project's original source code, recovered from the author's existing project folders. All six MATLAB files were copied byte-for-byte without modification. The experimental CSV inputs are also original project data; only their repository filenames were standardized. Existing upstream copyright and attribution notices remain in the source files.

### Supplementary Code

The following files were newly written by Codex during repository organization. They are supplementary code, not the author's original research implementation.

| File | Role |
|---|---|
| `classifier/load_classifier_experiment.m` | Reads a GNSS/IMU/CPT triplet, aligns GNSS time using an exact CPT GPS-time lookup, and reads IMU TimePoint. Defaults to requiring Predicted_Type. Explicitly selecting Type enables original-label replay with a warning. |
| `classifier/build_classifier_measurements.m` | Constructs a canonical 13-column matrix: pseudorange, range rate, satellite XYZ, satellite velocity XYZ, PRN, selected label, stored uncertainty, C/N0, elevation. |
| `classifier/TC_KF_Epoch_classifier.m` | Converts labels 0–3 into caller-supplied pseudorange/range-rate standard deviations and calls the existing TC filter. Supports GPS PRNs 1–32 with 17 states. |
| `classifier/verify_classifier_interfaces.m` | Checks four real triplets, rejection of missing predictions and ambiguous reference times, equivalence to fixed-noise filtering, and the effect of class-dependent covariance. |

`tools/audit_experiments.py` is also supplementary code written by Codex. It checks the original experimental inputs and generates `experiments/audit.json`; that JSON is a generated report, not an original experimental result. README documents and `.gitattributes` were added or updated during the same organization work.

**Supplementary modification to original source:** the root-level `TC_KF_Epoch.m` retains its original filter implementation, with a Codex-added block that accepts optional per-satellite `observation_pseudorange_SD` and `observation_range_rate_SD` values. This block is used by the supplementary classifier interface. Calls without those fields retain the original fixed-noise behavior. The root-level file is therefore original source with a documented supplementary modification, not an unchanged archival copy.

The supplementary interfaces do not impose physical class meanings or class-noise calibration: classifier and dataset versions have different label definitions. Their interface tests validate the supplementary code, not the performance of the original research system.

### Original Source Code: File Details

`classifier/archive/` preserves six original files byte-for-byte. Its [README](classifier/archive/README.md) records provenance and actual code behavior. Do not add the archive recursively to the MATLAB path: the saved calls do not form a compatible standalone version.

The original driver creates baseline/type/uncertainty trajectories and calls variant WLS/TC functions. In recovered copies, the class-noise lookup in TC_KF_Epoch_1 is commented out and observation columns differ between callers and callees. Its old WLS_positioning_1 dependency was later repurposed; that incompatible implementation is excluded. The archive therefore does not establish an AI-assisted navigation gain.

### Reading an experiment

From the repository root:

```matlab
addpath(pwd, fullfile(pwd, 'classifier'));
[CPT_data, GNSS_data, IMU_data, info] = ...
    load_classifier_experiment(fullfile('experiments', 'new'), 'Type');
measurements = build_classifier_measurements(GNSS_data(1,:));
```

This selects original-label replay. The loader retains stored scores and does not relabel, clamp, crop, resample, or overwrite CSVs. Both returned IMU timing columns use TimePoint; GNSS timestamps come from CPT Time rather than an epoch counter. Sensor axes, velocity frames, static initialization, and evaluation intervals still need experiment-specific confirmation.

For actual model-assisted replay, join model outputs as a Predicted_Type column by GPS time and PRN, retaining satellite order, and use the default loader call. Match the classifier's exact feature normalization and label schema. The archived inference prototype uses [Elevation, Azimuth, CNR], but its referenced predict_gnss.py module has not been recovered. Its pred_Uncertainty output merely duplicates the class ID.

### Class-conditioned filter update

Provide positive four-element configuration vectors `class_pseudorange_SD` (metres) and `class_range_rate_SD` (m/s), indexed by numeric label plus one. Initialize GPS covariance with bds_tag = 0 and call TC_KF_Epoch_classifier with the same state arguments as TC_KF_Epoch. The adapter reads labels from observation column 10 and forms diagonal variances as squared standard deviations. It does not use the stored uncertainty column as a probability or change the inertial propagation.

Calibration must be established for the selected label schema; commented historical lookup values are not treated as validated settings. This is an epoch-level integration point, not a replacement for full experiment initialization and evaluation.

```matlab
verify_classifier_interfaces
```

The interface checks passed on MATLAB R2024b Update 6. No trajectory accuracy or end-to-end classifier performance is implied by these checks.

## 14. Curated Experimental Data

Six original triplets are included under experiments/. Their [data guide](experiments/README.md) describes each file, source, and readiness. The [metadata audit](experiments/audit.json) records SHA-256 hashes, label counts, uncertainty ranges, and timing issues.

| Experiment | GNSS rows | GNSS epochs | IMU rows | CPT rows | Input status |
|---|---:|---:|---:|---:|---|
| aaa_test | 1,996 | 287 | 23,463 | 326 | Passes structural loading checks |
| new | 1,187 | 181 | 14,300 | 211 | Passes structural loading checks |
| data1 | 1,341 | 209 | 27,552 | 209 | Passes structural loading checks |
| data3 | 671 | 119 | 27,552 | 250 | Missing CPT/IMU timestamps need repair |
| data5 | 2,509 | 264 | 23,271 | 264 | Passes structural loading checks |
| data6 | 2,465 | 273 | 25,015 | 2,465 | Missing CPT time/position and repeated GPS times |

None contains Predicted_Type. Several uncertainty columns contain continuous or negative values of unconfirmed provenance. Structural readiness does not establish synchronization accuracy, prediction provenance, classifier accuracy, or trajectory performance.
