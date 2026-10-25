# Original Classifier Fusion Source

These six MATLAB files are the project's **Original Source Code**, recovered from the author's existing project folders and copied byte-for-byte without modification. This README is supplementary documentation written by Codex. The new interfaces in the parent classifier folder are **Supplementary Code** written by Codex, not recovered original research code. Add only the parent classifier folder and repository root when using those supplementary interfaces.

| File | Source relative to the original workspace | Actual behavior and limitations |
|---|---|---|
| KF_main_lite_1.m | System_ongoing/input/AAAtest/ | Reads CSVs, creates C/N0/elevation features and Type/Uncertainty cells, calls LCKF and the parallel TC driver, and compares trajectories. Embedded paths and timing assumptions remain. |
| TCKF_INS_GNSS_GH_1.m | AAAPNT/The Past/TCKF/ | Propagates baseline, type, and uncertainty histories. Constructs a 13-column observation matrix. Its original compatible classifier-only WLS_positioning_1 dependency was not recovered. |
| TC_KF_Epoch_1.m | AAAPNT/The Past/TCKF/ | Uses a 19-state clock model, but class-noise lookups are commented out. Reads columns 11/12 as C/N0/elevation, unlike the driver. GPS-only observation matrices may remain 17-column, conflicting with the covariance. |
| TC_KF_Epoch_2.m | AAAPNT/The Past/TCKF/ | Subtracts 50 times column-11 score from pseudorange innovation above threshold 0.01. Also reads column 11 as C/N0. This is an experimental score correction with inconsistent layout and incomplete clock handling. |
| WLS_positioning_2.m | AAAPNT/The Past/TCKF/ | Subtracts 40 times uncertainty from pseudorange. Its initial uncertainty-dependent covariance is overwritten by signal-strength weight compression. |
| predict_gnss_quality.m | System_ongoing/ | Normalizes Elevation/Azimuth/CNR and calls py.predict_gnss.predict_single_epoch. The Python module was not recovered. Returns the class ID again as pred_Uncertainty. No active caller was located in the recovered classifier-only drivers. |

These files document original experiments, not a validated standalone application. The maintained parent-folder adapter explicitly resolves the observation layout and accepts caller-provided class-noise calibration without silently altering historical code.
