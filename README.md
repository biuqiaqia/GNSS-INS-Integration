# GNSS/INS Integration

MATLAB research implementation of loosely coupled GNSS/INS integration.

The processing chain combines IMU mechanization in ECEF, GNSS weighted least-squares position and velocity solutions, and a Kalman filter update. The main experimental entry point is `LCKF_main_lite.m`.

## Main components

- `Nav_equations_ECEFdelete1.m`: inertial navigation propagation.
- `WLS_positioning.m` and `WLS_velocity.m`: GNSS measurement solutions.
- `LC_KF_Epoch.m`: loosely coupled Kalman-filter update.
- `Initialize_LC_P_matrix.m`: state-covariance initialization.
- `LCKF_INS_GNSS_GH.m`: alternate integration workflow.

Raw GNSS, IMU, CPT reference trajectories, and generated outputs are intentionally not stored in this repository.

