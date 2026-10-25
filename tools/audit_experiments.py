"""Read-only CSV checks; write a reproducible metadata report, not cleaned data."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    'aaa_test': ('System_ongoing/input/AAAtest', 'GNSS_data.csv', 'IMU_data.csv', 'CPT_data.csv'),
    'new': ('System_ongoing/input/NEW', 'GNSS_data.csv', 'IMU_data.csv', 'CPT_data.csv'),
    'data1': ('AAAPNT/The Past/data1', 'GNSS_data1.csv', 'IMU_data1.csv', 'CPT_data1.csv'),
    'data3': ('AAAPNT/The Past/data3', 'GNSS_data3.csv', 'IMU_data3.csv', 'CPT_data.csv'),
    'data5': ('AAAPNT/The Past/data5', 'GNSS_data5.csv', 'IMU_data5.csv', 'CPT_data5.csv'),
    'data6': ('AAAPNT/The Past/TCKF', 'GNSS_data6.csv', 'IMU_data6.csv', 'CPT_data6.csv'),
}


def audit(name, source):
    folder = ROOT / 'experiments' / name
    tables = {kind: pd.read_csv(folder / f'{kind}.csv') for kind in ('GNSS', 'IMU', 'CPT')}
    files = {}
    for index, (kind, table) in enumerate(tables.items(), 1):
        path = folder / f'{kind}.csv'
        files[kind] = {
            'source': f'{source[0]}/{source[index]}',
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'rows': len(table), 'columns': list(table.columns),
            'missing_cells': int(table.isna().sum().sum()),
        }
    g, i, c = (tables[k] for k in ('GNSS', 'IMU', 'CPT'))
    issues = []
    duplicate_gnss = int(g.duplicated(['GPS_Time', 'PRN']).sum())
    duplicate_cpt = int(c.GPS_Time.dropna().duplicated().sum())
    missing_cpt_times = int(c.GPS_Time.isna().sum())
    missing_reference_times = int(c.Time.isna().sum())
    missing_imu_times = int(i.TimePoint.isna().sum())
    bad_imu_steps = int((np.diff(i.TimePoint.to_numpy()) <= 0).sum())
    unmatched = int((~g.GPS_Time.isin(c.GPS_Time)).sum())
    if duplicate_gnss:
        issues.append('Duplicate GNSS epoch/PRN keys')
    if duplicate_cpt:
        issues.append('Repeated CPT GPS_Time makes exact reference-time lookup ambiguous')
    if missing_cpt_times or missing_reference_times or missing_imu_times:
        issues.append('Missing CPT Time/GPS_Time or IMU TimePoint values')
    if c[['Latitude','Longitude','Height']].isna().any().any():
        issues.append('Missing reference positions')
    if bad_imu_steps:
        issues.append('IMU TimePoint is not strictly increasing')
    if unmatched:
        issues.append('GNSS records without an exact CPT GPS_Time match')
    # Do not infer model provenance from an original Type column.
    return {
        'files': files, 'gnss_epochs': int(g.GPS_Time.nunique()),
        'type_counts': {str(k): int(v) for k, v in g.Type.value_counts().sort_index().items()},
        'uncertainty': {'min': float(g.Uncertainty.min()), 'max': float(g.Uncertainty.max()),
                        'unique_values': int(g.Uncertainty.nunique()),
                        'negative_rows': int((g.Uncertainty < 0).sum())},
        'has_predicted_type': 'Predicted_Type' in g,
        'duplicate_gnss_keys': duplicate_gnss, 'repeated_reference_times': duplicate_cpt,
        'missing_reference_gps_times': missing_cpt_times, 'missing_imu_times': missing_imu_times,
        'missing_reference_times': missing_reference_times,
        'nonpositive_imu_time_steps': bad_imu_steps,
        'unmatched_gnss_rows': unmatched,
        'imu_time_range': [float(i.TimePoint.min()), float(i.TimePoint.max())],
        'cpt_time_range': [float(c.Time.min()), float(c.Time.max())],
        'input_readiness': 'needs repair before strict loading' if issues else 'passes structural checks for explicit Type replay',
        'issues': issues,
    }


def main():
    report = {'note': 'Structural checks only. No AI prediction provenance or navigation accuracy is established.',
              'experiments': {name: audit(name, source) for name, source in SOURCES.items()}}
    output = ROOT / 'experiments' / 'audit.json'
    output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    for name, item in report['experiments'].items():
        counts = [item['files'][kind]['rows'] for kind in ('GNSS','IMU','CPT')]
        print(name, counts, item['gnss_epochs'], item['input_readiness'], item['issues'])


if __name__ == '__main__':
    main()
