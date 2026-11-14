"""Recompute counts and method agreement from the published detection coordinates."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import numpy as np
from config import ROOT, RESULTS, MATCH_RADIUS_UM, DENSE_RADIUS_UM
from spatial import match_points, within_circle


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def detections(section, method='model', common_only=True):
    records = read_csv(RESULTS / 'detections' / f'{method}_section_{section:02d}.csv')
    if common_only:
        records = [r for r in records if r['inside_atlas'] == 'True'
                   and r['inside_tissue_envelope'] == 'True']
    points = np.array([[float(r['y_px']), float(r['x_px'])] for r in records]).reshape(-1, 2)
    return points, records


def analyze():
    metadata = read_csv(ROOT / 'data/slices.csv')
    expected = {int(r['section']): r for r in read_csv(RESULTS / 'section_counts.csv')}
    regional = read_csv(RESULTS / 'regional_comparison.csv')
    neighborhoods = read_csv(RESULTS / 'dense_neighborhoods.csv')
    members = read_csv(RESULTS / 'neighborhood_members.csv')
    sections = []
    circle_counts = Counter()
    for rec in metadata:
        sec, um = int(rec['section']), float(rec['um_per_pixel'])
        a, ar = detections(sec)
        b, br = detections(sec, 'alternative')
        pairs = match_points(a, b, um, MATCH_RADIUS_UM)
        shared = int((pairs >= 0).sum())
        result = dict(section=sec, full=len(a), alternative=len(b), shared=shared,
                      full_only=len(a)-shared, alternative_only=len(b)-shared)
        for key, value in result.items():
            if key != 'section':
                assert value == int(expected[sec][key]), (sec, key)
        for method, records in [('full', ar), ('alternative', br)]:
            counts = Counter(int(r['structure_id']) for r in records)
            saved = [r for r in regional if int(r['section']) == sec]
            assert sum(int(r[method]) for r in saved) == len(records)
            for row in saved:
                sid, area = int(row['structure_id']), float(row['common_area_mm2'])
                assert int(row[method]) == counts[sid]
                assert np.isclose(float(row[method+'_per_mm2']), counts[sid]/area)
        for method, points, other in [('full', a, b), ('alternative', b, a)]:
            owned = []
            for row in neighborhoods:
                if int(row['section']) != sec or row['selected_by'] != method:
                    continue
                center = [float(row['center_y_px']), float(row['center_x_px'])]
                idx = within_circle(points, center, um, DENSE_RADIUS_UM)
                other_count = len(within_circle(other, center, um, DENSE_RADIUS_UM))
                assert len(idx) == int(row['count'])
                other_name = 'alternative' if method == 'full' else 'full'
                assert other_count == int(row[other_name])
                saved = [r for r in members if int(r['section']) == sec
                         and r['method'] == method and r['neighborhood'] == row['id']]
                saved_idx = np.array([int(r['common_detection_index']) for r in saved])
                assert np.array_equal(idx, saved_idx)
                saved_points = np.array([[float(r['y_px']), float(r['x_px'])] for r in saved])
                assert np.array_equal(points[idx], saved_points)
                owned.extend(idx.tolist())
                circle_counts[method] += 1
            assert len(owned) == len(set(owned))
        sections.append(result)
    validation = json.loads((RESULTS / 'validation.json').read_text())
    scored = read_csv(RESULTS / 'validation_by_crop.csv')
    for method, metrics in validation['metrics'].items():
        subset = [r for r in scored if r['method'] == method and r['scope'] == 'all'
                  and float(r['tolerance_um']) == MATCH_RADIUS_UM]
        assert len(subset) == 6
        for key in ['hand', 'predicted', 'matched', 'missed', 'extra']:
            assert sum(int(r[key]) for r in subset) == metrics[key]
    return dict(sections=sections, totals={k: sum(r[k] for r in sections)
                for k in ['full', 'alternative', 'shared', 'full_only', 'alternative_only']},
                neighborhoods=dict(circle_counts), matching_radius_um=MATCH_RADIUS_UM,
                comparison_is_accuracy=False,
                validation_note='Verifies published score arithmetic; does not rerun hand-count evaluation.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / '.local/recomputed_summary.json')
    args = parser.parse_args()
    result = analyze()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['totals'], indent=2))
    print('Verified all 16 sections, regional totals, neighborhood memberships and validation arithmetic.')
