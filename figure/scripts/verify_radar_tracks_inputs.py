"""Verify frozen radar/track plotting inputs without R, models, or third-party packages.

This verifies plotting-table identities and recomputes local correlations and
Jaccard from stored counts. It does not recompute global Pearson/AUROC from reads.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

SPECIES = ('arabidopsis', 'maize', 'wheat')
COUNTS = dict(arabidopsis=90132, maize=99572, wheat=197312)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def correlation(x, y):
    require(len(x) == len(y) and len(x) > 1, 'Invalid local correlation sample')
    mx, my = math.fsum(x) / len(x), math.fsum(y) / len(y)
    dx, dy = [v - mx for v in x], [v - my for v in y]
    return math.fsum(a*b for a, b in zip(dx, dy)) / math.sqrt(
        math.fsum(a*a for a in dx) * math.fsum(b*b for b in dy))


def indexed(rows):
    result = {(row['target_species'], row['source_model']): row for row in rows}
    require(len(result) == len(rows) == 9, 'Expected nine unique species/model pairs')
    require(set(result) == {(a, b) for a in SPECIES for b in SPECIES}, 'Pair coverage differs')
    return result


def verify(data_dir, rendered_dir=None):
    manifest = json.loads((data_dir / 'data_manifest.json').read_text(encoding='utf-8'))
    require(len(manifest['files']) == 8, 'Expected eight original data files')
    for entry in manifest['files']:
        path = (data_dir / entry['path']).resolve()
        require(data_dir.resolve() in path.parents, 'Input escapes the data directory')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'],
                'Input SHA-256 differs: ' + entry['path'])
    corr_rows = table(data_dir / 'correlation.csv')
    require(all(r['strand'] == 'pooled' for r in corr_rows), 'Pearson is not pooled')
    corr = indexed(corr_rows)
    auc = indexed([r for r in table(data_dir / 'auroc.csv') if r['strand'] == 'micro'])
    peak = indexed([r for r in table(data_dir / 'peak_overlap.csv') if float(r['percentile']) == 85])
    metrics = []
    for key in corr:
        c, a, p = corr[key], auc[key], peak[key]
        n = COUNTS[key[0]]
        require(int(c['n_observations']) == int(a['n_observations']) == int(p['n_strand_window_pairs']) == n,
                'Pooled sample sizes differ')
        require(int(c['n_windows']) * 2 == n, 'Strand-window pair count differs')
        tp, fp, fn, tn = [int(p[k]) for k in ('tp', 'fp', 'fn', 'tn')]
        require(tp + fp + fn + tn == n, 'Confusion counts do not cover the samples')
        require(tp == int(p['intersection_n']), 'Intersection differs from TP')
        require(tp + fp + fn == int(p['union_n']), 'Jaccard denominator differs')
        require(abs(tp / (tp + fp + fn) - float(p['micro_jaccard'])) < 1e-12, 'Jaccard differs')
        require(int(a['positive_n']) + int(a['negative_n']) == n, 'AUROC sample coverage differs')
        scores = (float(c['pearson_r']), float(a['HiLearn_AUROC']), float(p['micro_jaccard']))
        require(all(math.isfinite(v) and 0 <= v <= 1 for v in scores), 'Unexpected plotted metric range')
        metrics.append(dict(target_species=key[0], source_model=key[1],
                            pearson_r=scores[0], HiLearn_AUROC=scores[1], micro_jaccard=scores[2]))
    regions = table(data_dir / 'track_selected.csv')
    require(len(regions) == 3 and {r['target_species'] for r in regions} == set(SPECIES), 'Selected regions differ')
    local = []
    for region in regions:
        require(int(region['display_rank']) == 1, 'Unexpected selected region rank')
        rows = table(data_dir / 'regions' / (region['candidate_id'] + '.tracks.csv.gz'))
        require(len(rows) == int(region['n_total_windows']), 'Track row count differs')
        for row in rows:
            start = int(row['input_start0'])
            require(int(row['input_end0']) - start == 1024, 'Input width differs')
            require(int(row['center_start0']) == start + 448 and int(row['center_end0']) == start + 576,
                    'Displayed central 128-bp span differs')
        valid = [r for r in rows if r['valid'].lower() in ('1', 'true')]
        require(len(valid) == int(region['n_valid_windows']), 'Valid track count differs')
        strand = region['primary_strand']
        observed = [math.log2(1 + float(r['mean_cpm_1024_' + strand])) for r in valid]
        for source, field in ((region['target_species'], 'r_self'), (region['cross_model1'], 'r_model1'),
                              (region['cross_model2'], 'r_model2')):
            raw = [float(r['prediction_' + source + '_' + strand]) for r in valid]
            require(all(math.isfinite(v) and 0 <= v <= 1 for v in raw), 'Invalid prediction score')
            value = correlation(observed, [math.log2(1 + v) for v in raw])
            require(abs(value - float(region[field])) < 1e-12, 'Local correlation differs')
            local.append(dict(target_species=region['target_species'], source_model=source,
                              strand=strand, n_valid=len(valid), pearson_r=value))
    if rendered_dir is not None:
        rendered = indexed(table(rendered_dir / 'validation/radar_values.csv'))
        for row in metrics:
            actual = rendered[(row['target_species'], row['source_model'])]
            for key in ('pearson_r', 'HiLearn_AUROC', 'micro_jaccard'):
                require(abs(float(actual[key]) - row[key]) < 1e-12, 'Rendered metric differs')
        tracks = table(rendered_dir / 'validation/tracks_rendered.csv')
        require(len(tracks) == 9, 'Expected nine rendered track rows')
        lookup = {(r['target_species'], r['source_model']): r for r in local}
        for row in tracks:
            if row['track'] != 'experimental':
                require(abs(float(row['r']) - lookup[(row['target'], row['track'])]['pearson_r']) < 1e-12,
                        'Rendered local correlation differs')
    return dict(passed=True, input_files=8, plotted_metric_values=27,
                local_correlations_recomputed=local, selected_regions=3,
                rendered_tables_checked=rendered_dir is not None,
                scope='Frozen plotting inputs; global Pearson/AUROC read from original summary tables, not recomputed from raw reads.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data/radar_tracks')
    parser.add_argument('--rendered-dir', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = verify(args.data_dir.resolve(), args.rendered_dir)
    output = json.dumps(report, indent=2) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(output, encoding='utf-8')
    print(output, end='')


if __name__ == '__main__':
    main()
