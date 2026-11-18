#!/usr/bin/env python3
"""Check pretrained CPU inference and figure inputs; optionally reproduce four final figures."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable or absolute path')
    parser.add_argument('--font', default='Arial', help='Legacy option; final figures require Arial when --figures is used')
    parser.add_argument('--figures', action='store_true', help='Also reproduce Figure 1, Figure 2, Figure S1 and Figure S2')
    parser.add_argument('--figure-python', default=sys.executable, help='Python executable for the separate final-figure environment')
    parser.add_argument('--font-dir', type=Path, help='Directory containing licensed Arial TTF files for final figures')
    parser.add_argument('--output-dir', type=Path, help='Default: <root>/reproducibility/acceptance')
    parser.add_argument('--require-isolated-env', action='store_true', help='Require a venv with system-site-packages disabled')
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root/'reproducibility/acceptance').resolve()
    output.mkdir(parents=True, exist_ok=True)
    cfg = Path(sys.prefix)/'pyvenv.cfg'
    config_text = cfg.read_text(encoding='utf-8').lower() if cfg.is_file() else ''
    isolated = sys.prefix != sys.base_prefix and 'include-system-site-packages = false' in config_text
    report = {'schema_version': 2, 'scope': 'Published checkpoint CPU inference and bundled figure-input verification; final-figure reproduction is opt-in. Excludes raw-read processing and retraining.',
              'python': platform.python_version(), 'platform': platform.system(),
              'venv_without_system_site_packages': isolated, 'figures_requested': args.figures,
              'font_family': 'Arial' if args.figures else None,
              'packages': {}, 'steps': [], 'passed': False}
    for name in ['plant-pausing-bowtie2','pip','setuptools','torch','numpy','scipy','matplotlib','hilearn','reportlab','pypdf']:
        try:
            report['packages'][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            report['packages'][name] = None
    summary = output/'summary.json'

    def save():
        summary.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')

    def fail(message):
        report['error'] = message
        save()
        print('FAIL: '+message, file=sys.stderr)
        return 1

    if args.require_isolated_env and not isolated:
        return fail('Use a virtual environment created without --system-site-packages.')
    if args.figures and args.font != 'Arial':
        return fail('The four final figures require Arial; provide licensed Arial via --font-dir and install it for R/Cairo.')
    env = os.environ.copy()
    env['PYTHONPATH'] = str(root/'model/src')
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['CUDA_VISIBLE_DEVICES'] = ''
    env['MPLBACKEND'] = 'Agg'
    env['OMP_NUM_THREADS'] = '2'
    env['MKL_NUM_THREADS'] = '2'
    python = sys.executable
    commands = [
        ('dependency_consistency', [python,'-m','pip','check']),
        ('software_tests', [python,'-m','unittest','discover','-s',str(root/'model/tests'),'-v']),
        ('complete_publication_inputs', [python,str(root/'figure/scripts/verify_publication_inputs.py')]),
        ('selected_final_input_identities', [python,str(root/'figure/scripts/verify_selected_final_inputs.py')]),
        ('figure2_radar_track_inputs', [python,str(root/'figure/scripts/verify_radar_tracks_inputs.py')]),
        ('figureS1_regional_inputs', [python,str(root/'figure/plotting/figureS1/verify_region_inputs.py')]),
        ('decile_inputs', [python,str(root/'figure/scripts/verify_decile_inputs.py')]),
        ('pretrained_downloads_and_checksums', [python,str(root/'model/scripts/download_pretrained.py'),'--species','all']),
        ('all_real_dna_cpu_predictions', [python,str(root/'model/scripts/check_pretrained_example.py'),'--species','all','--device','cpu','--threads','2','--force']),
    ]
    if args.figures:
        figure_command = [python, str(root/'figure/plotting/reproduce_all.py'),
                          '--python', args.figure_python, '--rscript', args.rscript,
                          '--report-dir', str(output/'final_figures')]
        if args.font_dir:
            figure_command += ['--font-dir', str(args.font_dir.resolve())]
        commands.append(('four_selected_final_figures', figure_command))
    started_all = time.monotonic()
    for name, command in commands:
        print('RUN: '+name, flush=True)
        started = time.monotonic()
        result = subprocess.run(command, cwd=root, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
        log = output/(name+'.log')
        log.write_text(result.stdout, encoding='utf-8')
        report['steps'].append({'name': name, 'returncode': result.returncode,
                                'seconds': round(time.monotonic()-started,3),
                                'log': log.name, 'passed': result.returncode == 0})
        save()
        if result.returncode:
            return fail(name+' failed. See '+log.name+' for details.')
        print('PASS: '+name, flush=True)
    report['inference'] = {}
    for species in ['arabidopsis_thaliana','triticum_aestivum','zea_mays']:
        item = json.loads((root/'model/results/pretrained_example'/(species+'.check.json')).read_text(encoding='utf-8'))
        report['inference'][species] = item['comparison']
        if not item['comparison']['passed']:
            return fail('CPU comparison did not pass for '+species)
    if args.figures:
        figure_summary = json.loads((output/'final_figures/summary.json').read_text(encoding='utf-8'))
        if not figure_summary['passed']:
            return fail('Final figure reproduction did not pass.')
        report['final_figures'] = figure_summary
    else:
        report['final_figures'] = {'requested': False, 'status': 'not run'}
    report['seconds'] = round(time.monotonic()-started_all,3)
    report['passed'] = True
    save()
    print('PASS: model and input checks completed; 72 real DNA sequences / 144 strand scores checked.')
    print('Final figures: '+('all four reproduced and structurally checked.' if args.figures else 'not requested; use --figures to reproduce them.'))
    print('Summary: '+str(summary))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
