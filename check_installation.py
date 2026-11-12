#!/usr/bin/env python3
"""Run the distributed-model and figure acceptance checks; no training is started."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable or absolute path')
    parser.add_argument('--font', default='Arial', help='An installed font; the independent Linux check uses DejaVu Sans')
    parser.add_argument('--output-dir', type=Path, help='Default: <root>/reproducibility/acceptance')
    parser.add_argument('--require-isolated-env', action='store_true', help='Require a venv with system-site-packages disabled')
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root/'reproducibility/acceptance').resolve()
    output.mkdir(parents=True, exist_ok=True)
    cfg = Path(sys.prefix)/'pyvenv.cfg'
    config_text = cfg.read_text(encoding='utf-8').lower() if cfg.is_file() else ''
    isolated = sys.prefix != sys.base_prefix and 'include-system-site-packages = false' in config_text
    report = {'schema_version': 1, 'scope': 'Published checkpoint inference and all supplied figure workflows; excludes raw-read processing and retraining.',
              'python': platform.python_version(), 'platform': platform.system(),
              'venv_without_system_site_packages': isolated, 'font_family': args.font,
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
    rscript = shutil.which(args.rscript)
    if rscript is None:
        return fail('Rscript was not found. Install R with Cairo, or provide --rscript /path/to/Rscript.')
    env = os.environ.copy()
    env['PATH'] = str(Path(rscript).resolve().parent) + os.pathsep + env.get('PATH', '')
    env['PYTHONPATH'] = str(root/'model/src')
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['CUDA_VISIBLE_DEVICES'] = ''
    env['MPLBACKEND'] = 'Agg'
    env['GROSEQ_FONT_FAMILY'] = args.font
    env['OMP_NUM_THREADS'] = '2'
    env['MKL_NUM_THREADS'] = '2'
    python = sys.executable
    commands = [
        ('dependency_consistency', [python,'-m','pip','check']),
        ('software_tests', [python,'-m','unittest','discover','-s',str(root/'model/tests'),'-v']),
        ('complete_publication_inputs', [python,str(root/'figure/scripts/verify_publication_inputs.py')]),
        ('decile_and_frequency_inputs', [python,str(root/'figure/scripts/verify_decile_inputs.py')]),
        ('pretrained_downloads_and_checksums', [python,str(root/'model/scripts/download_pretrained.py'),'--species','all']),
        ('all_real_dna_cpu_predictions', [python,str(root/'model/scripts/check_pretrained_example.py'),'--species','all','--device','cpu','--threads','2','--force']),
        ('R_cairo_available', [rscript,'--vanilla','-e','stopifnot(capabilities("cairo")); cat(R.version.string, " Cairo=TRUE\\n")']),
        ('six_R_figure_families', [rscript,'--vanilla',str(root/'figure/plotting/run_all.R')]),
        ('HiLearn_ROC_PR', [python,str(root/'figure/scripts/plot_publication_roc_pr.py'),'--font-family',args.font]),
    ]
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
    report['generated_files'] = sorted(p.relative_to(root).as_posix() for p in (root/'figure/results').rglob('*') if p.is_file())
    report['seconds'] = round(time.monotonic()-started_all,3)
    report['passed'] = True
    save()
    print('PASS: all acceptance checks completed; 72 real DNA sequences / 144 strand scores checked.')
    print('Summary: '+str(summary))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
