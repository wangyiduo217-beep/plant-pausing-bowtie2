#!/usr/bin/env python3
"""Reproduce the four selected final figures and panels; no model download or inference."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


ENVIRONMENT_PROBE = """
import importlib.metadata as m, json, os, platform, sys
names = ['numpy','scipy','matplotlib','reportlab','pypdf','pdfplumber','PyMuPDF']
versions = {name:m.version(name) for name in names}
assert versions['pypdf'].split('.')[0] == '6', 'Final figures require pypdf >=6,<7'
assert versions['reportlab'].split('.')[0] == '4', 'Final figures require reportlab >=4,<5'
print('__GROSEQ_REPORT__'+json.dumps(dict(executable=sys.executable,python=platform.python_version(),packages=versions,
      pythonpath=os.environ.get('PYTHONPATH',''),venv=sys.prefix != sys.base_prefix)))
"""

PDF_CHECK = """
import json, pathlib, sys, pymupdf as fitz
records = []
for filename in sys.argv[1:]:
    path = pathlib.Path(filename)
    assert path.is_file() and path.stat().st_size > 0, str(path)
    with fitz.open(path) as doc:
        assert len(doc) == 1, str(path)
        page = doc[0]
        mm = [page.rect.width * 25.4 / 72, page.rect.height * 25.4 / 72]
        assert abs(mm[0]-210) < .1 and abs(mm[1]-297) < .1, (str(path),mm)
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        assert spans, 'No extractable text: '+str(path)
        fonts = sorted({s['font'] for s in spans})
        assert all('Arial' in f for f in fonts), (str(path),fonts)
        records.append(dict(path=str(path),bytes=path.stat().st_size,pages=1,page_mm=mm,fonts=fonts))
print('__GROSEQ_REPORT__'+json.dumps(dict(passed=True,documents=records,scope='A4 page/text/font checks; not visual or pixel-equivalence certification.')))
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', default=sys.executable, help='Python executable with the final-figure dependencies')
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable; R with Cairo is required for Figure 2')
    parser.add_argument('--font-dir', type=Path, default=Path(os.environ.get('GROSEQ_FONT_DIR', 'C:/Windows/Fonts')),
                        help='Licensed Arial TTF directory; Arial must also be installed for R/Cairo')
    parser.add_argument('--output-dir', type=Path, help='Parent directory for figure1/figure2/figureS1/figureS2; default figure/results')
    parser.add_argument('--report-dir', type=Path, help='Default: <output-dir>/reproduction_validation')
    parser.add_argument('--dry-run', action='store_true', help='Print commands without executing or creating files')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    figure = root/'figure'
    plotting = figure/'plotting'
    output = (args.output_dir or figure/'results').resolve()
    report_dir = (args.report_dir or output/'reproduction_validation').resolve()
    font_dir = args.font_dir.resolve()
    python = shutil.which(args.python) or args.python
    rscript = shutil.which(args.rscript) or args.rscript
    out1, out2, outs1, outs2 = [output/name for name in ('figure1','figure2','figureS1','figureS2')]
    source = out1/'source_panels'
    font_args = ['--font-dir', str(font_dir)]
    a4s = [out1/'GROseq_combined_v2_A4.pdf',
           out2/'figures/GROseq_radar_tracks_with_legend_A4.pdf',
           outs1/'region_signal_six_panels_A4.pdf',
           outs2/'GROseq_CPM_deciles_grouped_strands_A4.pdf']
    commands = [
        ('figure_environment', [python,'-c',ENVIRONMENT_PROBE]),
        ('figure_dependency_consistency', [python,'-m','pip','check']),
        ('selected_input_identities', [python,str(figure/'scripts/verify_selected_final_inputs.py')]),
        ('figureS2_decile_inputs', [python,str(figure/'scripts/verify_decile_inputs.py')]),
        ('R_cairo', [rscript,'--vanilla','-e','stopifnot(capabilities("cairo")); cat(R.version.string)']),
        ('figure1_source_panels', [python,str(plotting/'figure1/build_source_panels.py'),'--output-dir',str(source)]+font_args),
        ('figure1_composite', [python,str(plotting/'figure1/build_combined_v2.py'),
                               '--source-figure2',str(source/'GROseq_v1_tracks_OLS_correlations_A4.pdf'),
                               '--output-dir',str(out1)]+font_args),
        ('figure2_composite', [rscript,'--vanilla',str(plotting/'figure2/plot_radar_tracks_composite.R'),'--output-dir',str(out2)]),
        ('figure2_components', [rscript,'--vanilla',str(plotting/'figure2/plot_components.R'),'--component','all','--output-dir',str(out2)]),
        ('figure2_rendered_values', [python,str(figure/'scripts/verify_radar_tracks_inputs.py'),
                                    '--rendered-dir',str(out2),'--report',str(report_dir/'figure2_values.json')]),
        ('figureS1_regional', [python,str(plotting/'figureS1/build_region_signal.py'),'--output-dir',str(outs1)]+font_args),
        ('figureS2_deciles', [python,str(plotting/'figureS2/run_deciles_final.py'),'--panel','all','--output-dir',str(outs2)]+font_args),
        ('final_A4_structure', [python,'-c',PDF_CHECK]+[str(p) for p in a4s]),
    ]
    if args.dry_run:
        print(json.dumps([{'name':name,'command':command} for name,command in commands],indent=2))
        return 0
    for executable in (python,rscript):
        if not shutil.which(executable) and not Path(executable).is_file():
            parser.error('Executable not found: '+executable)
    for name in ('arial.ttf','arialbd.ttf','ariali.ttf','arialbi.ttf'):
        if not (font_dir/name).is_file():
            parser.error('Missing licensed Arial font: '+str(font_dir/name))
    report_dir.mkdir(parents=True,exist_ok=True)
    summary = report_dir/'summary.json'
    report = dict(schema_version=1,passed=False,scope='Four selected final figures and component exports; no model inference or downloads.',
                  root=str(root),output=str(output),font_directory=str(font_dir),steps=[])
    def save():
        summary.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    env = os.environ.copy()
    env.update(CUDA_VISIBLE_DEVICES='',MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1',
               OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',GROSEQ_FONT_DIR=str(font_dir))
    env['PATH'] = str(Path(rscript).resolve().parent)+os.pathsep+env.get('PATH','')
    started_all = time.monotonic()
    for name, command in commands:
        print('RUN: '+name,flush=True)
        started = time.monotonic()
        result = subprocess.run(command,cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                text=True,encoding='utf-8',errors='replace')
        log = report_dir/(name+'.log')
        log.write_text(result.stdout,encoding='utf-8')
        report['steps'].append(dict(name=name,command=command,returncode=result.returncode,
                                   passed=result.returncode==0,seconds=round(time.monotonic()-started,3),log=str(log)))
        save()
        if result.returncode:
            print('FAIL: '+name+'; see '+str(log),file=sys.stderr)
            return result.returncode
        if name in ('figure_environment','final_A4_structure'):
            payloads = [line.removeprefix('__GROSEQ_REPORT__') for line in result.stdout.splitlines()
                        if line.startswith('__GROSEQ_REPORT__')]
            if len(payloads) != 1:
                report['error'] = 'Expected exactly one structured report from '+name
                save()
                print('FAIL: '+report['error']+'; see '+str(log),file=sys.stderr)
                return 1
            report[name] = json.loads(payloads[0])
        print('PASS: '+name,flush=True)
    report['seconds'] = round(time.monotonic()-started_all,3)
    report['passed'] = True
    save()
    print('PASS: Figure 1, Figure 2, Figure S1 and Figure S2 reproduced; '+str(summary))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
