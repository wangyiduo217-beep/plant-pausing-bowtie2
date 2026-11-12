"""Configuration-driven command line; paths are relative to the JSON file."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

from .workflow import run_sample, _execute


def resolve_path(value, base):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Paths must be nonempty strings')
    path = Path(value).expanduser()
    return str((base / path).resolve() if not path.is_absolute() else path.resolve())


def load_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(config, dict) or not isinstance(config.get('settings'), dict):
        raise ValueError('Configuration must contain a settings object')
    if not isinstance(config.get('samples'), list) or not config['samples']:
        raise ValueError('Configuration must contain a nonempty samples list')
    settings = copy.deepcopy(config['settings'])
    for key in ('reference_index', 'rrna_index'):
        if settings.get(key):
            settings[key] = resolve_path(settings[key], path.parent)
    output = Path(resolve_path(config.get('output', '../results'), path.parent))
    runs = set()
    prepared = []
    for entry in config['samples']:
        if not isinstance(entry, dict):
            raise ValueError('Each sample must be an object')
        sample = copy.deepcopy(entry)
        run = sample.get('run')
        if not isinstance(run, str) or run in runs:
            raise ValueError('Every sample must have a distinct run identifier')
        runs.add(run)
        if 'reads' in sample:
            if not isinstance(sample['reads'], list):
                raise ValueError('reads must be a list of one or two FASTQ paths')
            sample['reads'] = [resolve_path(x, path.parent) for x in sample['reads']]
        if sample.get('sra'):
            sample['sra'] = resolve_path(sample['sra'], path.parent)
        individual = copy.deepcopy(settings)
        if 'trimming' in sample:
            overrides = sample.pop('trimming')
            if not isinstance(overrides, dict):
                raise ValueError('Per-sample trimming must be an object')
            individual['trimming'] = {**individual.get('trimming', {}), **overrides}
        prepared.append((sample, individual))
    return output, prepared


def index_command(fasta, prefix, threads):
    if threads < 1:
        raise ValueError('Index threads must be positive')
    return ['bowtie2-build', '--threads', str(threads), str(fasta), str(prefix)]


def build_index(args):
    fasta = Path(args.fasta).expanduser().resolve()
    prefix = Path(args.prefix).expanduser().resolve()
    command = index_command(fasta, prefix, args.threads)
    if args.dry_run:
        return {'command': command, 'note': 'Execution builds in a temporary sibling directory before installing index shards.'}
    if sys.platform != 'linux':
        raise RuntimeError('Execution requires Linux; use --dry-run to inspect commands elsewhere')
    if not fasta.is_file():
        raise ValueError('Reference FASTA does not exist: ' + str(fasta))
    if shutil.which('bowtie2-build') is None:
        raise RuntimeError('bowtie2-build is not on PATH')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    import fcntl
    with prefix.with_name(prefix.name + '.build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        existing = [p for p in prefix.parent.iterdir()
                    if p.name.startswith(prefix.name + '.') and p.name.endswith(('.bt2', '.bt2l'))]
        if existing or prefix.with_name(prefix.name + '.build.json').exists():
            raise FileExistsError('Index prefix already has output; use a new prefix')
        staging = prefix.parent / ('.' + prefix.name + '.building-' + uuid.uuid4().hex)
        staging.mkdir()
        temporary_prefix = staging / 'index'
        log = staging / 'build.log'
        executed_command = index_command(fasta, temporary_prefix, args.threads)
        def interrupted(signum, frame):
            raise RuntimeError('Index build interrupted by signal ' + str(signum))
        previous = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGTERM, signal.SIGINT)}
        try:
            _execute({'stage': 'build_index', 'argv': executed_command}, log)
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
        shard_names = ['1', '2', '3', '4', 'rev.1', 'rev.2']
        suffix = next((suffix for suffix in ('bt2', 'bt2l')
                       if all((staging / ('index.' + name + '.' + suffix)).is_file()
                              and (staging / ('index.' + name + '.' + suffix)).stat().st_size > 0
                              for name in shard_names)), None)
        if suffix is None:
            raise RuntimeError('Index build did not produce all six shards; inspect ' + str(log))
        outputs = []
        for name in shard_names:
            target = prefix.with_name(prefix.name + '.' + name + '.' + suffix)
            (staging / ('index.' + name + '.' + suffix)).rename(target)
            outputs.append(str(target))
        stat = fasta.stat()
        result = {'fasta': str(fasta), 'fasta_size': stat.st_size, 'fasta_mtime_ns': stat.st_mtime_ns,
                  'prefix': str(prefix), 'command': command, 'executed_command': executed_command, 'outputs': outputs,
                  'log': str(log), 'finished_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        record = prefix.with_name(prefix.name + '.build.json')
        temporary_record = record.with_suffix('.json.tmp')
        temporary_record.write_text(json.dumps(result, indent=2))
        temporary_record.replace(record)
        return result


def main(argv=None):
    parser = argparse.ArgumentParser(description='Plant GRO-seq trimming, Bowtie2 alignment and QC')
    subs = parser.add_subparsers(dest='command', required=True)
    for command in ('plan', 'run'):
        sub = subs.add_parser(command, help='Print commands without executing' if command == 'plan' else 'Process configured samples')
        sub.add_argument('config', help='JSON configuration file')
        sub.add_argument('--sample', action='append', help='Run identifier to select; repeat for multiple samples')
    index = subs.add_parser('build-index', help='Build a Bowtie2 index without overwriting an existing index')
    index.add_argument('--fasta', required=True)
    index.add_argument('--prefix', required=True)
    index.add_argument('--threads', type=int, default=8)
    index.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.command == 'build-index':
            print(json.dumps(build_index(args), indent=2, ensure_ascii=False))
            return 0
        output, prepared = load_config(args.config)
        if args.sample:
            requested = set(args.sample)
            unknown = requested - {s['run'] for s, _ in prepared}
            if unknown:
                raise ValueError('Unknown sample identifiers: ' + ', '.join(sorted(unknown)))
            prepared = [(s, c) for s, c in prepared if s['run'] in requested]
        # Validate every selected command plan before the first sample writes anything.
        plans = [run_sample(sample, settings, output, dry_run=True) for sample, settings in prepared]
        if args.command == 'plan':
            print(json.dumps({'output': str(output), 'samples': plans}, indent=2, ensure_ascii=False, default=str))
            return 0
        for sample, settings in prepared:
            result = run_sample(sample, settings, output)
            print(json.dumps(result, ensure_ascii=False, default=str), flush=True)
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError, KeyError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
