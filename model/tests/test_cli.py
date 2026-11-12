import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from plant_pausing_bowtie2.cli import load_config, main


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.path = self.root / 'config.json'
        self.config = {'output': 'results', 'settings': {'reference_index': 'reference/genome'},
                       'samples': [{'run': 'sample1', 'reads': ['reads/a.fastq']}]}
        self.write_config()

    def tearDown(self):
        self.temp.cleanup()

    def write_config(self):
        self.path.write_text(json.dumps(self.config))

    def test_relative_paths_are_relative_to_config(self):
        output, [(sample, settings)] = load_config(self.path)
        self.assertEqual(output, (self.root / 'results').resolve())
        self.assertEqual(sample['reads'][0], str((self.root / 'reads/a.fastq').resolve()))
        self.assertEqual(settings['reference_index'], str((self.root / 'reference/genome').resolve()))

    def test_duplicate_run_rejected(self):
        self.config['samples'].append(dict(self.config['samples'][0]))
        self.write_config()
        with self.assertRaises(ValueError):
            load_config(self.path)

    def test_per_sample_protocol_does_not_change_other_sample(self):
        self.config['settings']['trimming'] = {'poly_a': True, 'random_clip': 0}
        self.config['samples'][0]['trimming'] = {'poly_a': False}
        self.config['samples'].append({'run': 'sample2', 'reads': ['other.fastq']})
        self.write_config()
        _, prepared = load_config(self.path)
        self.assertFalse(prepared[0][1]['trimming']['poly_a'])
        self.assertTrue(prepared[1][1]['trimming']['poly_a'])

    def test_plan_and_selection_do_not_create_outputs(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.assertEqual(main(['plan', str(self.path), '--sample', 'sample1']), 0)
        result = json.loads(stream.getvalue())
        self.assertEqual(len(result['samples']), 1)
        self.assertFalse((self.root / 'results').exists())
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(['plan', str(self.path), '--sample', 'missing']), 1)

    def test_index_preview_does_not_require_fasta_or_create_directory(self):
        output = self.root / 'not-created' / 'genome[1]'
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.assertEqual(main(['build-index', '--fasta', str(self.root/'missing.fa'),
                                   '--prefix', str(output), '--threads', '2', '--dry-run']), 0)
        self.assertEqual(json.loads(stream.getvalue())['command'][0], 'bowtie2-build')
        self.assertFalse(output.parent.exists())


if __name__ == '__main__':
    unittest.main()
