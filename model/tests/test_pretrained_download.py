"""Integrity checks for interrupted and range-capable checkpoint downloads."""
import hashlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from plant_pausing_bowtie2.pretrained import PretrainedError, checkpoint_path, download_checkpoint


class Response(io.BytesIO):
    def __init__(self, body, content_range=None, status=206):
        super().__init__(body)
        self.status=status
        self.headers={'Content-Range':content_range}
    def geturl(self):
        return 'https://example.org/cdn/model.pt'


class RangeDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root=Path(self.temporary.name)
        self.body=b'abcdefgh'
        self.item={'url':'https://example.org/model.pt','bytes':len(self.body),
                   'sha256':hashlib.sha256(self.body).hexdigest(),
                   'checkpoint_relative_path':'results/model/best.pt'}
        self.chunk=patch('plant_pausing_bowtie2.pretrained.DOWNLOAD_CHUNK_BYTES',4)
        self.chunk.start(); self.addCleanup(self.chunk.stop)

    def assert_clean_failure(self, responses):
        target=checkpoint_path(self.root,self.item)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(b'previous')
        with patch('urllib.request.urlopen',side_effect=responses), patch('time.sleep'):
            with self.assertRaises(PretrainedError): download_checkpoint(self.root,self.item,force=True)
        self.assertEqual(target.read_bytes(),b'previous')
        self.assertEqual(list(target.parent.glob('*.partial')),[])

    def test_exact_ranges_are_assembled_and_hash_checked(self):
        with patch('urllib.request.urlopen',side_effect=[Response(b'abcd','bytes 0-3/8'),Response(b'efgh','bytes 4-7/8')]) as request:
            result=download_checkpoint(self.root,self.item)
        self.assertEqual(result['transport'],'verified_http_ranges')
        self.assertEqual(checkpoint_path(self.root,self.item).read_bytes(),self.body)
        self.assertEqual([call[0][0].get_header('Range') for call in request.call_args_list],['bytes=0-3','bytes=4-7'])

    def test_truncated_range_retries_same_offset_without_duplicate_bytes(self):
        replies=[Response(b'ab','bytes 0-3/8'),Response(b'abcd','bytes 0-3/8'),Response(b'efgh','bytes 4-7/8')]
        with patch('urllib.request.urlopen',side_effect=replies) as request, patch('time.sleep'):
            download_checkpoint(self.root,self.item)
        self.assertEqual(checkpoint_path(self.root,self.item).read_bytes(),self.body)
        self.assertEqual([call[0][0].get_header('Range') for call in request.call_args_list],['bytes=0-3','bytes=0-3','bytes=4-7'])

    def test_unsupported_range_at_start_accepts_only_complete_full_response(self):
        with patch('urllib.request.urlopen',return_value=Response(self.body,status=200)):
            result=download_checkpoint(self.root,self.item)
        self.assertEqual(result['transport'],'full_response')

    def test_wrong_content_range_is_rejected(self):
        self.assert_clean_failure([Response(b'abcd','bytes 1-4/8')])

    def test_full_response_after_a_range_is_not_appended(self):
        self.assert_clean_failure([Response(b'abcd','bytes 0-3/8'),Response(self.body,status=200)])

    def test_repeated_network_failure_is_bounded(self):
        with patch('urllib.request.urlopen',side_effect=TimeoutError('network timeout')) as request, patch('time.sleep'):
            with self.assertRaises(PretrainedError): download_checkpoint(self.root,self.item)
        self.assertEqual(request.call_count,3)
        self.assertFalse(checkpoint_path(self.root,self.item).exists())

    def test_corrupt_ranges_never_replace_existing_checkpoint(self):
        self.assert_clean_failure([Response(b'xxxx','bytes 0-3/8'),Response(b'efgh','bytes 4-7/8')])


if __name__=='__main__': unittest.main()
