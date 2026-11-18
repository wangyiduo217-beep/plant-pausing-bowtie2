"""Check exact input/reference bytes for the three selected final layouts."""
from pathlib import Path
import hashlib,json

def main():
    root=Path(__file__).resolve().parents[2]
    manifest=json.loads((root/'figure/docs/selected_final_input_manifest.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        path=root/item['path']
        if path.stat().st_size!=item['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
            raise SystemExit('Input differs from selected final version: '+item['path'])
    print(json.dumps(dict(passed=True,files=len(manifest['files']),
        scope='Byte integrity only. Scientific and PDF-render comparisons are documented in each figure module.')))

if __name__=='__main__':main()
