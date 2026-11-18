"""Verify frozen source bytes and reproduce the 48 summaries from all windows."""
from pathlib import Path
import csv,gzip,hashlib,json
import numpy as np
from scipy.stats import pearsonr,spearmanr

def verify(figure,output=None):
    figure=Path(figure);root=figure.parent;data=figure/'data/figureS1'
    manifest=json.loads((data/'data_manifest.json').read_text())
    for item in manifest['inputs']:
        path=root/item['path']
        assert path.stat().st_size==item['bytes']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'],path
    with (data/'regional_correlation_and_mean_signal.tsv').open(encoding='utf-8') as f:
        records=list(csv.DictReader(f,delimiter='\t'))
    assert len(records)==48
    errors={};counts={};results=[]
    for species,n in [('arabidopsis',45066),('wheat',98656),('maize',49786)]:
        with gzip.open(figure/'data/regional'/f'{species}.regional_plot_data.tsv.gz','rt',encoding='utf-8') as f:
            windows=list(csv.DictReader(f,delimiter='\t'))
        assert len(windows)==n
        assert len({w['window_index'] for w in windows})==n
        counts[species]=n
        region=np.array([w['region'] for w in windows])
        for strand in ('plus','minus'):
            x=np.array([float(w['mean_cpm_'+strand]) for w in windows]);y=np.array([float(w['prediction_'+strand]) for w in windows])
            assert np.all(np.isfinite(x)) and np.all(x>=0)
            assert np.all(np.isfinite(y)) and np.all((y>=0)&(y<=1))
            subset=[r for r in records if r['species']==species and r['strand']==strand]
            assert len(subset)==8 and sum(int(r['n_windows']) for r in subset)==n
            for r in subset:
                take=region==r['region'];a=x[take];b=y[take];lx=np.log2(1+a);ly=np.log2(1+b)
                actual=dict(n_windows=len(a),pearson_r=pearsonr(lx,ly).statistic,
                    spearman_rho=spearmanr(lx,ly).statistic,mean_cpm=np.mean(a),sd_cpm=np.std(a,ddof=1),
                    median_cpm=np.median(a),zero_count=int(np.sum(a==0)),zero_fraction=np.mean(a==0),
                    mean_log2p1_cpm=np.mean(lx),log2p1_mean_cpm=np.log2(1+np.mean(a)))
                for k,v in actual.items():
                    err=abs(v-float(r[k]));errors[k]=max(errors.get(k,0),float(err))
                    assert np.isclose(v,float(r[k]),rtol=1e-12,atol=1e-12),(species,strand,r['region'],k,v,r[k])
                results.append(dict(species=species,strand=strand,region=r['region'],**actual))
    report=dict(passed=True,inputs=manifest['inputs'],rows_checked=48,total_test_windows=sum(counts.values()),
                window_counts=counts,max_absolute_error=errors,recomputed=results)
    if output:
        output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    figure=Path(__file__).resolve().parents[2]
    report=verify(figure,figure/'results/figureS1/validation/input_validation.json')
    print(json.dumps({k:v for k,v in report.items() if k not in ('inputs','recomputed')}))
