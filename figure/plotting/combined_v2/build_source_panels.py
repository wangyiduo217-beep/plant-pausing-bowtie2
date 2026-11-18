"""Both axes log2(1+value): full-test grayscale scatter, zero-origin axes, OLS straight lines."""
from pathlib import Path
import argparse,csv,hashlib,io,json,math,shutil,subprocess,sys,textwrap

# Portability changes only: input/output paths and discovery of licensed Arial.
# Scientific transforms, density estimation, fitting, and plotting match the
# archived 2026-11-10 original in provenance/original_build_final.py.
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parents[2]/'data/combined_v2',help='Directory containing test_windows/ and tracks/ original data.')
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[2]/'results/combined_v2_source',help='Output directory; input data are never modified.')
parser.add_argument('--font-dir',type=Path,default=None,help='Directory containing licensed Arial TTFs (Arial regular, bold, italic, bold italic). Otherwise use installed Arial.')
ARGS=parser.parse_args()

import numpy as np
from scipy.stats import pearsonr,spearmanr,linregress
from scipy.ndimage import gaussian_filter,map_coordinates
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator,FuncFormatter
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
import pdfplumber

DATA=ARGS.data_dir.resolve()
OUT=ARGS.output_dir.resolve()
CORR=DATA/'test_windows'
TRACK=DATA/'tracks'
MM=72/25.4
W,H=183.,170.
if ARGS.font_dir is not None:
    if not ARGS.font_dir.is_dir():
        raise FileNotFoundError(f'Font directory does not exist: {ARGS.font_dir}')
    for font_path in ARGS.font_dir.glob('*.ttf'):
        fm.fontManager.addfont(str(font_path))
    for font_path in ARGS.font_dir.glob('*.TTF'):
        fm.fontManager.addfont(str(font_path))
FONT_PATHS={}
for name,weight,style in [('Arial','normal','normal'),('Arial-Bold','bold','normal'),('Arial-Italic','normal','italic'),('Arial-BoldItalic','bold','italic')]:
    prop=fm.FontProperties(family='Arial',weight=weight,style=style)
    try:
        font_path=fm.findfont(prop,fallback_to_default=False)
    except ValueError as exc:
        raise RuntimeError('Exact reproduction requires licensed Arial fonts. Install Arial or supply --font-dir; fonts are not redistributed by this repository.') from exc
    fm.fontManager.addfont(font_path)
    FONT_PATHS[name]=font_path
    pdfmetrics.registerFont(TTFont(name,font_path))
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':7,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'figure.titlesize':7,
    'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,
    'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','axes.formatter.useoffset':False,'savefig.facecolor':'white'})
SPECS=[('arabidopsis','Arabidopsis thaliana',14,45066,'3'),
       ('wheat','Triticum aestivum',5,98656,'3A / 3B / 3D'),
       ('maize','Zea mays',3,49786,'6')]
COLORS=['#21844B','#78B799','#3278A5','#8BAFD2']
sources={};checks=[];artists=[];font_spans=[]
Y_LIM=(0.,1.)
HIST_BINS=(256,256)
GREY=LinearSegmentedColormap.from_list('PublicationGrey',['#c4c4c4','#929292','#4b4b4b','#111111'])
FIT_COLOR='#555555'
METRIC_COLOR='#D55E00'
PANEL_CACHE={}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source(path):
    path=Path(path).resolve();sources[str(path)]=sha(path);return path
def load(path):
    with np.load(source(path),allow_pickle=False) as z:return {n:z[n] for n in z.files}
def dump(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def put(fig,x,y,text,bold=False,italic=False,size=7,ha='left',va='baseline',**kw):
    return fig.text(x/W,y/H,text,ha=ha,va=va,fontsize=size,fontfamily='Arial',
        fontweight='bold' if bold else 'normal',fontstyle='italic' if italic else 'normal',**kw)
def axis(fig,x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
def clean_axis(ax):
    ax.spines[['top','right']].set_visible(False)
    for spine in ax.spines.values():spine.set_linewidth(.5)
    ax.tick_params(pad=1.3,length=2,width=.5,labelsize=7)
    ax.grid(False)
def edges_inside(ax):
    labels=ax.get_xticklabels()
    if labels:labels[0].set_ha('left');labels[-1].set_ha('right')
def tick_y_inside(ax):
    labels=ax.get_yticklabels()
    if labels:labels[0].set_va('bottom');labels[-1].set_va('top')

def experimental_x(data,index):
    values=np.log2(1.+data['mean_cpm'][:,index].astype(float))
    # The cached x was computed on another runtime; log2 can differ by one ulp.
    np.testing.assert_allclose(values,data['log2_1plus_mean_cpm'][:,index],atol=1e-14,rtol=0)
    return values

def prepare_panel(xx,raw,key,strand):
    token=(key,strand)
    if token in PANEL_CACHE:return PANEL_CACHE[token]
    assert len(xx)==len(raw) and np.isfinite(xx).all() and np.isfinite(raw).all()
    assert np.all((xx>=0)&(xx<=12)&(raw>=0)&(raw<=1))
    yy=np.log2(1.+raw.astype(float))
    assert np.all((yy>=Y_LIM[0])&(yy<=Y_LIM[1]))
    # A binned Gaussian density estimate controls point color only; coordinates stay exact.
    hist,xe,ye=np.histogram2d(xx,yy,bins=HIST_BINS,range=((0,12),Y_LIM))
    assert int(hist.sum())==len(xx)
    smoothed=gaussian_filter(hist,sigma=(2.,2.),mode='reflect')
    coords=np.vstack(((xx-xe[0])/(xe[1]-xe[0])-.5,(yy-ye[0])/(ye[1]-ye[0])-.5))
    density=map_coordinates(smoothed,coords,order=1,mode='nearest',prefilter=False)
    density_color=np.log1p(density)/np.log1p(float(density.max()))
    order=np.argsort(density_color,kind='stable')
    base=linregress(xx,yy)
    lx=np.linspace(float(xx.min()),float(xx.max()),401)
    ly=base.intercept+base.slope*lx
    info={'n':len(xx),'n_drawn':len(xx),'y_transform':'log2(1 + raw predicted score)',
        'x_transform':'log2(1 + equal-weight mean per-SRX CPM)',
        'pearson_r':float(pearsonr(xx,yy).statistic),'spearman_rho':float(spearmanr(xx,yy).statistic),
        'raw_pearson_r':float(pearsonr(xx,raw).statistic),'raw_spearman_rho':float(spearmanr(xx,raw).statistic),
        'all_input_windows_drawn':True,'zero_scores':int(np.sum(raw==0)),
        'min_log_score':float(yy.min()),'max_log_score':float(yy.max()),
        'density_method':'256x256 all-point histogram, Gaussian sigma=2 bins per axis, bilinear evaluation; log1p relative color per panel',
        'fit_method':'Ordinary least squares on all current transformed coordinates',
        'fit_drawn_x_range':[float(lx.min()),float(lx.max())],
        'OLS_slope':float(base.slope),'OLS_intercept':float(base.intercept),
        'axis_x_limits':[0,12.2],'axis_y_limits':[0,1.02],
        'both_axis_minima_zero':True,
        'identity_line':False,'fit_line_color':FIT_COLOR,'fit_line_style':'dashed'}
    np.savez_compressed(OUT/'data'/f'{key}.{strand}.log_scatter.npz',
        x=xx,y=yy,raw_prediction=raw,density=density,density_color=density_color,draw_order=order,
        histogram=hist,hist_x_edges=xe,hist_y_edges=ye,fit_x=lx,fit_y=ly,
        OLS_slope=base.slope,OLS_intercept=base.intercept)
    PANEL_CACHE[token]=(yy,density_color,order,lx,ly,info)
    return PANEL_CACHE[token]

def scatter_panel(ax,xx,raw,key,strand,layout):
    yy,col,order,lx,ly,info=prepare_panel(xx,raw,key,strand)
    dot=ax.scatter(xx[order],yy[order],c=col[order],cmap=GREY,norm=Normalize(0,1),
        s=.5,alpha=.8,edgecolors='none',rasterized=True,zorder=2)
    np.testing.assert_array_equal(np.asarray(dot.get_offsets()),np.column_stack((xx[order],yy[order])))
    ax.plot(lx,ly,color=FIT_COLOR,ls='--',lw=.8,zorder=3)
    ax.set_xlim(0,12.2);ax.set_ylim(0,1.02)
    ax.set_xticks([0,6,12]);ax.set_yticks([0,.5,1],labels=['0','0.5','1'])
    clean_axis(ax);edges_inside(ax)
    return dict(info,layout=layout)

def count_legend(fig,x=148,y=3.5,w=34):
    cax=axis(fig,x,y,w,1.6)
    cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=GREY),cax=cax,orientation='horizontal')
    cb.set_ticks([0,1],labels=['Low','High'])
    cb.ax.tick_params(labelsize=7,length=0,pad=1)
    cb.ax.xaxis.set_minor_locator(FixedLocator([]))
    cb.outline.set_linewidth(.5)
    labs=cb.ax.get_xticklabels()
    for lab in labs:lab.set_fontfamily('Arial')
    labs[0].set_ha('left');labs[-1].set_ha('right')

def make_a4_proof(mainpdf,target,title,caption):
    packet=io.BytesIO();can=canvas.Canvas(packet,pagesize=(210*MM,297*MM))
    can.setTitle(title);can.setAuthor('Research figure')
    x=13.5*MM;y=96*MM;line_width=183*MM
    can.setFont('Arial-Bold',7);can.drawString(x,y,title);y-=11
    can.setFont('Arial',7);line='';lines=[]
    for word in caption.split():
        candidate=(line+' '+word).strip()
        if pdfmetrics.stringWidth(candidate,'Arial',7)>line_width:
            lines.append(line);line=word
        else:line=candidate
    if line:lines.append(line)
    for line in lines:can.drawString(x,y,line);y-=10
    assert y>20*MM
    can.showPage();can.save();packet.seek(0)
    page=PdfReader(packet).pages[0]
    page.merge_transformed_page(PdfReader(mainpdf).pages[0],
        Transformation().translate(tx=13.5*MM,ty=106*MM),expand=False)
    writer=PdfWriter();writer.add_page(page)
    writer.add_metadata({'/Title':title,'/Subject':'Original 183x170mm figure translated onto A4 without scaling; 7pt Arial, 8pt panel letters'})
    with target.open('wb') as f:writer.write(f)

def standalone(stats):
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    put(fig,18,169,'GRO-seq signal and sequence-model predictions',bold=True,va='top')
    results=[]
    for col,(key,name,n_srx,n,chroms) in enumerate(SPECS):
        x=18+58*col;middle=x+23.5
        put(fig,middle,160,name,italic=True,ha='center')
        put(fig,middle,155,f'chr{chroms}; n = {n:,}',ha='center')
        original=load(CORR/f'{key}.all_test_windows.npz')
        for si,strand in enumerate(['plus','minus']):
            top=148-si*66
            put(fig,x-7,top,'abcdef'[si*3+col],bold=True,size=8)
            put(fig,middle,top,'Plus strand' if si==0 else 'Minus strand',ha='center')
            xx=experimental_x(original,si);raw=original['predicted'][:,si].astype(float)
            yy=np.log2(1.+raw)
            r=float(pearsonr(xx,yy).statistic);rho=float(spearmanr(xx,yy).statistic)
            ref=stats[(key,strand)]
            assert abs(float(pearsonr(xx,raw).statistic)-float(ref['pearson_r']))<1e-12
            assert abs(rho-float(ref['spearman_rho']))<1e-12
            put(fig,middle,top-4,f'r = {r:.3f}; ρ = {rho:.3f}',ha='center',color=METRIC_COLOR)
            ax=axis(fig,x,top-53,47,46)
            info=scatter_panel(ax,xx,raw,key,strand,'standalone')
            if col==0:ax.set_ylabel('log2(1 + predicted score)',fontsize=7,labelpad=1.5)
            results.append({'species':key,'strand':strand,'pearson_r':r,'spearman_rho':rho,**info})
    put(fig,99.5,18,'GRO-seq signal (log2[1 + mean CPM])',ha='center')
    fig.add_artist(plt.Line2D([18/W,25/W],[9/H,9/H],transform=fig.transFigure,ls='--',color=FIT_COLOR,lw=.8))
    put(fig,27,9,'OLS fit',va='center')
    put(fig,111,9,'Within-panel density',ha='right',va='center')
    count_legend(fig,x=116,y=8.2,w=65)
    put(fig,18,2.5,'strand_seiplant_v1; repeats retained; Arabidopsis / wheat / maize: 14 / 5 / 3 SRX')
    fig.canvas.draw()
    for t in fig.findobj(match=Text):
        if t.get_visible() and t.get_text():
            assert t.get_fontfamily()==['Arial'] and t.get_fontsize() in [7,8]
    path=OUT/'figures/GROseq_v1_OLS_correlations_183x170mm.pdf'
    for ext in ['pdf','svg','png']:fig.savefig(path.with_suffix('.'+ext),dpi=600)
    plt.close(fig)
    caption=(
        '(a-c) Plus-strand and (d-f) minus-strand associations between GRO-seq signal and predicted scores '
        'for Arabidopsis thaliana, Triticum aestivum and Zea mays, respectively. '
        'All formal held-out 1,024-bp test windows are used: 45,066, 98,656 and 49,786 per species. '
        'Experimental signal is log2(1 + mean SRX CPM), with each SRX normalized before equal-weight averaging; '
        '14, 5 and 3 SRX are used, respectively. The y-axis is log2(1 + original strand_seiplant_v1 '
        'window score). Both axes use log2(1 + value); the experimental mean CPM is transformed once. Repeats are retained. '
        'Every point retains its original window identity; grayscale represents a binned Gaussian '
        'density estimate scaled within each panel. Dark-gray dashed straight lines are ordinary least-squares '
        'fits to all current transformed coordinates, drawn over the observed x range. '
        'Both axes start at zero; point coordinates are unchanged. '
        'r is Pearson correlation recalculated on the displayed transformed coordinates; rho is Spearman '
        'correlation. These are descriptive test-set associations. '
        'The GRO-seq measurements also supplied the training labels.'
    )
    make_a4_proof(path,OUT/'GROseq_v1_OLS_correlations_A4.pdf',
        'Sequence-based GRO-seq predictions: all held-out test windows.',caption)
    (OUT/'correlations_figure_legend_EN.txt').write_text(caption+'\n',encoding='utf-8')
    return results

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for folder in ['figures','validation','provenance','data']:(OUT/folder).mkdir(exist_ok=True)
    with source(TRACK/'selected_regions.tsv').open(encoding='utf-8') as f:
        regions={r['species']:r for r in csv.DictReader(f,delimiter='\t')}
    with source(CORR/'correlation_and_OLS.tsv').open(encoding='utf-8') as f:
        stats={(r['species'],r['strand']):r for r in csv.DictReader(f,delimiter='\t')}
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    put(fig,21,169.5,'Representative genomic region',bold=True,va='top')
    put(fig,113,169.5,'All held-out test windows',bold=True,va='top')
    all_panel_letters=[]
    for ri,(key,name,n_srx,n,chroms) in enumerate(SPECS):
        row=regions[key];track=load(TRACK/f'{key}.selected_track.npz')
        original=load(CORR/f'{key}.all_test_windows.npz')
        assert original['predicted'].shape==(n,2) and len(track['srx'])==n_srx
        np.testing.assert_allclose(track['cpm'],(track['counts']*1e6/track['denominator'][:,None,None]).mean(0),rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(track['log2p1_cpm'],np.log2(1+track['cpm']),atol=1e-12,rtol=0)
        assert np.all(track['end0']-track['start0']==128)
        assert np.isfinite(track['log2p1_cpm']).all()
        valid=np.asarray(track['valid'],dtype=bool)
        if valid.ndim==1:valid=np.repeat(valid[:,None],2,axis=1)
        assert valid.shape==track['pred'].shape
        assert np.isnan(track['pred'][~valid]).all()
        assert np.isfinite(track['pred'][valid]).all()
        top=160-ri*51
        letter='abcdefghi'[ri*3];all_panel_letters.append(letter)
        put(fig,0,top,letter,bold=True,size=8)
        put(fig,21,top,name,italic=True)
        coord=f"chr{row['chrom']}:{int(row['start0'])+1:,}-{int(row['end0']):,}"
        put(fig,21,top-3.5,coord)
        put(fig,147.75,top+3.5,f"chr{chroms}; n = {n:,}",ha='center')
        arrays=[track['log2p1_cpm'][:,0],track['pred'][:,0],track['log2p1_cpm'][:,1],track['pred'][:,1]]
        labels=['GRO-seq +','Prediction +','GRO-seq -','Prediction -']
        for ki,(values,label,color) in enumerate(zip(arrays,labels,COLORS)):
            lower=top-8-6.4-ki*8.5
            ax=axis(fig,21,lower,77,6.4)
            missing=np.zeros(len(values),dtype=bool)
            if ki in [1,3]:
                missing=~valid[:,0 if ki==1 else 1]
                for left,right in zip(track['start0'][missing],track['end0'][missing]):
                    ax.axvspan(left/1e6,right/1e6,facecolor='#E0E0E0',edgecolor='none',zorder=0)
            bars=ax.bar(track['start0']/1e6,values,width=(track['end0']-track['start0'])/1e6,
                        align='edge',color=color,edgecolor='none',linewidth=0)
            np.testing.assert_array_equal(np.asarray([b.get_height() for b in bars]),values)
            ax.set_xlim(int(row['start0'])/1e6,int(row['end0'])/1e6)
            ymax=1. if ki in [1,3] else max(.5,math.ceil(float(np.max(values))*1.08*2)/2)
            ax.set_ylim(0,ymax);ax.set_yticks([0,ymax],labels=['0',f'{ymax:g}'])
            clean_axis(ax);tick_y_inside(ax)
            put(fig,14.3,lower+3.2,label,ha='right',va='center')
            tick_step=10000 if int(row['end0'])-int(row['start0'])<=50000 else 20000
            ticks=np.arange(int(row['start0'])+tick_step//2,int(row['end0']),tick_step)/1e6
            ax.set_xticks(ticks)
            if ki<3:
                ax.tick_params(axis='x',labelbottom=False,length=0)
                ax.spines['bottom'].set_color('#888888')
            else:
                ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.3f}'))
                ax.set_xlabel(f"Position on chr{row['chrom']} (Mb)",fontsize=7,labelpad=1.7)
            artists.append({'species':key,'kind':'track','strand':'plus' if ki<2 else 'minus',
                'axis_label':label,'n':len(values),'missing_bins':int(missing.sum()),
                'bar_values_exactly_equal_source':True,'score_axis':[0,1] if ki in [1,3] else None})
        for ci,strand in enumerate(['plus','minus']):
            letter='abcdefghi'[ri*3+ci+1];all_panel_letters.append(letter)
            x=113+ci*39.0
            put(fig,x-9.5,top,letter,bold=True,size=8)
            put(fig,x+15.25,top,'Plus strand' if ci==0 else 'Minus strand',ha='center')
            raw=original['predicted'][:,ci].astype(float);yy=np.log2(1.+raw);xx=experimental_x(original,ci)
            r=float(pearsonr(xx,yy).statistic);rho=float(spearmanr(xx,yy).statistic)
            ref=stats[(key,strand)]
            assert abs(float(pearsonr(xx,raw).statistic)-float(ref['pearson_r']))<1e-12
            assert abs(rho-float(ref['spearman_rho']))<1e-12
            put(fig,x+15.25,top-3.5,f'r = {r:.3f}; ρ = {rho:.3f}',ha='center',color=METRIC_COLOR)
            ax=axis(fig,x,top-38,30.5,30)
            info=scatter_panel(ax,xx,raw,key,strand,'composite')
            if ci==0:ax.set_ylabel('log2(1 + predicted score)',fontsize=7,labelpad=1.5)
            clean_axis(ax);edges_inside(ax)
            checks.append({'species':key,'strand':strand,'n':n,'pearson_r':r,'spearman_rho':rho,
                'original_values_unchanged':True,**info})
    put(fig,147.75,12.8,'GRO-seq signal (log2[1 + mean CPM])',ha='center')
    put(fig,0,8.3,'strand_seiplant_v1; repeats retained; Arabidopsis / wheat / maize: 14 / 5 / 3 SRX')
    fig.add_artist(plt.Line2D([0,7/W],[3.7/H,3.7/H],transform=fig.transFigure,ls='--',color=FIT_COLOR,lw=.8))
    put(fig,9,3.7,'OLS fit',va='center')
    fig.add_artist(plt.Rectangle((55/W,2.8/H),4/W,1.8/H,transform=fig.transFigure,
        facecolor='#E0E0E0',edgecolor='none'))
    put(fig,61,3.7,'No valid DNA input',va='center')
    put(fig,113,4.2,'Within-panel density',va='center')
    count_legend(fig)
    # Assert physical font settings on every visible text artist before authoring exports.
    fig.canvas.draw()
    font_audit=[];bounds=[];renderer=fig.canvas.get_renderer()
    for t in fig.findobj(match=Text):
        if not t.get_visible() or not t.get_text():continue
        size=float(t.get_fontsize());assert size in [7.,8.],(t.get_text(),size)
        assert t.get_fontfamily()==['Arial'],(t.get_text(),t.get_fontfamily())
        if size==8:assert t.get_text() in all_panel_letters and t.get_fontweight()=='bold'
        bb=t.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        bounds.append({'text':t.get_text(),'bounds':list(bb.extents)})
        font_audit.append({'text':t.get_text(),'font_family':t.get_fontfamily(),'size_pt':size,'weight':t.get_fontweight()})
    mainpdf=OUT/'figures/GROseq_v1_tracks_OLS_correlations_183x170mm.pdf'
    fig.savefig(mainpdf,dpi=600)
    fig.savefig(mainpdf.with_suffix('.svg'),dpi=600)
    fig.savefig(mainpdf.with_suffix('.png'),dpi=600)
    plt.close(fig)
    # A4 proof: translate the original vector figure only. No scaling; 7pt remains 7pt.
    title='GRO-seq signal and sequence-model predictions across three plant species.'
    caption=(
        '(a, d, g) Representative regions from held-out chromosomes of Arabidopsis thaliana (40 kb), '
        'Triticum aestivum (80 kb) and Zea mays (80 kb), respectively. The wheat and maize regions '
        'were extended 40 kb to the right while retaining their original left boundary. '
        'Four tracks show GRO-seq and predicted scores '
        'on the plus and minus strands. Experimental signal counts aligned SE/R1 5-prime endpoints '
        'in 128-bp bins; each SRX is normalized to CPM before equal-weight averaging and log2(1 + mean CPM) transformation. '
        'The 14, 5 and 3 SRX used for Arabidopsis, wheat and maize match the old-model data sources. '
        'Predictions are unchanged 0-1 window scores from 1,024-bp DNA inputs advanced by 128 bp, '
        'assigned to their central 128 bp for display. Repeats are retained. '
        'Gray prediction-track bins indicate unavailable predictions where a 1,024-bp DNA input '
        'contains a non-ATCG reference base; experimental signal is retained and missing scores are not zero-filled. '
        '(b, c, e, f, h, i) Signal-prediction associations across all formal 1,024-bp test windows '
        '(45,066, 98,656 and 49,786 windows, respectively), with no subsampling. '
        'Correlation panels use log2(1 + original predicted score) on the y-axis and log2(1 + mean CPM) '
        'on the x-axis. The prediction tracks remain on their original 0-1 scale. r is Pearson correlation '
        'recalculated on the displayed transformed coordinates; rho is Spearman correlation. '
        'Dark-gray dashed straight lines are ordinary least-squares fits to all transformed test-window '
        'coordinates, drawn over the observed x range. Both x- and y-axes start at zero; '
        'point coordinates are unchanged. All scatter points remain visible. '
        'Gray shading encodes relative binned Gaussian density within '
        'each panel; no subsampling or outlier removal is used. '
        'Track regions were selected for experimental peaks '
        'and agreement and are illustrative examples. Experimental signal and predictions have different units; '
        'GRO-seq measurements also supplied the training labels.'
    )
    packet=io.BytesIO();can=canvas.Canvas(packet,pagesize=(210*MM,297*MM))
    can.setTitle('GRO-seq v1: tracks and correlations, A4 proof, Arial 7pt / panel 8pt')
    can.setAuthor('Research figure')
    x=13.5*MM;y=96*MM;line_width=183*MM
    can.setFont('Arial-Bold',7);can.drawString(x,y,title);y-=11
    can.setFont('Arial',7)
    words=caption.split();line=''
    lines=[]
    for word in words:
        candidate=(line+' '+word).strip()
        if pdfmetrics.stringWidth(candidate,'Arial',7)>line_width:
            lines.append(line);line=word
        else:line=candidate
    if line:lines.append(line)
    for line in lines:can.drawString(x,y,line);y-=10
    assert y>20*MM
    can.showPage();can.save();packet.seek(0)
    background=PdfReader(packet).pages[0]
    main_page=PdfReader(mainpdf).pages[0]
    assert abs(float(main_page.mediabox.width)/MM-W)<1e-5 and abs(float(main_page.mediabox.height)/MM-H)<1e-5
    background.merge_transformed_page(main_page,Transformation().translate(tx=13.5*MM,ty=106*MM),expand=False)
    writer=PdfWriter();writer.add_page(background)
    writer.add_metadata({'/Title':'GRO-seq v1: A4 composite','/Subject':'7pt Arial annotations; 8pt Arial Bold panel letters; original 183x170 mm figure, no scaling'})
    a4=OUT/'GROseq_v1_tracks_OLS_correlations_A4.pdf'
    with a4.open('wb') as f:writer.write(f)
    standalone_results=standalone(stats)
    document_checks={}
    for pdf in [mainpdf,a4,OUT/'figures/GROseq_v1_OLS_correlations_183x170mm.pdf',OUT/'GROseq_v1_OLS_correlations_A4.pdf']:
        with pdfplumber.open(pdf) as doc:
            assert len(doc.pages)==1
            page=doc.pages[0];fonts=sorted({c['fontname'] for c in page.chars})
            # PDFMiner reports rotated glyph advance as `size`; the physical font height is its width.
            actual_size=lambda c: c['size'] if c['upright'] else c['width']
            sizes=sorted({round(actual_size(c),3) for c in page.chars})
            assert sizes==[7.,8.],sizes
            assert all('Arial' in f for f in fonts),fonts
            letters=[c for c in page.chars if abs(actual_size(c)-8)<1e-5]
            expected_letters='abcdefghi' if 'tracks_' in pdf.name else 'abcdef'
            assert sorted(c['text'] for c in letters)==list(expected_letters)
            assert all('Bold' in c['fontname'] for c in letters)
            expected=(210,297) if pdf.name.endswith('_A4.pdf') else (183,170)
            np.testing.assert_allclose([page.width/MM,page.height/MM],expected,atol=1e-4,rtol=0)
            document_checks[pdf.name]={'pages':1,'size_mm':list(expected),'embedded_font_names':fonts,
                'actual_text_sizes_pt':sizes,'ordinary_annotation_size_pt':7,'panel_size_pt':8,
                'panel_letters':expected_letters,'text_editable':True,'character_count':len(page.chars)}
    for path,value in sources.items():assert sha(path)==value,'Existing source changed'
    dump(OUT/'validation/build_validation.json',{'state':'passed','source_files_unchanged':True,
        'fonts':font_audit,'text_bounds_normalized':bounds,'documents':document_checks,
        'scatter_statistics':checks,'standalone_statistics':standalone_results,
        'track_artists':artists,'source_files_sha256':sources,
        'reference_guidelines':['https://www.nature.com/nature/for-authors/final-submission',
            'https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/',
            'https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/']})
    (OUT/'figure_legend_EN.txt').write_text(title+'\n\n'+caption+'\n',encoding='utf-8')
    (OUT/'README_排版与数据说明.txt').write_text(
        '三物种最终图：两轴从零开始、log2(1+值)、OLS直线及向右扩展Track\n\n'
        'A4竖版210×297mm，1页。图主体183×170mm，不缩放放入A4；图注置于下方。\n'
        '全部普通标注实际7pt Arial；面板字母a–i实际8pt Arial Bold、正体小写。科学物种名为7pt Arial Italic。\n'
        'PDF文字保留可编辑并嵌入TrueType字体。SVG为可编辑文字；PNG仅用于预览。\n'
        '每行对应拟南芥、小麦、玉米；左侧4轨为实验+、预测+、实验−、预测−；右侧为全测试集的正、负链相关性。\n'
        '全部使用strand_seiplant_v1，重复序列保留；小麦仍为旧五文库，不是新三文库训练结果。\n'
        '拟南芥Track保持40kb。小麦chr3A:292140000–292220000、玉米chr6:74920000–75000000为80kb（以上0-based半开坐标）。\n'
        '两者保留原左边界，向右扩展40kb；所有扩展段均使用真实实验与原模型预测，不平移或补空白伪造信号。\n'
        '小麦新增输入若含参考N等非ATCG碱基，预测为NaN，图中标浅灰缺失区域；实验信号保留，不将缺失预测填0。\n'
        'Track实验信号统计128bp，模型输入1024bp/步长128bp；散点实验信号统计完整1024bp。\n'
        '横轴log2(1+mean CPM)，纵轴log2(1+原始预测分数)。两轴均直接从各自原始数值先加1再取log2，横轴没有重复取log。\n'
        '全部原始正式测试窗逐点展示，未抽样、未移除异常点。Pearson r按变换后的坐标重新计算，Spearman rho不变。\n'
        '综合图每物种一行；独立相关性图每物种一列、上正链下负链。Track的原始预测0-1刻度保留，相关性纵轴使用加1再取log2。\n'
        '灰度为全点256×256直方图经sigma=2格Gaussian平滑后在原坐标双线性取值；只控制颜色。深浅在每面板内归一化。\n'
        '深灰虚线为OLS直线：使用全部当前log2(1+值)坐标重新计算斜率与截距，未直接复用旧的原始分数拟合系数。\n'
        '六张相关性散点图的x、y轴显示下限均为0；原始点、相关系数与OLS拟合保持不变。\n'
        '直线画在有观测值的x范围内，没有增加y=x参考线或假定独立窗口的置信区间。\n'
        '按PDF的实际大小100%查看/打印，避免“适合页面”再次缩放。独立投稿优先使用figures/内183×170mm矢量PDF。\n'
        '依据Nature官方final-submission和research-figure-guide设置：普通字体5–7pt、面板字母8pt、双栏宽183mm、主图高不超过170mm。\n',encoding='utf-8')
    shutil.copyfile(Path(__file__),OUT/'provenance/build_final.py')
    print(json.dumps({'state':'created_and_structurally_validated','output':str(OUT),'pdfs':document_checks},ensure_ascii=False))

if __name__=='__main__':main()
