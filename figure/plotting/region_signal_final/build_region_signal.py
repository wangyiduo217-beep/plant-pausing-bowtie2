"""Six regional correlation/mean-signal panels on one A4 without font scaling."""
from pathlib import Path
import csv, hashlib, io, json, shutil, sys
import argparse, os
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path)
parser.add_argument('--font-dir',type=Path,default=os.environ.get('GROSEQ_FONT_DIR'))
args=parser.parse_args()
FIGURE=Path(__file__).resolve().parents[2]
BASE=FIGURE.parent
FONT_DIR=args.font_dir or Path(os.environ.get('WINDIR','C:/Windows'))/'Fonts'
for name in ('arial.ttf','arialbd.ttf','ariali.ttf'):
    if not (FONT_DIR/name).is_file():
        parser.error('Arial is required; use --font-dir or GROSEQ_FONT_DIR (font binaries are not bundled).')
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,font_manager as fm
from matplotlib.text import Text
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
import pdfplumber

DATA=FIGURE/'data/region_signal_final'
OUT=args.output_dir or FIGURE/'results/region_signal_final'
for d in ['data','figures','provenance','validation']:(OUT/d).mkdir(parents=True,exist_ok=True)
W,H,MM=183.,170.,72/25.4
PLOT_X=[12.,110.];PLOT_WIDTH=64.;PLOT_HEIGHT=28.5
REGIONS=['upstream_1_3kb','upstream_lt1kb','five_prime_utr','CDS','intron','three_prime_utr','downstream_le1kb','distal']
LABELS=['Up\n1-3kb','Up\n<1kb',"5'UTR\n",'CDS\n','Intron\n',"3'UTR\n",'Down\n<=1kb','Distal\n']
SPECS=[('arabidopsis','Arabidopsis thaliana','chr3',45066,14),
       ('wheat','Triticum aestivum','chr3A / 3B / 3D',98656,5),
       ('maize','Zea mays','chr6',49786,3)]
COLORS={'plus':'#21844B','minus':'#3278A5'};ORANGE='#D55E00'
for f in ['arial.ttf','arialbd.ttf','ariali.ttf']:fm.fontManager.addfont(str(FONT_DIR/f))
for name,f in [('Arial','arial.ttf'),('Arial-Bold','arialbd.ttf')]:pdfmetrics.registerFont(TTFont(name,str(FONT_DIR/f)))
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,
    'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.5,'xtick.major.width':.5,
    'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,
    'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def txt(fig,x,y,s,size=7,**kw):return fig.text(x/W,y/H,s,fontsize=size,fontfamily='Arial',**kw)

# Verify frozen input bytes and recalculate all 48 summaries before drawing.
from verify_region_inputs import verify
verification=verify(FIGURE,OUT/'validation/input_validation.json')
SOURCES={item['path']:item['sha256'] for item in verification['inputs']}
with (DATA/'regional_correlation_and_mean_signal.tsv').open(encoding='utf-8') as f:
    STATS=list(csv.DictReader(f,delimiter='\t'))
assert len(STATS)==48
shutil.copy2(DATA/'region_definition.json',OUT/'provenance/region_definition.json')

fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
txt(fig,0,169,'Regional correlation and mean GRO-seq signal',fontweight='bold',va='top')
fig.legend(handles=[Patch(facecolor='#E0E6E8',edgecolor='#657984',label='Pearson r'),
    Line2D([0],[0],color=ORANGE,marker='o',ms=2.5,lw=.8,label='Mean CPM')],
    loc='upper right',bbox_to_anchor=(1.,1.),ncol=2,frameon=False,
    handlelength=1.15,columnspacing=1.0,borderaxespad=0)
for ci,strand in enumerate(['plus','minus']):
    txt(fig,PLOT_X[ci]+PLOT_WIDTH/2,157,'Plus strand (+)' if ci==0 else 'Minus strand (-)',
        ha='center',fontweight='bold',color=COLORS[strand])
ARTISTS=[];CATEGORY_TEXT=[];PANELS=[]
for ri,(key,name,chrom,n,srx) in enumerate(SPECS):
    header=149-ri*50;bottom=header-34
    txt(fig,0,header,name,fontstyle='italic')
    txt(fig,W,header,f'{chrom}; n = {n:,}; {srx} SRX',ha='right')
    for ci,strand in enumerate(['plus','minus']):
        color=COLORS[strand];px=PLOT_X[ci]
        txt(fig,0 if ci==0 else 98,header-6,'abcdef'[ri*2+ci],size=8,fontweight='bold')
        ax=fig.add_axes([px/W,bottom/H,PLOT_WIDTH/W,PLOT_HEIGHT/H]);other=ax.twinx()
        sub=[next(s for s in STATS if s['species']==key and s['strand']==strand and s['region']==reg) for reg in REGIONS]
        xs=np.arange(8);rs=np.array([float(s['pearson_r']) for s in sub]);means=np.array([float(s['mean_cpm']) for s in sub])
        bars=ax.bar(xs,rs,width=.55,facecolor=to_rgba(color,.16),edgecolor=color,lw=.65,zorder=2)
        line,=other.plot(xs,means,color=ORANGE,marker='o',ms=2.5,lw=.8,zorder=3)
        np.testing.assert_array_equal([b.get_height() for b in bars],rs)
        np.testing.assert_array_equal(line.get_ydata(),means)
        ax.set_xlim(-.6,7.6);other.set_xlim(-.6,7.6)
        ax.set_ylim(0,1);other.set_ylim(0,6)
        ax.set_yticks([0,.25,.5,.75,1],labels=['0','0.25','0.50','0.75','1'])
        other.set_yticks([0,1.5,3,4.5,6],labels=['0','1.5','3','4.5','6'])
        ax.set_xticks([]);other.set_xticks([])
        ax.set_ylabel('Pearson r',color=color,labelpad=1.4)
        other.set_ylabel('Mean CPM',color=ORANGE,labelpad=1.4)
        ax.tick_params(axis='y',colors=color,length=1.8,pad=1.2)
        other.tick_params(axis='y',colors=ORANGE,length=1.8,pad=1.2)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        other.spines['top'].set_visible(False);other.spines['left'].set_visible(False)
        ax.spines['left'].set_color(color);other.spines['right'].set_color(ORANGE)
        ax.grid(False);other.grid(False)
        points=means/6
        for j,(r,m,s,label) in enumerate(zip(rs,means,sub,LABELS)):
            gx=px+(j+.6)/8.2*PLOT_WIDTH
            label_y=r+.025
            local=[points[j]]
            if j:local.append(points[j]+.35*(points[j-1]-points[j]))
            if j<7:local.append(points[j]+.35*(points[j+1]-points[j]))
            low,high=min(local)-.02,max(local)+.02
            if label_y < high+.035 and label_y+.10 > low-.01:
                label_y=max(label_y,high+.045)
            txt(fig,gx,bottom+label_y*PLOT_HEIGHT,f'{r:.2f}',ha='center',va='bottom',color=color,
                bbox={'facecolor':'white','edgecolor':'none','alpha':.95,'pad':.2})
            t=txt(fig,gx,bottom-1.2,label+'\n'+str(int(s['n_windows'])),ha='center',va='top',linespacing=1.12)
            CATEGORY_TEXT.append((ri,ci,j,t))
            ARTISTS.append(dict(species=key,strand=strand,region=s['region'],n_windows=int(s['n_windows']),
                pearson_r=float(r),mean_cpm=float(m),left_bar=float(r),right_line=float(m)))
        assert sum(int(s['n_windows']) for s in sub)==n
        assert ax.get_ylim()==(0.,1.) and other.get_ylim()==(0.,6.)
        PANELS.append(dict(species=key,strand=strand,bounds_mm=[px,bottom,PLOT_WIDTH,PLOT_HEIGHT]))
txt(fig,0,1.3,'Counts below categories denote n; Up/Down: upstream/downstream; mean CPM includes zeros.')
fig.canvas.draw();renderer=fig.canvas.get_renderer();FONT_AUDIT=[]
for t in fig.findobj(match=Text):
    if not t.get_visible() or not t.get_text():continue
    assert t.get_fontsize() in [7.,8.] and t.get_fontfamily()==['Arial']
    bb=t.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    assert bb.x0>=-.002 and bb.x1<=1.002 and bb.y0>=-.002 and bb.y1<=1.002,(t.get_text(),bb.extents)
    FONT_AUDIT.append(dict(text=t.get_text(),size_pt=t.get_fontsize(),bounds=list(bb.extents)))
# Neighbouring category labels, including exact counts, must not overlap.
for ri in range(3):
    for ci in range(2):
        labels=[t for r,c,j,t in CATEGORY_TEXT if r==ri and c==ci]
        for a,b in zip(labels,labels[1:]):
            assert a.get_window_extent(renderer).x1 < b.get_window_extent(renderer).x0, (a.get_text(),b.get_text())
main=OUT/'figures/region_signal_six_panels_183x170mm.pdf'
for ext in ['pdf','svg','png']:fig.savefig(main.with_suffix('.'+ext),dpi=600)
plt.close(fig)
TITLE='Regional correlation and experimental mean signal across plant species and strands.'
CAPTION=(
    '(a, b) Arabidopsis thaliana; (c, d) Triticum aestivum; (e, f) Zea mays. '
    'Left column: plus strand; right column: minus strand. Green/blue bars use the left axis (Pearson r); '
    'orange points and connecting lines use the right axis (mean experimental CPM). '
    'Pearson r is calculated across all 1,024-bp test windows in each region between log2(1+mean SRX CPM) '
    'and log2(1+predicted score). The orange values are arithmetic means of the original per-window mean SRX CPM, '
    'including zeros; the right-axis means are not log-transformed. Each SRX is normalized independently '
    'to CPM before equal-weight averaging (14, 5 and 3 SRX for Arabidopsis, wheat and maize). '
    'Both axes use common limits across all six panels: 0-1 for r and 0-6 CPM for mean signal. '
    'Up means upstream and Down means downstream, defined relative to the original canonical protein-coding '
    'transcript reference. Distances are in kb. Counts below category names are the full regional window counts n, '
    'shared by the two strands; no data are omitted. Total test windows are 45,066, 98,656 and 49,786, respectively. '
    'The same eight region assignments, values and models as the two-page version are retained; only the layout '
    'and category-label abbreviations have changed. The line connects category summaries as a visual guide. '
    'No score threshold, outlier removal, invented replicate/error bars or significance tests are used. '
    'These summaries describe covariation and do not establish a causal effect of signal intensity. '
    'strand_seiplant_v1; repeats retained. Ordinary text is 7-pt Arial; panel letters are 8-pt Arial Bold.'
)
buffer=io.BytesIO();c=canvas.Canvas(buffer,pagesize=(210*MM,297*MM));c.setTitle(TITLE)
px=13.5*MM;py=96*MM;c.setFont('Arial-Bold',7);c.drawString(px,py,TITLE);py-=11;c.setFont('Arial',7)
line='';lines=[]
for word in CAPTION.split():
    candidate=(line+' '+word).strip()
    if pdfmetrics.stringWidth(candidate,'Arial',7)>183*MM:lines.append(line);line=word
    else:line=candidate
if line:lines.append(line)
for line in lines:c.drawString(px,py,line);py-=10
assert py>20*MM
c.showPage();c.save();buffer.seek(0)
page=PdfReader(buffer).pages[0]
page.merge_transformed_page(PdfReader(main).pages[0],Transformation().translate(tx=13.5*MM,ty=106*MM))
writer=PdfWriter();writer.add_page(page)
with (OUT/'region_signal_six_panels_A4.pdf').open('wb') as f:writer.write(f)
(OUT/'figure_legend_EN.txt').write_text(TITLE+'\n\n'+CAPTION+'\n',encoding='utf-8')
DOCS={}
for p,mm in [(main,[183,170]),(OUT/'region_signal_six_panels_A4.pdf',[210,297])]:
    with pdfplumber.open(p) as doc:
        assert len(doc.pages)==1;pg=doc.pages[0]
        chars=[ch for ch in pg.chars if ch['text'].strip()]
        sizes=sorted({round(ch['size'] if ch['upright'] else ch['x1']-ch['x0'],3) for ch in chars})
        assert sizes==[7.,8.] and all('Arial' in ch['fontname'] for ch in chars)
        assert ''.join(ch['text'] for ch in chars if round(ch['size'] if ch['upright'] else ch['x1']-ch['x0'],3)==8)=='abcdef'
        np.testing.assert_allclose([pg.width/MM,pg.height/MM],mm,atol=1e-4,rtol=0)
        assert all(ch['x0']>=-.1 and ch['x1']<=pg.width+.1 and ch['top']>=-.1 and ch['bottom']<=pg.height+.1 for ch in chars)
        DOCS[str(p.relative_to(OUT))]=dict(n_pages=1,size_mm=mm,font_sizes_pt=sizes,all_Arial=True,sha256=sha(p))
for p,h in SOURCES.items():assert sha(BASE/p)==h
shutil.copy2(DATA/'region_definition.json',OUT/'provenance/region_definition.json')
shutil.copy2(__file__,OUT/'provenance/build_one_page.py')
dump(OUT/'validation/build_validation.json',dict(state='passed',layout='three species rows by two strand columns',
    ordinary_font_pt=7,panel_font_pt=8,all_Arial=True,no_graphics_scaling=True,all_category_labels_nonoverlapping=True,
    data_source_sha256=SOURCES,all_48_values_unchanged=True,artists=ARTISTS,panels=PANELS,fonts=FONT_AUDIT,documents=DOCS))
print(json.dumps(dict(state='passed',A4_pages=1,panels=6,stats=48,output=str(OUT)),ensure_ascii=False))

# Separate panel PDFs preserve the numeric drawing and original text as PDF objects.
# The parent figure is freshly drawn above; these are not exports of a screenshot.
import fitz
component_dir=OUT/'components';component_dir.mkdir(exist_ok=True)
with fitz.open(main) as full:
    panel_map=[]
    for ri,(key,name,chrom,n,srx) in enumerate(SPECS):
        header=149-ri*50;bottom=header-34
        for ci,strand in enumerate(('plus','minus')):
            letter='abcdef'[ri*2+ci]
            clip=fitz.Rect((0 if ci==0 else 98)*MM,(H-header-3)*MM,
                           (85 if ci==0 else 183)*MM,(H-bottom+10)*MM)
            one=fitz.open();pg=one.new_page(width=clip.width,height=clip.height)
            pg.show_pdf_page(pg.rect,full,0,clip=clip)
            file=component_dir/f'{letter}_{key}_{strand}.pdf';one.save(file);one.close()
            panel_map.append(dict(panel=letter,species=key,strand=strand,file=file.name))
dump(OUT/'validation/panel_map.json',panel_map)
