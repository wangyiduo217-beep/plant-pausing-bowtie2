"""Reassemble the selected final GRO-seq A4 figure without changing data.

Default: exact preserved vector/PDF assets for workflow, model, tracks and
scatter; saved numerical vertices for ROC/PR. --source-figure2 permits a
fresh numerical rebuild from build_source_panels.py. No training or inference.
"""
from pathlib import Path
import argparse, os, sys, io, json, hashlib, re, shutil, subprocess
import pymupdf as fitz
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, default=REPO/'figure/results/figure1')
parser.add_argument('--font-dir', type=Path, default=Path(os.environ.get('GROSEQ_FONT_DIR', 'C:/Windows/Fonts')),
                    help='Directory with user-installed arial.ttf, arialbd.ttf and ariali.ttf; fonts are not redistributed.')
parser.add_argument('--source-figure2', type=Path, default=HERE/'assets/source_figure2_tracks_correlations.pdf',
                    help='Optional newly rebuilt, zero-origin Figure2 A4 PDF with the original physical layout.')
parser.add_argument('--poppler-preview', action='store_true', help='Also render a preview with pdftoppm from PATH.')
parser.add_argument('--no-panels', action='store_true', help='Skip individual a-n panel exports.')
args = parser.parse_args()
OUT = args.output_dir.resolve()
FONT_DIR = args.font_dir.resolve()
SOURCES = [HERE/'assets/source_figure1_workflow_model.pdf', args.source_figure2.resolve()]
DATA = REPO/'figure/data/figure1/curves'
for p in (OUT, OUT/'validation', OUT/'source_data'):
    p.mkdir(parents=True, exist_ok=True)
for name, file in [('Arial','arial.ttf'),('Arial-Bold','arialbd.ttf'),('Arial-Italic','ariali.ttf')]:
    font_path = FONT_DIR/file
    if not font_path.is_file():
        parser.error(f'Missing {font_path}. Supply --font-dir with licensed Arial fonts.')
    pdfmetrics.registerFont(TTFont(name, str(font_path)))
for path in SOURCES:
    if not path.is_file():
        parser.error(f'Missing source asset: {path}')
trapezoid = getattr(np, 'trapezoid', None)
if trapezoid is None:
    trapezoid = np.trapz
MM = 72/25.4
PW, PH = 210*MM, 297*MM
src1,src2 = [fitz.open(p) for p in SOURCES]
source_spans = [[s for b in d[0].get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']] for d in (src1,src2)]

# Remove old panel letters and replace orange correlation annotations with black editable text.
p=src2[0]
reinsert=[]
for s in source_spans[1]:
    if (s['size']>7.5 and s['text'] in 'abcdefghi') or s['text'].startswith('r = '):
        p.add_redact_annot(fitz.Rect(s['bbox']),fill=False)
        if s['text'].startswith('r = '): reinsert.append(s)
p.apply_redactions(images=0,graphics=0)
for s in reinsert:
    p.insert_text(s['origin'],s['text'],fontsize=7,fontname='mergeArial',fontfile=str(FONT_DIR/'arial.ttf'),color=(0,0,0))

doc=fitz.open(); page=doc.new_page(width=PW,height=PH)
placements=[]
def place(source, clip, x, y, scale, label):
    cr=fitz.Rect(clip)
    dest=fitz.Rect(x*MM,y*MM,x*MM+cr.width*scale,y*MM+cr.height*scale)
    page.show_pdf_page(dest,source,0,clip=cr)
    placements.append(dict(label=label,source=1 if source is src1 else 2,clip=list(cr),destination=list(dest),scale=scale))
    return dest

# Overview and model preserve original vector drawing and all substantive labels.
overview=place(src1,(38,51,554,154), (210-516*.75/MM)/2,10.5,.75,'a: overview')
model_x=(210-541*.84/MM)/2
upper=place(src1,(27,174,568,266),model_x,42.2,.84,'b: model, upper row')
lower=place(src1,(27,285,568,385),model_x,73,.84,'b: model, lower row')

# The three rows keep tracks and scatter point layers directly from the user's PDF.
row_tops=[112,152,192]
track_x,scatter_x,roc_x=15,90,151
for i,y in enumerate(row_tops):
    dy=i*144.567
    place(src2,(38,80+dy,317,220+dy),track_x,y+10*.75/MM,.75,f'row {i+1}: tracks')
    place(src2,(331,70+dy,556,206+dy),scatter_x,y,.75,f'row {i+1}: plus and minus scatter')

buf=io.BytesIO(); c=canvas.Canvas(buf,pagesize=(PW,PH),initialFontName='Arial')
c.setTitle('Plant GRO-seq model and evaluation: combined A4 figure')
text_records=[]
def text(x,y,t,size=6,font='Arial',align='left',color=(0,0,0)):
    c.setFillColorRGB(*color);c.setFont(font,size)
    xx,yy=x*MM,PH-y*MM
    getattr(c,{'left':'drawString','center':'drawCentredString','right':'drawRightString'}[align])(xx,yy,t)
    text_records.append(dict(text=t,size=size,x_mm=x,y_mm=y,font=font))
def line(x1,y1,x2,y2,width=.5,color=(0,0,0),dash=None):
    c.setStrokeColorRGB(*color);c.setLineWidth(width);c.setDash(dash or [])
    c.line(x1*MM,PH-y1*MM,x2*MM,PH-y2*MM)
    c.setDash([])
def panel(x,y,s):text(x,y,s,size=8,font='Arial-Bold')

panel(15,13,'a');panel(15,44.6,'b')
# Clear only the old schematic return-line segment so it can be routed with padding.
# This does not touch graph data or labels.
old_left=model_x+(34.8199997-27)*.84/MM
return_y=73+(337.4573364-285)*.84/MM
c.setFillColorRGB(1,1,1)
c.rect((old_left-.25)*MM,PH-(return_y+.1)*MM,.5*MM,(return_y-73+.15)*MM,stroke=0,fill=1)
# Six consistent module outlines and a black enclosure follow the user's sketch.
# White interiors keep the original vector elements legible; colour indicates grouping.
def box(x,y,w,h,color,radius=2.5,lw=.65):
    c.setStrokeColorRGB(*color);c.setLineWidth(lw);c.setDash([])
    c.roundRect(x*MM,PH-(y+h)*MM,w*MM,h*MM,radius*MM,stroke=1,fill=0)
box(24.0,41.1,18.0,28.7,(.56,.59,.62))
box(50.0,41.1,32.0,28.7,(.29,.63,.78))
box(84.0,41.1,82.8,28.7,(.80,.60,.16))
box(27.5,72.2,43.5,31.4,(.15,.64,.39))
box(77.4,72.2,21.4,31.4,(.55,.38,.69))
box(104.0,72.2,85.0,31.4,(.78,.34,.33))
c.setStrokeColorRGB(0,0,0);c.setLineWidth(.75)
c.rect(20*MM,PH-105*MM,174.3*MM,(105-8.7)*MM,stroke=1,fill=0)
line(20,39,194.3,39,width=.65)
# Join the original model's retained right and left connector segments.
right=model_x+(566.453003-27)*.84/MM
left=22.8
line(right,upper.y1/MM,right,71.1,width=.5)
line(right,71.1,left,71.1,width=.5)
line(left,71.1,left,return_y,width=.5)
line(left,return_y,old_left+.32,return_y,width=.5)

text(51.5,110,'Representative genomic region',6,'Arial-Bold','center')
text(120,110,'All held-out test windows',6,'Arial-Bold','center')
text(175,110,'ROC / PR (micro)',6,'Arial-Bold','center')
metrics={'arabidopsis':(.8361396751,.7334185919),'wheat':(.9070427095,.7950871811),'maize':(.8776964811,.7482158014)}
roc_checks=[]
for i,(species,top) in enumerate(zip(metrics,row_tops)):
    for xx,letter in zip([15,90,119.25,151],'cdefghijklmn'[i*4:(i+1)*4]):panel(xx,top+4.8,letter)
    with np.load(DATA/f'{species}.HiLearn_curves.npz') as z:
        vals={k:z[k].copy() for k in z.files}
    target = OUT/'source_data'/f'{species}.HiLearn_curves.npz'
    if target.resolve() != (DATA/f'{species}.HiLearn_curves.npz').resolve():
        shutil.copy2(DATA/f'{species}.HiLearn_curves.npz',target)
    # Same HiLearn plotted vertices, including artificial PR endpoints; no new evaluation.
    auroc=-float(trapezoid(vals['roc_tpr'],vals['roc_fpr']))
    prauc=-float(trapezoid(vals['pr_plot_precision'],vals['pr_plot_recall']))
    assert abs(auroc-metrics[species][0])<1e-9 and abs(prauc-metrics[species][1])<1e-9
    roc_checks.append(dict(species=species,auroc=auroc,pr_auc=prauc,vertices_unchanged=True))
    x0,y0,w,h=163,top+10.75,24,24
    line(x0,y0,x0,y0+h);line(x0,y0+h,x0+w,y0+h)
    for v in [0,.5,1]:
        xx=x0+v*w;yy=y0+(1-v)*h
        line(xx,y0+h,xx,y0+h+.7,width=.4)
        line(x0-.7,yy,x0,yy,width=.4)
        text(xx,y0+h+2.9,f'{v:g}',5.25,align='center')
        text(x0-1.2,yy+.6,f'{v:g}',5.25,align='right')
    line(x0,y0+h,x0+w,y0,width=.5,color=(.5,.5,.5),dash=[2,1.5])
    for x,y,color,dash in [(vals['roc_fpr'],vals['roc_tpr'],(.9,.25,.13),[]),(vals['pr_plot_recall'],vals['pr_plot_precision'],(.20,.32,.57),[2.5,1.5])]:
        c.setStrokeColorRGB(*color);c.setLineWidth(.75);c.setDash(dash)
        path=c.beginPath();path.moveTo((x0+x[0]*w)*MM,PH-(y0+(1-y[0])*h)*MM)
        for a,b in zip(x[1:],y[1:]):path.lineTo((x0+a*w)*MM,PH-(y0+(1-b)*h)*MM)
        c.drawPath(path);c.setDash([])
    # All four spines remain visible, including the top/right edges requested by the user.
    c.setStrokeColorRGB(0,0,0);c.setLineWidth(.5);c.setDash([])
    c.rect(x0*MM,PH-(y0+h)*MM,w*MM,h*MM,stroke=1,fill=0)
    text(175,top+7.5,'AUROC %.3f; PR AUC %.3f'%(auroc,prauc),5.25,align='center')
    text(175,y0+h+6,'FPR (ROC) / Recall (PR)',5.25,align='center')
    c.saveState();c.translate((x0-5)*MM,PH-(y0+h/2)*MM);c.rotate(90);c.setFillColorRGB(0,0,0);c.setFont('Arial',5.25)
    c.drawCentredString(0,0,'TPR (ROC) / Precision (PR)');c.restoreState()

# Shared keys avoid repeated legends in the small ROC/PR plots.
text(120,234.5,'GRO-seq signal (log2[1 + mean CPM])',5.5,align='center')
line(15,238,21,238,dash=[2,1.5]);text(22.5,238.7,'OLS fit',5.5)
c.setFillColorRGB(.87,.87,.87);c.rect(42*MM,PH-239*MM,4*MM,2*MM,fill=1,stroke=0)
text(47.5,238.7,'No valid DNA input',5.5)
text(92,238.7,'Within-panel density',5.5)
for j in range(80):
    v=float(np.interp(j/79,[0,1/3,2/3,1],[196,146,75,17]))/255
    c.setFillColorRGB(v,v,v);c.rect((114+j*.25)*MM,PH-239*MM,.26*MM,1.5*MM,stroke=0,fill=1)
text(114,241.6,'Low',5.25);text(134,241.6,'High',5.25,align='right')
line(151,238,156,238,color=(.9,.25,.13),width=.75);text(157,238.7,'ROC',5.5)
line(166,238,171,238,color=(.20,.32,.57),width=.75,dash=[2.5,1.5]);text(172,238.7,'PR',5.5)
line(180,238,185,238,color=(.5,.5,.5),dash=[2,1.5]);text(186,238.7,'Random',5.25)
text(15,244,'strand_seiplant_v1; repeats retained; Arabidopsis / wheat / maize: 14 / 5 / 3 SRX',5.75)

caption=(
'a, Plant GRO-seq datasets and strand-resolved nascent transcription. b, Sequence model: 1,024-bp one-hot DNA inputs produce two 0-1 window scores. '
'c, g, k, Representative held-out regions of Arabidopsis thaliana, Triticum aestivum and Zea mays (40, 80 and 80 kb). Experimental tracks show log2(1 + mean CPM) in 128-bp bins, counting aligned SE/read-1 5-prime endpoints; CPM is calculated independently per SRX and averaged equally. Prediction tracks retain raw window scores, with a 128-bp sliding step and display in the central 128 bp. Gray bins mark unavailable DNA inputs and are not zero-filled. '
'd, e, h, i, l, m, Plus- and minus-strand correlations over all 45,066, 98,656 and 49,786 formal test windows, respectively. Both mean CPM and original predicted scores are transformed by log2(1 + value) before Pearson r; rho denotes Spearman correlation. Dashed lines are OLS fits. All source scatter points, zero-origin axes and relative within-panel density shading are retained. '
'f, j, n, HiLearn 0.2.2 ROC and PR curves pooling both strands (micro). Positives overlap same-strand consensus intervals supported by at least two SRX by at least 1 bp. Curves use raw scores; PR AUC is the native trapezoidal area including artificial endpoints, not average precision. '
'Track examples were selected for experimental peaks and agreement. Evaluation uses held-out genomic sequences from the experiments that supplied training labels; it is not independent-library validation.'
)
title='Plant GRO-seq sequence model and evaluation across three species.'
text(15,250,title,7,'Arial-Bold')
def wrap(t,width,size):
    lines=[];cur=''
    for word in t.split():
        nxt=(cur+' '+word).strip()
        if pdfmetrics.stringWidth(nxt,'Arial',size)>width:lines.append(cur);cur=word
        else:cur=nxt
    if cur:lines.append(cur)
    return lines
caption_lines=wrap(caption,180*MM,6)
for i,t in enumerate(caption_lines):text(15,254+i*2.75,t,6)
assert 254+(len(caption_lines)-1)*2.75<289
c.showPage();c.save();overlay=fitz.open(stream=buf.getvalue(),filetype='pdf')
page.show_pdf_page(page.rect,overlay,0)
doc.set_metadata({'title':title,'author':'Plant GRO-seq research','subject':'Combined original source figures; embedded Arial; editable vector graphics; 210 x 297 mm A4'})
doc.subset_fonts()
a4=OUT/'GROseq_combined_v2_A4.pdf'
doc.save(a4,garbage=4,deflate=True)

cropped=fitz.open();q=cropped.new_page(width=180*MM,height=238*MM)
q.show_pdf_page(q.rect,doc,0,clip=fitz.Rect(15*MM,8*MM,195*MM,246*MM))
cropped.set_metadata({'title':title,'subject':'180 x 238 mm artwork only; no caption; exact same scale as A4'})
cropped.save(OUT/'GROseq_combined_v2_artwork_180x238mm.pdf',garbage=4,deflate=True)
(OUT/'figure_legend_EN.txt').write_text(title+'\n\n'+caption+'\n',encoding='utf8')
(OUT/'original_Figure2_legend_EN.txt').write_text('\n'.join(s['text'] for s in source_spans[1] if s['bbox'][1]>560),encoding='utf8')
(OUT/'validation/placements.json').write_text(json.dumps(placements,ensure_ascii=False,indent=2),encoding='utf8')
final=fitz.open(a4);spans=[s for b in final[0].get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
visible=[s for s in spans if s['text'].strip()]
sizes=sorted(set(round(s['size'],3) for s in visible))
report={'page_count':len(final),'page_mm':[final[0].rect.width/MM,final[0].rect.height/MM],'visible_font_sizes_pt':sizes,'font_names':sorted(set(s['font'] for s in visible)),'min_font_pt':min(s['size'] for s in visible),'max_font_pt':max(s['size'] for s in visible),'roc_checks':roc_checks,'sources_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES},'caption_words':len(caption.split()),'caption_lines':len(caption_lines),'source_scatter_and_tracks_reused_as_pdf_objects':True,'data_recalculated':False}
assert report['min_font_pt']>=5-1e-3
assert report['max_font_pt']<=8+1e-3
assert len(final)==1
(OUT/'validation/validation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
final[0].get_pixmap(dpi=160).save(OUT/'GROseq_combined_A4_preview.png')
final[0].get_pixmap(dpi=600).save(OUT/'GROseq_combined_A4_600dpi.png')
if args.poppler_preview:
    subprocess.run(['pdftoppm','-r','160','-png','-singlefile',str(a4),str(OUT/'validation/poppler_preview')],check=True)


# Standalone exports have separate, non-overlapping scientific panel regions.
# Track groups contain four tracks; each scatter export contains ONE output.
# Shared species/n and x-axis labels are supplied to standalone scatters.
panel_exports = []
if not args.no_panels:
    panel_dir = OUT/'panels'
    panel_dir.mkdir(exist_ok=True)
    def export_clip(label, source_doc, rect, kind, species=None, n=None, scale=1.0, shared_x=False):
        rect = fitz.Rect(rect)
        header_height = 24 if species else 0
        footer_height = 15 if shared_x else 0
        padding = 4 if species else 0
        width = rect.width*scale + 2*padding
        if shared_x:
            width = max(width, 99.5)
        left_pad = width - rect.width*scale - padding
        height = rect.height*scale + header_height + footer_height
        result = fitz.open()
        dest = result.new_page(width=width, height=height)
        dest.show_pdf_page(fitz.Rect(left_pad,header_height,width-padding,height-footer_height),source_doc,0,clip=rect)
        if species:
            annot = io.BytesIO()
            ac = canvas.Canvas(annot,pagesize=(width,height),initialFontName='Arial')
            ac.setFillColorRGB(0,0,0)
            ac.setFont('Arial-Bold',8)
            ac.drawString(3,height-9,label)
            ac.setFont('Arial-Italic',5.5)
            ac.drawCentredString(width/2,height-8,species)
            if n is not None:
                ac.setFont('Arial',5.25)
                ac.drawCentredString(width/2,height-17,f'n = {n:,} test windows')
            if shared_x:
                ac.setFont('Arial',5.25)
                ac.drawCentredString(width/2,4,'GRO-seq signal (log2[1 + mean CPM])')
                if kind == 'minus_scatter':
                    ac.saveState()
                    ac.translate(8,footer_height+(rect.height*scale)/2)
                    ac.rotate(90)
                    ac.drawCentredString(0,0,'log2(1 + predicted score)')
                    ac.restoreState()
            ac.showPage(); ac.save()
            layer = fitz.open(stream=annot.getvalue(),filetype='pdf')
            dest.show_pdf_page(dest.rect,layer,0)
        result.set_metadata({'title':f'Panel {label}: {kind}', 'subject':'Exact source-graph objects; standalone labels added where shared in composite'})
        result.subset_fonts()
        target=panel_dir/f'panel_{label}_{kind}.pdf'
        result.save(target,garbage=4,deflate=True)
        panel_exports.append({'panel':label,'kind':kind,'file':target.relative_to(OUT).as_posix(),'source_clip_points':list(rect),'content_scale':scale,'width_mm':width/MM,'height_mm':height/MM})
    export_clip('a',final,(14*MM,8.4*MM,195*MM,39.3*MM),'workflow')
    export_clip('b',final,(14*MM,39.5*MM,195*MM,105.3*MM),'architecture')
    for i,(species,n,letters) in enumerate([('A. thaliana',45066,'cdef'),('T. aestivum',98656,'ghij'),('Z. mays',49786,'klmn')]):
        dy=i*144.567
        # Original panel letters and orange r labels have already been removed
        # from src2 above; the same r annotations were reinserted in black.
        export_clip(letters[0],src2,(38,80+dy,317,220+dy),'tracks',species=species,scale=.75)
        export_clip(letters[1],src2,(331,80+dy,453,206+dy),'plus_scatter',species=species,n=n,scale=.75,shared_x=True)
        export_clip(letters[2],src2,(453,80+dy,556,206+dy),'minus_scatter',species=species,n=n,scale=.75,shared_x=True)
        top=112+i*40
        export_clip(letters[3],final,(150.5*MM,(top+1.5)*MM,194.5*MM,(top+41.5)*MM),'roc_pr')
    (OUT/'validation/panel_exports.json').write_text(json.dumps(panel_exports,indent=2)+'\n',encoding='utf-8')

report['portable_rebuild'] = True
report['single_panel_exports'] = len(panel_exports)
report['source_figure2_override'] = SOURCES[1] != (HERE/'assets/source_figure2_tracks_correlations.pdf').resolve()
report['font_directory'] = str(FONT_DIR)
report['reproduction_scope'] = 'Original vector assets for a/b; source Figure2 objects for tracks/scatter; numerical saved vertices for ROC/PR.'
(OUT/'validation/validation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
