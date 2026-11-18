"""Display the existing strand-specific decile estimates side by side; no pooling."""
from pathlib import Path
import argparse, csv, hashlib, io, json, shutil, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.text import Text
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
import pdfplumber

W,H=183.,115.
MM=72/25.4
COLORS={'plus':'#21844B','minus':'#3278A5'}
OFFSETS={'plus':-.18,'minus':.18}
SPECS=[('arabidopsis','Arabidopsis thaliana'),('wheat','Triticum aestivum'),('maize','Zea mays')]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def put(fig,x,y,s,size=7,bold=False,italic=False,ha='left',**kw):
    return fig.text(x/(fig.get_figwidth()*25.4),y/(fig.get_figheight()*25.4),s,fontsize=size,fontfamily='Arial',fontweight='bold' if bold else 'normal',
                    fontstyle='italic' if italic else 'normal',ha=ha,**kw)

def a4_copy(mainfile,out,caption):
    packet=io.BytesIO(); c=canvas.Canvas(packet,pagesize=(210*MM,297*MM))
    c.setTitle('CPM decile correlations: grouped plus and minus strands')
    x=13.5*MM;y=153*MM
    c.setFont('Arial-Bold',7);c.drawString(x,y,'Strand-specific correlations displayed together within each decile.');y-=11
    c.setFont('Arial',7);line=''
    for word in caption.split():
        trial=(line+' '+word).strip()
        if pdfmetrics.stringWidth(trial,'Arial',7)>183*MM:
            c.drawString(x,y,line);y-=10;line=word
        else:line=trial
    if line:c.drawString(x,y,line)
    c.showPage();c.save();packet.seek(0)
    page=PdfReader(packet).pages[0]
    page.merge_transformed_page(PdfReader(mainfile).pages[0],Transformation().translate(tx=13.5*MM,ty=166*MM),expand=False)
    writer=PdfWriter();writer.add_page(page)
    writer.add_metadata({'/Title':'Grouped strand-specific decile correlations','/Subject':'Unchanged estimates; separate strand-specific rank deciles; 7 pt Arial; 8 pt panels'})
    file=out/'GROseq_CPM_deciles_grouped_strands_A4.pdf'
    with file.open('wb') as f:writer.write(f)
    return file

def draw_species_panel(fig,col,species,name,lookup,artists,count_artists,left,standalone=False):
    fw,fh=fig.get_figwidth()*25.4,fig.get_figheight()*25.4
    mid=left+23.75
    put(fig,left-6,99,'abc'[col],size=8,bold=True)
    put(fig,mid,99,name,italic=True,ha='center')
    ax=fig.add_axes([left/fw,27/fh,47.5/fw,42/fh])
    ax.set_xlim(.45,10.55);ax.set_ylim(-.25,.50)
    ax.set_xticks(range(1,11),labels=[f'Q{i}' for i in range(1,11)])
    ax.set_yticks([-.2,0,.2,.4],labels=['-0.2','0','0.2','0.4'])
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=.5,length=2,pad=1.6,labelsize=7)
    if col==0 or standalone:ax.set_ylabel('Pearson r',labelpad=2)
    ax.axhline(0,color='#B0B0B0',lw=.6,ls='--',zorder=0)
    for si,strand in enumerate(['plus','minus']):
        color=COLORS[strand];offset=OFFSETS[strand];n_y=1.22 if strand=='plus' else 1.03
        for gi in range(1,11):
            r=lookup[(species,strand,gi)]
            x=gi+offset;y=float(r['pearson_r']);lo=float(r['ci95_low']);hi=float(r['ci95_high']);n=int(r['n_windows'])
            ax.vlines(x,lo,hi,color=color,lw=.9,zorder=2)
            ax.hlines([lo,hi],x-.07,x+.07,color=color,lw=.9,zorder=2)
            ax.plot(x,y,'o',color=color,ms=3.2,mec='white',mew=.35,zorder=3)
            t=ax.text(gi,n_y,f'{n:,}',transform=ax.get_xaxis_transform(),ha='center',va='bottom',rotation=90,
                      fontsize=7,color=color)
            count_artists.append(t)
            artists.append(dict(species=species,strand=strand,group_id=gi,x=x,x_offset=offset,
                pearson_r=y,ci95_low=lo,ci95_high=hi,n_windows=n))
        ax.text(-.015,n_y,'n (+)' if strand=='plus' else 'n (-)',transform=ax.transAxes,
                ha='right',va='bottom',fontsize=7,color=color)
        zero=lookup[(species,strand,0)]
        assert zero['pearson_r']=='' and zero['spearman_rho']==''
        put(fig,mid,17 if strand=='plus' else 12,
            f"CPM = 0 ({'+' if strand=='plus' else '-'}): n = {int(zero['n_windows']):,}; r = NA",ha='center',color=color)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    module_dir=Path(__file__).resolve().parent
    figure_dir=module_dir.parents[1]
    ap.add_argument('--data-dir',type=Path,default=figure_dir/'data/figureS2')
    ap.add_argument('--output-dir',type=Path,default=figure_dir/'results/figureS2')
    ap.add_argument('--font-dir',type=Path,help='Directory containing Arial regular, bold and italic TTF files.')
    ap.add_argument('--panel',choices=['all','a','b','c'],default='all',help='all generates the complete figure plus all three standalone panels.')
    args=ap.parse_args();source=args.data_dir.resolve();out=args.output_dir.resolve()
    for sub in ['data','code','validation']:(out/sub).mkdir(parents=True,exist_ok=True)
    hashes={}
    for name in ['group_statistics.csv','zero_signal_prediction_summary.csv','overall_correlation_reference.csv']:
        src=source/name;dest=out/'data'/name;hashes[name]=sha(src)
        if src.resolve()!=dest.resolve():shutil.copy2(src,dest)
        assert sha(dest)==hashes[name]
    with (out/'data/group_statistics.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    assert len(rows)==66
    lookup={(r['species'],r['strand'],int(r['group_id'])):r for r in rows}
    assert len(lookup)==66
    # Arial must be available locally; do not silently substitute a different font.
    font_paths={}
    for key,file,weight,style in [('Arial','arial.ttf','normal','normal'),('Arial-Bold','arialbd.ttf','bold','normal'),('Arial-Italic','ariali.ttf','normal','italic')]:
        if args.font_dir:
            path=args.font_dir/file
            if not path.is_file(): raise FileNotFoundError(path)
            fm.fontManager.addfont(str(path))
        else:
            path=Path(fm.findfont(fm.FontProperties(family='Arial',weight=weight,style=style),fallback_to_default=False))
        font_paths[key]=str(path)
        pdfmetrics.registerFont(TTFont(key,str(path)))
    plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.linewidth':.5,
        'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
    panel_values=[]
    for col,(species,name) in enumerate(SPECS):
        letter='abc'[col]
        if args.panel!='all' and args.panel!=letter: continue
        panel=plt.figure(figsize=(70/25.4,105/25.4),facecolor='white')
        values=[];counts=[]
        draw_species_panel(panel,col,species,name,lookup,values,counts,left=15,standalone=True)
        handles=[Line2D([0],[0],marker='o',linestyle='none',color=COLORS[s],markersize=3.5,
                        label='Plus strand (+)' if s=='plus' else 'Minus strand (-)') for s in ['plus','minus']]
        panel.legend(handles=handles,loc='center',bbox_to_anchor=(.55,91/105),ncol=2,frameon=False,
                     fontsize=7,handlelength=1.,handletextpad=.4,columnspacing=1.)
        put(panel,37.5,5,'Nonzero CPM decile (Q1 lowest; Q10 highest)',ha='center')
        panel.canvas.draw()
        stem=out/f'panel_{letter}_{species}'
        for ext in ['pdf','svg','png']: panel.savefig(str(stem)+'.'+ext,dpi=300)
        plt.close(panel)
        panel_values.extend(values)
    dump(out/'validation/panel_artist_values.json',panel_values)
    if args.panel!='all':
        dump(out/'validation/panel_run.json',dict(panel=args.panel,n_points=len(panel_values),input_sha256=hashes,
             statistics_recomputed=False,strands_pooled=False))
        return
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    put(fig,15,109,'Prediction-experiment association across nonzero CPM deciles',bold=True)
    put(fig,15,104,'strand_seiplant_v1 | 1,024-bp test windows | strand-specific estimates')
    handles=[Line2D([0],[0],marker='o',linestyle='none',color=COLORS[s],markersize=3.5,
                    label='Plus strand (+)' if s=='plus' else 'Minus strand (-)') for s in ['plus','minus']]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.54,91/H),ncol=2,frameon=False,
               fontsize=7,handlelength=1.,handletextpad=.4,columnspacing=2.5)
    artists=[];count_artists=[]
    for col,(species,name) in enumerate(SPECS):
        draw_species_panel(fig,col,species,name,lookup,artists,count_artists,left=15+57.5*col)
    put(fig,W/2,5,'Nonzero CPM decile (Q1 lowest; Q10 highest)',ha='center')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();font_records=[]
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text():continue
        box=t.get_window_extent(renderer)
        if box.width<=0 or box.height<=0:continue
        assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.width+1 and box.y1<=fig.bbox.height+1,(t.get_text(),box)
        font_records.append(dict(text=t.get_text(),size_pt=t.get_fontsize()))
    # Within each panel the two rows of exact n labels must not collide.
    for i,a in enumerate(count_artists):
        ba=a.get_window_extent(renderer)
        for b in count_artists[i+1:]:
            bb=b.get_window_extent(renderer)
            assert not ba.overlaps(bb),(a.get_text(),b.get_text())
    stem=out/'GROseq_CPM_deciles_grouped_strands_183x115mm'
    fig.savefig(str(stem)+'.pdf');fig.savefig(str(stem)+'.svg');fig.savefig(str(stem)+'.png',dpi=600)
    plt.close(fig)
    caption=('a-c, Arabidopsis, wheat and maize. Within each Q1-Q10 group, green (left) and blue (right) points show the previously calculated plus- and minus-strand Pearson correlations, respectively. '
        'Whiskers are the original 95% confidence intervals from 1,000 chromosome-stratified, paired 100-kb block-bootstrap replicates. '
        'Two color-coded rows above each panel give the exact strand-specific window counts. '
        'The strands were neither pooled nor re-ranked jointly: positive-CPM windows retain their original ten equal-count rank groups within each species and strand. '
        'Thus matching Q labels indicate relative rank within each strand, not necessarily identical absolute CPM intervals. '
        'Both experiment (mean SRX CPM) and original predicted score were transformed by log2(1 + value) before calculating correlations, matching the overall-test analysis. '
        'Zero-CPM windows remain separate, with undefined correlations (NA). No statistics or input values were changed. '
        'strand_seiplant_v1; repeats retained; Arabidopsis/wheat/maize use 14/5/3 SRX, including the original five-library wheat model. Arial 7 pt; panel letters 8 pt.')
    (out/'figure_legend_EN.txt').write_text(caption+'\n',encoding='utf-8')
    a4=a4_copy(Path(str(stem)+'.pdf'),out,caption)
    font_checks=[]
    for file in [Path(str(stem)+'.pdf'),a4]:
        with pdfplumber.open(file) as doc:
            p=doc.pages[0];chars=p.chars
            sizes=sorted(set(round(c['size'] if c['upright'] else c['width'],3) for c in chars))
            assert sizes==[7.,8.],sizes
            assert all('Arial' in c['fontname'] for c in chars)
            assert all(c['x0']>=-.2 and c['x1']<=p.width+.2 and c['top']>=-.2 and c['bottom']<=p.height+.2 for c in chars)
            font_checks.append(dict(file=file.name,page_pt=[p.width,p.height],sizes_pt=sizes,fontnames=sorted(set(c['fontname'] for c in chars))))
    dump(out/'validation/PDF_font_and_page_checks.json',font_checks)
    dump(out/'validation/artist_values.json',artists)
    dump(out/'validation/font_sizes.json',font_records)
    dump(out/'validation/source_identity_and_plot_checks.json',dict(state='passed',source_dir=str(source),source_sha256=hashes,
        copied_tables_byte_identical=True,pooled_strands=False,recomputed_statistics=False,grouping_still_strand_specific=True,
        x_offsets=OFFSETS,n_plot_points=len(artists),exact_counts_in_two_color_rows=True))
    (out/'README_说明.txt').write_text(
        '本图仅将同一物种的正负链相关性用group/dodge形式放在同一面板。绿色左点=正链，蓝色右点=负链。\n'
        '正负链数据不合并、不共同重新排序、不重算相关性或bootstrap。r、置信区间及n与原十分位结果逐值一致。\n'
        'Q1-Q10仍是各条链内部的非零CPM十分位，同名Q不保证两链的绝对CPM区间完全相同。\n'
        '两行竖排n分别为绿色正链、蓝色负链。零CPM组保留数量，r为NA。\n'
        '实验mean CPM及预测原分数都先log2(1+值)再计算Pearson，与整体相关性方法一致。\n'
        'data/group_statistics.csv是原统计表字节级拷贝，包括组内CPM上下界、r、CI、n。\n'
        '原模型strand_seiplant_v1保留重复序列，拟南芥/小麦/玉米为14/5/3SRX，小麦为原5文库。\n'
        '普通7pt Arial，面板字母8pt；根目录A4 PDF页面精确210x297mm，图183x115mm放置但不缩放。\n'
        'run_deciles_final.py读取已有统计表重画完整图和单面板，不重复bootstrap。\n',encoding='utf-8')
    target=out/'code/run_deciles_final.py'
    if Path(__file__).resolve()!=target.resolve():shutil.copy2(Path(__file__),target)
    assert all(sha(source/name)==h for name,h in hashes.items())
    print(json.dumps(dict(state='passed',output=str(out),points=len(artists),data_pooled=False,statistics_recomputed=False)))

if __name__=='__main__':main()
