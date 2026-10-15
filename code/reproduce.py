"""Publication figures from stored aggregate inputs; no fitting or simulation."""
from pathlib import Path
import argparse,csv,hashlib,json,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
import warnings
import numpy as np
from PIL import Image,ImageDraw,ImageFont

def configure_font():
 # Select an available serif font for portable rendering.
 import os
 names=['Times New Roman','Liberation Serif','STIX Two Text','STIXGeneral','DejaVu Serif']
 if os.environ.get('FIGURE_FORCE_SERIF_FALLBACK')=='1':names=names[1:]
 for name in names:
  try:font_manager.findfont(name,fallback_to_default=False)
  except ValueError:continue
  if name!='Times New Roman':warnings.warn('Times New Roman is unavailable; using '+name+'. Layout is measured in the selected font. Scientific source values are unchanged.',RuntimeWarning)
  return name
 raise RuntimeError('No approved serif font is available')
FONT=configure_font()
plt.rcParams.update({'font.family':FONT,'font.size':10,'axes.labelsize':10,'axes.titlesize':10.5,'xtick.labelsize':9.5,'ytick.labelsize':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'black','axes.labelcolor':'black','text.color':'black','xtick.color':'black','ytick.color':'black','axes.linewidth':0.7,'svg.fonttype':'none','svg.hashsalt':'exact-marginal-variance-figures','ps.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','mathtext.fontset':'custom','mathtext.rm':FONT,'mathtext.it':FONT+':italic','mathtext.bf':FONT+':bold','mathtext.sf':FONT,'mathtext.tt':FONT,'mathtext.fallback':'stix'})

def measured_wrap(fig,text,width_inches,size,weight='normal'):
 # Font metrics, rather than character counts, determine line wrapping.
 from matplotlib.font_manager import FontProperties
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 prop=FontProperties(family=FONT,size=size,weight=weight)
 limit=width_inches*fig.dpi
 result=[]
 for paragraph in text.split('\n'):
  line=''
  for word in paragraph.split():
   candidate=(line+' '+word).strip()
   if line and renderer.get_text_width_height_descent(candidate,prop,False)[0]>limit:result.append(line);line=word
   else:line=candidate
  result.append(line)
 return '\n'.join(result)
METHODS=['oracle','separate_R1_upper','separate_R2_upper','unrestricted_upper']
LABEL={'oracle':'Retained correspondence','separate_R1_upper':'R1: exact records','separate_R2_upper':'R2: complete tables','unrestricted_upper':'Unrestricted rearrangement'}
COLOR={'oracle':'#111111','separate_R1_upper':'#0072B2','separate_R2_upper':'#D55E00','unrestricted_upper':'#666666'}
MARK={'oracle':'o','separate_R1_upper':'s','separate_R2_upper':'^','unrestricted_upper':'D'}
LINE={'oracle':'-','separate_R1_upper':'--','separate_R2_upper':':','unrestricted_upper':'-.'}
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def f(r,k):return float(r[k])
def grid(ax):ax.grid(axis='y',color='#dddddd',linewidth=.5);ax.set_axisbelow(True)
def panel(ax,letter,title):ax.set_title(letter+'. '+title,loc='left',pad=10,fontweight='bold')
def setup(nrows,ncols,height):
 fig,axs=plt.subplots(nrows,ncols,figsize=(7.1,height),squeeze=False)
 return fig,axs
def main(data,output):
 data=Path(data);output=Path(output);output.mkdir(parents=True,exist_ok=True)
 definition=json.loads((data/'scenario_definitions.json').read_text(encoding='utf-8'))
 cells={c['id']:c for c in definition['cells']}
 summaries=json.loads((data/'method_summaries.json').read_text())
 lookup={(r['cell'],r['method']):r for r in summaries}
 trace=[];inventory=[];layout=[]
 def bind(fig,panel_id,source,records,fields):
  for row in records:
   key=';'.join(f'{k}={row[k]}' for k in ['cell','method','anchor','n','metric','mode'] if k in row)
   for field in fields:
    if field in row and row[field]!='':trace.append({'figure':fig,'panel':panel_id,'canonical_file':source,'row_key':key,'field':field,'exact_source_value':str(row[field]),'display_rule':'Stored value; linear coordinate mapping only'})
 def save(fig,name,panels,source_files):
  where='main' if not name.startswith('S') else 'supplementary'
  dest=output/where;dest.mkdir(exist_ok=True)
  base=dest/('Figure_'+name)
  fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.get_size_inches()*fig.dpi
  nondrawn=set()
  for axes in fig.axes:
   for axis in [axes.xaxis,axes.yaxis]:
    lo,hi=sorted(axis.get_view_interval())
    for tick in axis.get_major_ticks()+axis.get_minor_ticks():
     if not lo<=tick.get_loc()<=hi:nondrawn.update([tick.label1,tick.label2])
  outside=[]
  for text in fig.findobj(matplotlib.text.Text):
   if not text.get_visible() or not text.get_text() or text in nondrawn:continue
   b=text.get_window_extent(renderer)
   if b.x0<-.5 or b.y0<-.5 or b.x1>w+.5 or b.y1>h+.5:outside.append(text.get_text())
  # Tight export measures the actual font extents and expands the margins if
  # another approved serif needs more room than the nominal canvas.
  export_layout={'bbox_inches':'tight','pad_inches':.10}
  fig.savefig(base.with_suffix('.svg'),**export_layout,metadata={'Date':None,'Creator':None})
  t=base.with_suffix('.svg').read_text(encoding='utf-8')
  t=re.sub(r'<metadata>.*?</metadata>','',t,flags=re.S)
  base.with_suffix('.svg').write_text(t,encoding='utf-8',newline='\n')
  fig.savefig(base.with_suffix('.eps'),**export_layout,metadata={'Creator':''})
  eps_path=base.with_suffix('.eps')
  eps_lines=eps_path.read_text(encoding='utf-8').splitlines(keepends=True)
  eps_path.write_text(''.join('%%Title: Comes; Figure '+name+'; top is the upper page edge\n' if line.startswith('%%Title:') else line for line in eps_lines),encoding='utf-8')
  fig.savefig(base.with_suffix('.png'),dpi=200,**export_layout,metadata={'Software':None})
  fig.savefig(base.with_suffix('.tiff'),dpi=600,**export_layout,pil_kwargs={'compression':'tiff_lzw'})
  with Image.open(base.with_suffix('.tiff')) as original:
   rgba=original.convert('RGBA');rgb=rgba.convert('RGB');assert rgb.tobytes()==rgba.convert('RGB').tobytes()
   from PIL import TiffImagePlugin
   tags=TiffImagePlugin.ImageFileDirectory_v2()
   tags[270]='Comes; Figure '+name+'; top is the upper image edge'
   tags[274]=1
   tags[315]='Calin-Adrian Comes'
   rgb.save(base.with_suffix('.tiff'),compression='tiff_lzw',dpi=(600,600),tiffinfo=tags)
  with Image.open(base.with_suffix('.tiff')) as im:
   pixels=im.size;dpi=im.info.get('dpi');assert dpi and min(dpi)>=600
  fonts=sorted({text.get_fontfamily()[0] for text in fig.findobj(matplotlib.text.Text) if text.get_text()})
  layout.append({'axes':[{'xlim':list(a.get_xlim()),'ylim':list(a.get_ylim())} for a in fig.axes], 'legend_primitives':[{'markers':[h.get_marker() for h in leg.legend_handles], 'linestyles':[h.get_linestyle() for h in leg.legend_handles]} for leg in fig.legends], 'TIFF_mode':'RGB', 'figure':name,'panels':panels,'fonts':fonts,'font_size_min':min(text.get_fontsize() for text in fig.findobj(matplotlib.text.Text) if text.get_visible() and text.get_text()),'text_inside_canvas':not outside,'fallback_tight_export':FONT!='Times New Roman','all_text_retained_in_export':True,'inches':fig.get_size_inches().tolist(),'TIFF_pixels':pixels,'DPI':[{'numerator':v.numerator,'denominator':v.denominator} for v in dpi]})
  for suffix in ['eps','svg','tiff','png']:
   file=base.with_suffix('.'+suffix)
   inventory.append({'figure':'Figure '+name,'panels':panels,'canonical_data_source':';'.join(source_files),'dimensions_inches':f'{fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f}' if FONT=='Times New Roman' else f'{pixels[0]/600:.2f} x {pixels[1]/600:.2f}','raster_pixels':f'{pixels[0]} x {pixels[1]}' if suffix=='tiff' else '', 'raster_DPI':600 if suffix=='tiff' else 200 if suffix=='png' else '', 'fonts':';'.join(fonts),'final_filename':file.relative_to(output).as_posix(),'sha256':sha(file)})
  plt.close(fig);print('Rendered Figure',name,flush=True)
 def legend(fig,handles,labels,y=.015,ncol=3):
  clean=[]
  for handle in handles:
   h=handle.lines[0] if hasattr(handle,'lines') else handle
   marker=h.get_marker()
   if marker not in [None,'None','none','',' ']:
    clean.append(Line2D([],[],color=h.get_color(),marker=marker,linestyle='none',markersize=5.2))
   else:
    clean.append(Line2D([],[],color=h.get_color(),linestyle=h.get_linestyle(),linewidth=h.get_linewidth()))
  fig.legend(clean,labels,loc='lower center',bbox_to_anchor=(.5,y),ncol=ncol,frameon=False,fontsize=9.5,columnspacing=1.3,handlelength=1.0,handletextpad=.65)
 def coverage(ax,source,ids,xlabels,figure,panel_id,connect=False):
  handles=[];labels=[];rr=rows(data/source)
  available=[m for m in METHODS if any(r['method']==m for r in rr)]
  offsets=np.linspace(-.13,.13,len(available))
  for m,offset in zip(available,offsets):
   rec=[next(r for r in rr if int(r['cell'])==c and r['method']==m and r.get('coverage')) for c in ids]
   y=np.array([f(r,'coverage') for r in rec]);low=np.array([f(r,'Wilson_lower') for r in rec]);high=np.array([f(r,'Wilson_upper') for r in rec])
   x=np.arange(len(ids))+offset
   h=ax.plot(x,y,color=COLOR[m],marker=MARK[m],linestyle=LINE[m] if connect else 'none',linewidth=1.1,markersize=4.4,label=LABEL[m])[0]
   ax.vlines(x,low,high,color=COLOR[m],linewidth=1)
   ax.hlines(low,x-.025,x+.025,color=COLOR[m],linewidth=1);ax.hlines(high,x-.025,x+.025,color=COLOR[m],linewidth=1)
   handles.append(h);labels.append(LABEL[m]);bind(figure,panel_id,source,rec,['coverage','Wilson_lower','Wilson_upper','covered','method_completed'])
  ax.axhline(.95,color='#777777',ls='--',lw=.8,zorder=0)
  ax.set_xticks(range(len(ids)),xlabels);ax.set_ylabel('Coverage among completed intervals');ax.set_ylim(.88,1.012);grid(ax)
  return handles,labels
 def core_ids(rate):return sorted(c for c in range(100,112) if cells[c]['censor_rate']==rate)
 def core_labels(ids):return [r'$\lambda_W$='+f'{cells[c]["lambda_W"]:g}' for c in ids]
 def core_ticks(ax,ids):
  ax.set_xticks(range(len(ids)),core_labels(ids))
  for j in range(0,len(ids),2):
   rule_transform=ax.get_xaxis_transform()+ScaledTranslation(0,-23/72,ax.figure.dpi_scale_trans)
   label_transform=ax.get_xaxis_transform()+ScaledTranslation(0,-35/72,ax.figure.dpi_scale_trans)
   ax.plot([j-.38,j+1.38],[0,0],transform=rule_transform,color='#999999',lw=.55,clip_on=False)
   ax.text(j+.5,0,f'q={cells[ids[j]]["death_first_probability"]:g}',transform=label_transform,ha='center',va='center',fontsize=9.5,fontweight='bold' if cells[ids[j]]['death_first_probability']>0 else 'normal')


 # Main Figure 1: unchanged mean-per-realization gains and stored error bars.
 fig,axs=setup(1,2,3.5);fig.subplots_adjust(left=.085,right=.985,bottom=.26,top=.88,wspace=.32)
 gain=rows(data/'core_gain.csv');handles=[]
 for idx,rate in enumerate([.25,.75]):
  ax=axs[0,idx];ids=core_ids(rate)
  for m,off in [('separate_R1_upper',-.06),('separate_R2_upper',.06)]:
   rr=[next(r for r in gain if int(r['cell'])==c and r['method']==m) for c in ids]
   h=ax.errorbar(np.arange(6)+off,[f(r,'gain_percent') for r in rr],yerr=[f(r,'MC_halfwidth') for r in rr],color=COLOR[m],marker=MARK[m],linestyle=LINE[m],ms=4.5,capsize=2,lw=1.1,label=LABEL[m]);handles.append(h)
   bind('1',chr(65+idx),'core_gain.csv',rr,['gain_percent','MC_halfwidth','censor_rate','predeclared_main_family'])
  core_ticks(ax,ids);ax.set_ylim(0,30);ax.set_yticks([0,10,20,30]);ax.set_ylabel('Reduction in interval width (%)' if idx==0 else '');ax.axhline(0,color='black',lw=.6);grid(ax)
  panel(ax,chr(65+idx),r'Censoring hazard $\lambda_C$ = '+str(rate))
 legend(fig,handles[:2],[LABEL[m] for m in ['separate_R1_upper','separate_R2_upper']],ncol=2)
 save(fig,'1','A, B',['core_gain.csv','scenario_definitions.json'])

 fig,axs=setup(2,2,5.75);fig.subplots_adjust(left=.105,right=.985,bottom=.16,top=.92,wspace=.35,hspace=.78)
 coverage_data=rows(data/'core_coverage_cost.csv')
 for col,rate in enumerate([.25,.75]):
  ids=core_ids(rate);labels=core_labels(ids)
  h,l=coverage(axs[0,col],'core_coverage_cost.csv',ids,labels,'2',chr(65+col));core_ticks(axs[0,col],ids);panel(axs[0,col],chr(65+col),r'Coverage, $\lambda_C$ = '+str(rate))
  ax=axs[1,col]
  ratios=[m for m in METHODS if any(r['method']==m and r.get('mean_width_ratio') for r in coverage_data)]
  for m,off in zip(ratios,np.linspace(-.1,.1,len(ratios))):
   rr=[next(r for r in coverage_data if int(r['cell'])==c and r['method']==m and r.get('mean_width_ratio')) for c in ids]
   ax.errorbar(np.arange(6)+off,[f(r,'mean_width_ratio') for r in rr],yerr=[f(r,'ratio_MC_halfwidth') for r in rr],color=COLOR[m],marker=MARK[m],ls=LINE[m],ms=4.4,capsize=2,lw=1.1)
   bind('2',chr(67+col),'core_coverage_cost.csv',rr,['mean_width_ratio','ratio_MC_halfwidth'])
  core_ticks(ax,ids);ax.set_ylabel('Width / retained-correspondence width' if col==0 else '');ax.axhline(1,color='#777777',ls='--',lw=.8);ax.set_ylim(.95,2.20);ax.set_yticks([1,1.25,1.5,1.75,2]);grid(ax);panel(ax,chr(67+col),r'Width ratio, $\lambda_C$ = '+str(rate))
 legend(fig,h,l,y=.005,ncol=2 if len(l)==4 else 3);save(fig,'2','A, B, C, D',['core_coverage_cost.csv','scenario_definitions.json'])

 fig,axs=setup(1,2,3.85);fig.subplots_adjust(left=.105,right=.985,bottom=.27,top=.86,wspace=.32)
 for col,ids in enumerate([[112,104,113],[114,111,115]]):
  ax=axs[0,col];h,l=coverage(ax,'sample_size_coverage.csv',ids,['100','200','500'],'3',chr(65+col),True)
  panel(ax,chr(65+col),r'$q=0.25,\ \lambda_W=0.5,\ \lambda_C=0.25$' if col==0 else r'$q=0.5,\ \lambda_W=2,\ \lambda_C=0.75$')
  ax.set_xlabel('Sample size, n');ax.set_ylabel('Coverage among completed intervals' if col==0 else '')

  if col==0:ax.annotate('933/1,000',xy=(2-.13,.933),xytext=(.85,.903),arrowprops={'arrowstyle':'-','color':'black','lw':.7},fontsize=10)
 legend(fig,h,l,y=.01);save(fig,'3','A, B',['sample_size_coverage.csv','scenario_definitions.json'])

 fig,axs=setup(2,1,5.4);fig.subplots_adjust(left=.29,right=.98,bottom=.14,top=.91,hspace=.65)
 rr=rows(data/'rotterdam_intervals.csv');ids=['oracle','separate_R1_upper','separate_R2_upper','unrestricted_upper','independence'];labels=['Retained correspondence','R1: exact records','R2: complete tables','Unrestricted rearrangement','Zero-covariance benchmark']
 ax=axs[0,0]
 for j,m in enumerate(ids):
  r=next(r for r in rr if r['method']==m);color=COLOR.get(m,'#555555');ax.plot([f(r,'lower'),f(r,'upper')],[j,j],color=color,lw=2);ax.plot(f(r,'point'),j,MARK.get(m,'D'),color=color,ms=4)
  bind('4','A','rotterdam_intervals.csv',[r],['point','lower','upper','width'])
 ax.set_yticks(range(5),labels);ax.invert_yaxis();ax.set_xlabel('Restricted mean time alive after recurrence (days)');ax.axvline(f(rr[0],'point'),color='#aaaaaa',ls=':',lw=.8);grid(ax);panel(ax,'A','Same-center intervals')
 ax=axs[1,0];ranges=rows(data/'rotterdam_variance_endpoints.csv')
 for j,r in enumerate(ranges):
  color=['#0072B2','#D55E00','#666666'][j];ax.plot([f(r,'variance_lower'),f(r,'variance_upper')],[j,j],color=color,lw=2);ax.plot([f(r,'variance_lower'),f(r,'variance_upper')],[j,j],['s','^','D'][j],color=color,ms=4)
  bind('4','B','rotterdam_variance_endpoints.csv',[r],['variance_lower','variance_upper'])
 ax.set_yticks(range(3),['R1: exact records','R2: complete tables','Unrestricted rearrangement']);ax.invert_yaxis();ax.set_xlabel('Empirical variance (days²)');ax.set_ylim(2.5,-.5);panel(ax,'B','Attained variance endpoints')
 a=ax.axvline(f(rr[0],'variance'),color='black',ls='--',lw=.8,label='Retained correspondence');b=ax.axvline(f(next(r for r in rr if r['method']=='independence'),'variance'),color='#555555',ls=':',lw=.9,label='Zero-covariance benchmark')
 bind('4','B','rotterdam_intervals.csv',[rr[0],next(r for r in rr if r['method']=='independence')],['variance'])
 legend(fig,[a,b],['Retained correspondence','Zero-covariance benchmark'],y=.005,ncol=2);save(fig,'4','A, B',['rotterdam_intervals.csv','rotterdam_variance_endpoints.csv'])

 fig=plt.figure(figsize=(7.1,2.6));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off()
 boxes=[(.015,'Linked endpoint records','Marginal estimates and\ncase-weight derivatives\n\nPatient correspondence retained\n\nPaired empirical variance','#555555'),(.355,'R1: exact unlinked records','Exact marginal records\nCase-weight derivatives retained\n\nPatient correspondence removed\n\nCompatible variance endpoints','#0072B2'),(.695,'R2: complete event-risk tables','Marginal derivatives retained\nComplete event history\nand horizon checkpoints\n\nPatient correspondence absent\n\nCompatible variance endpoints','#D55E00')]
 for x,title,body,color in boxes:
  ax.add_patch(Rectangle((x,.26),.29,.51,facecolor='white',edgecolor='#b9bec3',lw=.7))
  ax.plot([x,x+.29],[.77,.77],color=color,lw=1.5)
  title=measured_wrap(fig,title,.27*7.1,10,'bold')
  body=measured_wrap(fig,body.replace('\n\n','\n'),.27*7.1,9.5)
  ax.text(x+.145,.700,title,ha='center',va='center',fontsize=10,fontweight='bold',linespacing=1.15)
  ax.text(x+.145,.480,body,ha='center',va='center',fontsize=9.5,linespacing=1.10)
 for x0,x1,label in [(.307,.351,'Remove patient\ncorrespondence'),(.647,.691,'Aggregate censoring information\nwithin event-risk cells')]:
  center=(x0+x1)/2
  ax.annotate('',xy=(x1,.52),xytext=(x0,.52),arrowprops={'arrowstyle':'->','lw':1,'color':'#333333'})
  label=measured_wrap(fig,label,.30*7.1,9.5)
  ax.text(center,.865,label,ha='center',va='center',fontsize=9.5,linespacing=1.2)
  ax.plot([center,center],[.81,.565],color='#777777',lw=.6)
 ax.text(.5,.140,measured_wrap(fig,'Sparse published Kaplan–Meier output does not satisfy either exact-release specification.',6.8,9.5),ha='center',va='center',fontsize=9.5)
 ax.set_xlim(0,1);ax.set_ylim(0,1)
 save(fig,'S1','Unlettered information-flow schematic',['main_information_table.csv'])

 fig,axs=setup(2,2,5.35);fig.subplots_adjust(left=.1,right=.985,bottom=.14,top=.91,wspace=.32,hspace=.58)
 scale=rows(data/'sample_size_scaling.csv')
 for col,anchor in enumerate(['Low','High']):
  rr=[r for r in scale if r['anchor']==anchor]
  for row,(metric,err,ylabel) in enumerate([('sqrt_n_RMSE','sqrt_n_RMSE_MC_halfwidth','√n × RMSE'),('sqrt_n_oracle_width','sqrt_n_oracle_width_MC_halfwidth','√n × mean paired interval width')]):
   ax=axs[row,col];ax.errorbar(range(3),[f(r,metric) for r in rr],yerr=[f(r,err) for r in rr],marker='o',color='black',linestyle='-',lw=1.1,ms=4.4,capsize=2);ax.set_xticks(range(3),['100','200','500']);ax.set_xlabel('Sample size, n');ax.set_ylabel(ylabel if col==0 else '');grid(ax)
   letter=chr(65+row*2+col);panel(ax,letter,r'$q=0.25,\ \lambda_W=0.5$' if col==0 else r'$q=0.5,\ \lambda_W=2$');bind('S2',letter,'sample_size_scaling.csv',rr,[metric,err,'n','completed'])
 for row,(metric,err) in enumerate([('sqrt_n_RMSE','sqrt_n_RMSE_MC_halfwidth'),('sqrt_n_oracle_width','sqrt_n_oracle_width_MC_halfwidth')]):
  lower=min(f(r,metric)-f(r,err) for r in scale);upper=max(f(r,metric)+f(r,err) for r in scale);pad=.06*(upper-lower)
  for col in range(2):axs[row,col].set_ylim(lower-pad,upper+pad)
 save(fig,'S2','A, B, C, D',['sample_size_scaling.csv'])

 fig,axs=setup(1,2,3.2);fig.subplots_adjust(left=.19,right=.985,bottom=.24,top=.86,wspace=.30)
 rr=rows(data/'precision_distributions.csv')
 for col,metric in enumerate(['width','width_oracle_ratio']):
  ax=axs[0,col]
  for j,m in enumerate(METHODS):
   r=next(r for r in rr if r['method']==m and r['metric']==metric);ax.plot([f(r,'min'),f(r,'max')],[j,j],color=COLOR[m],lw=.75);ax.plot([f(r,'p05'),f(r,'p95')],[j,j],color=COLOR[m],lw=3.5);ax.plot(f(r,'median'),j,MARK[m],color=COLOR[m],ms=4.8);bind('S3',chr(65+col),'precision_distributions.csv',[r],['min','p05','median','p95','max','n'])
  wrapped_labels={'oracle':'Retained\ncorrespondence','unrestricted_upper':'Unrestricted\nrearrangement'}
  ax.set_yticks(range(4),[wrapped_labels.get(m,LABEL[m]) for m in METHODS] if col==0 else ['']*4);ax.invert_yaxis();ax.set_xlabel('Interval width' if col==0 else 'Width / paired interval width');panel(ax,chr(65+col),'Interval width' if col==0 else 'Width ratio')
 save(fig,'S3','A, B',['precision_distributions.csv'])

 for name,ids,source in [('S4',list(range(116,123)),'shape_support_coverage.csv'),('S5',[123,124],'informative_censoring_coverage.csv')]:
  if name=='S4':
   fig,axs=setup(2,1,6.15);fig.subplots_adjust(left=.12,right=.985,bottom=.14,top=.93,hspace=.75);axA,axB=axs[:,0]
   labels=['Common\ncensoring\n'+r'$\lambda_C=0.25$','Common\ncensoring\n'+r'$\lambda_C=0.75$','Discrete\nh=0.1','Discrete\nh=0.25','Weibull\nshape=0.75','Weibull\nshape=1.5','Limited\nsupport\nn=100, τ=3']
  else:
   fig,axs=setup(1,2,3.55);fig.subplots_adjust(left=.105,right=.985,bottom=.26,top=.87,wspace=.38);axA,axB=axs[0]
   labels=[r'Base $\lambda_C=0.25$',r'Base $\lambda_C=0.75$']
  h,l=coverage(axA,source,ids,labels,name,'A');panel(axA,'A','Conditional coverage')
  for m,off in zip(METHODS[:3],[-.13,0,.13]):
   vals=[lookup[c,m]['covered_planned_operational_yield'] for c in ids];axB.plot(np.arange(len(ids))+off,vals,MARK[m],color=COLOR[m],ms=4.4,ls='none')
   bind(name,'B','method_summaries.json',[lookup[c,m] for c in ids],['covered_planned_operational_yield','covered','planned'])
  axB.set_xticks(range(len(ids)),labels);axB.set_ylabel('Covered / all planned samples');axB.set_ylim(0,1.012);grid(axB);panel(axB,'B','Covered / planned');legend(fig,h,l,y=.005);save(fig,name,'A, B',[source,'method_summaries.json'])

 fig,axs=setup(2,1,6.0);fig.subplots_adjust(left=.11,right=.985,bottom=.18,top=.93,hspace=.78)
 rr=rows(data/'bootstrap_operations.csv');ids=[1,102,104,111,119,122]
 labels=['Uncensored\nFGM, a=1','q=0\n'+r'$\lambda_W=2$'+'\n'+r'$\lambda_C=0.25$','q=0.25\n'+r'$\lambda_W=0.5$'+'\n'+r'$\lambda_C=0.25$','q=0.5\n'+r'$\lambda_W=2$'+'\n'+r'$\lambda_C=0.75$','Discrete\nh=0.25','Limited support\nn=100, τ=3']
 handles=[];legend_labels=[]
 modes=[('literal','Modified order-filtered','o','#111111'),('support_conditioned','Support-conditioned','s','#0072B2'),('paired','Paired bootstrap','^','#D55E00')]
 for row,metric in enumerate(['cap_rate','unsupported_fraction_of_accepted']):
  ax=axs[row,0]
  for (mode,label,marker,color),off in zip(modes,[-.14,0,.14]):
   selected=[next(r for r in rr if int(r['cell'])==c and r['mode']==mode) for c in ids];h=ax.plot(np.arange(6)+off,[f(r,metric) for r in selected],marker,ls='none',ms=4.7,color=color,label=label)[0]
   if row==0:handles.append(h);legend_labels.append(label)
   bind('S6',chr(65+row),'bootstrap_operations.csv',selected,[metric])
  ax.set_xticks(range(6),labels);ax.set_ylim(-.025,.75);ax.set_ylabel('Fraction of executed streams capped' if row==0 else 'Unsupported / accepted resamples');grid(ax);panel(ax,chr(65+row),'Attempt cap' if row==0 else 'Unsupported accepted resamples')
 legend(fig,handles,legend_labels,y=.005);save(fig,'S6','A, B',['bootstrap_operations.csv'])

 # Parameter labels are also bound, separately from outcome coordinates.
 groups={'1':{ 'A':core_ids(.25),'B':core_ids(.75)},'2':{'A':core_ids(.25),'B':core_ids(.75),'C':core_ids(.25),'D':core_ids(.75)},'3':{'A':[112,104,113],'B':[114,111,115]},'S2':{'A':[112,104,113],'B':[114,111,115],'C':[112,104,113],'D':[114,111,115]},'S3':{'A':[113],'B':[113]},'S4':{'A':list(range(116,123)),'B':list(range(116,123))},'S5':{'A':[123,124],'B':[123,124]},'S6':{'A':[1,102,104,111,119,122],'B':[1,102,104,111,119,122]}}
 for figure,parts in groups.items():
  for part,ids in parts.items():
   for cid in ids:
    bind(figure,part,'scenario_definitions.json',[{'cell':cid,**cells[cid]}],['n','tau','death_first_probability','lambda_P','lambda_W','censor_rate','base_censor_rate','a','weibull_shape_P','weibull_shape_W'])
 for part in ['A','B']:bind('4',part,'scenario_definitions.json',[definition['medical_demo']],['N','tau_days','late_deaths'])
 for file,rr in [('figure_sources.csv',trace),('figure_files.csv',inventory)]:
  with (output/file).open('w',encoding='utf-8',newline='') as dest:
   writer=csv.DictWriter(dest,fieldnames=list(rr[0]));writer.writeheader();writer.writerows(rr)
 (output/'layout_checks.json').write_text(json.dumps(layout,indent=2),encoding='utf-8')
 return inventory,trace,layout


def reproduce_tables(data,output):
 output=Path(output)/'tables';output.mkdir(parents=True,exist_ok=True)
 for table in json.loads((Path(data)/'main_table_metadata.json').read_text(encoding='utf-8')):
  with (Path(data)/table['csv']).open(encoding='utf-8',newline='') as f:rr=list(csv.reader(f))
  text='### Table %s. %s\n\n'%(table['table'],table['title'])
  text+='| '+' | '.join(rr[0])+' |\n| '+' | '.join(['---']*len(rr[0]))+' |\n'
  text+=''.join('| '+' | '.join(row)+' |\n' for row in rr[1:])
  text+='\n'+table['note']+'\n'
  (output/('Table_%s.md'%table['table'])).write_text(text,encoding='utf-8')
  with (output/('Table_%s.csv'%table['table'])).open('w',encoding='utf-8',newline='') as f:csv.writer(f).writerows(rr)

if __name__=='__main__':
 parser=argparse.ArgumentParser(description='Reproduce tables and figures from stored aggregates.')
 parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data')
 parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'outputs')
 a=parser.parse_args();reproduce_tables(a.data,a.output);main(a.data,a.output)
