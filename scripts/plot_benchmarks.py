#!/usr/bin/env python3
"""Create ONE runtime figure only from genuine, paired, successful measurements.
The plot has no built-in demonstration values. Missing conditions remain blank.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def summarize(path):
    df=pd.read_csv(path)
    required={'dataset','cohort','stage','cache_state','replicate','status','seconds','data_origin','selection_sha256','return_type','data_type','manifest'}
    if not required<=set(df):raise ValueError('Missing benchmark columns')
    if not df.data_origin.eq('LIVE_ATLAS').all():raise ValueError('Synthetic or unverified-origin timings cannot be used for the manuscript figure')
    passed=df.loc[df.status.eq('PASS')].copy()
    if passed.empty:raise ValueError('No successful live measurements. No figure created.')
    if not (pd.to_numeric(passed.seconds,errors='coerce')>0).all():raise ValueError('Runtimes must be finite, positive, measured values')
    if not np.isfinite(passed.seconds).all():raise ValueError('Nonfinite runtime')
    for _,g in df.groupby(['dataset','cohort']):
        for field in ('selection_sha256','return_type','data_type','manifest'):
            if g[field].nunique()!=1:raise ValueError('Do not combine different selections, formats, transformations or manifests')
    complete=[]
    for _,g in passed.groupby(['dataset','cohort','replicate','stage']):
        if g.cache_state.duplicated().any():raise ValueError('Duplicate benchmark conditions')
        if set(g.cache_state)=={'cold','warm'}:complete.append(g)
    if not complete:raise ValueError('No complete cold/warm pairs. No figure created.')
    passed=pd.concat(complete,ignore_index=True)
    passed['condition']=passed.stage+' / '+passed.cache_state
    med=passed.pivot_table(index=['dataset','cohort'],columns='condition',values='seconds',aggfunc='median')
    med=med.reindex(pd.MultiIndex.from_frame(df[['dataset','cohort']].drop_duplicates()))
    return med.reindex(columns=['load_data / cold','load_data / warm','fetch_data / cold','fetch_data / warm']),df

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('csv');ap.add_argument('--output',default='figure1.png');a=ap.parse_args()
    med,raw=summarize(a.csv)
    height=max(4,min(15,.24*len(med)+2))
    fig,ax=plt.subplots(figsize=(8,height),layout='constrained')
    vals=med.to_numpy();finite=vals[np.isfinite(vals)]
    lo=float(finite.min());hi=max(float(finite.max()),lo*1.01)
    im=ax.imshow(np.ma.masked_invalid(vals),aspect='auto',norm=LogNorm(vmin=lo,vmax=hi))
    ax.set_xticks(range(4),['load: cold','load: warm','fetch: cold','fetch: warm'])
    ax.set_yticks(range(len(med)),[d+' | '+c.replace('_',' ') for d,c in med.index],fontsize=7)
    for i in range(len(med)):
        for j in range(4):
            if not np.isfinite(vals[i,j]):ax.text(j,i,'not measured',ha='center',va='center',fontsize=6)
    ax.set_title('Median runtime (seconds); frozen cell/gene selections',fontsize=11)
    fig.colorbar(im,ax=ax,label='Seconds (logarithmic scale)',fraction=.035)
    fig.savefig(a.output,dpi=350)
    fig.savefig(Path(a.output).with_suffix('.svg'))
    med.to_csv(Path(a.output).with_suffix('.summary.csv'))
    print('Wrote',a.output,'; inspect coverage and Figure 1 legend before submission.')

if __name__=='__main__':main()
