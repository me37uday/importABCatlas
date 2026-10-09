#!/usr/bin/env python3
"""Paired cold/warm, fresh-R-process benchmark. NEVER creates invented results.

Requires R package installed, Python env configured, psutil, reviewed frozen
profiles and explicit download consent. Every replicate has its own NEW cache;
this script never deletes a cache. Source transfers can consume substantial disk.
"""
import argparse,csv,hashlib,json,os,shutil,subprocess,time,uuid
from pathlib import Path
import psutil

FIELDS=['dataset','cohort','stage','cache_state','replicate','status','seconds',
        'sampled_peak_process_tree_rss_mb','net_cache_growth_bytes','n_cells','n_genes',
        'return_type','data_type','manifest','selection_sha256','data_origin','error']

def cache_bytes(root):
    # File sizes in cache, NOT network traffic. Sum file lengths, not disk blocks.
    return sum(p.stat().st_size for p in Path(root).rglob('*') if p.is_file()) if Path(root).exists() else 0

def validate_profile(p):
    for key in ('dataset','cohort','cell_ids','genes','release','return_type','data_type'):
        if key not in p:raise ValueError('Profile missing '+key)
    if p['cohort'] not in ('neuronal','non_neuronal'):raise ValueError('Unreviewed cohort name')
    if p.get('cohort_reviewed') is not True:raise ValueError('Cohort must be explicitly reviewed')
    for key in ('cell_ids','genes'):
        if not isinstance(p[key],list) or not p[key] or len(p[key])!=len(set(p[key])):raise ValueError(key+' must be a nonempty unique list')
        if not all(isinstance(x,str) for x in p[key]):raise ValueError(key+' must contain exact strings')
    if p['return_type'] not in ('matrix','sce','seurat','spatial'):raise ValueError('Unsupported object format')
    if p['dataset']=='HMBA_PatchSeq' and p['cohort']=='non_neuronal':raise ValueError('Patch-seq non-neuronal cohort is not applicable')

def execute(req,worker,rscript):
    out=Path(req['out_dir']);out.mkdir(parents=True,exist_ok=False)
    request=out/'request.json';request.write_text(json.dumps(req,indent=2))
    peak=0;start=time.perf_counter()
    with (out/'R.log').open('w') as log:
        proc=subprocess.Popen([rscript,str(worker),str(request)],stdout=log,stderr=subprocess.STDOUT)
        parent=psutil.Process(proc.pid)
        while proc.poll() is None:
            flag=out/'stage.txt'
            if flag.exists() and flag.read_text().strip()=='running':
                try:
                    processes=[parent]+parent.children(recursive=True)
                    rss=sum(z.memory_info().rss for z in processes if z.is_running())
                    peak=max(peak,rss)
                except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            time.sleep(.05)
    result_path=out/'result.json'
    result=json.loads(result_path.read_text()) if result_path.exists() else {'status':'FAIL','error':'Worker failed; see R.log'}
    result['sampled_peak_process_tree_rss_mb']=peak/1024**2 if peak else None
    result['process_elapsed_seconds']=time.perf_counter()-start
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profiles',default='benchmark_profiles.json');ap.add_argument('--cache-root',required=True)
    ap.add_argument('--out',default='benchmark_results');ap.add_argument('--replicates',type=int,default=3)
    ap.add_argument('--rscript',default='Rscript');ap.add_argument('--allow-downloads',action='store_true')
    a=ap.parse_args()
    if not a.allow_downloads:ap.error('Full H5AD files may be downloaded. Inspect plans, then pass --allow-downloads.')
    if a.replicates<1:ap.error('replicates must be positive')
    if not shutil.which(a.rscript):ap.error('Rscript not found')
    profiles=json.loads(Path(a.profiles).read_text())['profiles']
    if not profiles:ap.error('No reviewed profiles')
    for p in profiles:validate_profile(p)
    keys=[(p['dataset'],p['cohort']) for p in profiles]
    if len(set(keys))!=len(keys):ap.error('Duplicate dataset/cohort profiles. Run object-format comparisons separately.')
    stamp=time.strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:8]
    root=Path(a.cache_root).expanduser().resolve()/stamp;root.mkdir(parents=True,exist_ok=False)
    out=Path(a.out).resolve()/stamp;out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(a.profiles,out/'frozen_profiles.json')
    results=[];worker=Path(__file__).with_name('benchmark_worker.R')
    for number,p in enumerate(profiles):
        digest=hashlib.sha256(json.dumps({'cell_ids':p['cell_ids'],'genes':p['genes']},sort_keys=True).encode()).hexdigest()
        for rep in range(1,a.replicates+1):
            cache=root/f'profile-{number}-rep-{rep}'
            cold_ready=True
            for state in ('cold','warm'):
                rds=out/f'profile-{number}-{rep}-{state}.rds'
                previous_ok=True
                for stage in ('load_data','fetch_data'):
                    before=cache_bytes(cache)
                    req={'profile':p,'stage':stage,'cache':str(cache),'metadata_rds':str(rds),
                         'out_dir':str(out/f'profile-{number}-{rep}-{state}-{stage}')}
                    result=execute(req,worker,a.rscript) if previous_ok and (state=='cold' or cold_ready) else {'status':'NOT_RUN','error':'No complete preceding cold run, or load_data failed'}
                    if state=='cold' and result['status']!='PASS':cold_ready=False
                    previous_ok=result['status']=='PASS'
                    row={k:p.get(k,'') for k in ('dataset','cohort','return_type','data_type')}
                    row.update(stage=stage,cache_state=state,replicate=rep,status=result['status'],
                      seconds=result.get('seconds'),sampled_peak_process_tree_rss_mb=result.get('sampled_peak_process_tree_rss_mb'),
                      net_cache_growth_bytes=cache_bytes(cache)-before,n_cells=result.get('n_cells'),n_genes=result.get('n_genes'),
                      manifest=p['release'],selection_sha256=digest,data_origin='LIVE_ATLAS',error=result.get('error',''))
                    results.append(row)
                    with (out/'benchmark_results.csv').open('w',newline='') as f:
                        writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader();writer.writerows(results)
                    print(p['dataset'],p['cohort'],rep,state,stage,result['status'],flush=True)
    print('Results:',out);print('Caches were not deleted:',root)

if __name__=='__main__':main()
