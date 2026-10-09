#!/usr/bin/env python3
"""Safely copy a reviewed candidate into a CLEAN Git feature branch.
Default is dry run. A complete worktree backup (excluding .git) precedes writes.
Never commits, pushes, changes remotes, removes .git, or rewrites Git history.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import time
import uuid
import zipfile

OWNED_DIRS = ('R','man','inst','tests','scripts','docs','validation','.github')
OWNED_FILES = ('DESCRIPTION','NAMESPACE','README.md','NEWS.md','LICENSE',
               '.gitignore','.Rbuildignore','importABCatlas.Rproj')
OBSOLETE = ('src','README.Rmd','README.html','environment.yml','current_state.md','.README.Rmd.swp')

def git(root,*args):
    p=subprocess.run(['git','-C',str(root),*args],text=True,capture_output=True)
    if p.returncode:raise ValueError(p.stderr.strip() or 'Git command failed')
    return p.stdout.strip()

def migrate(candidate,target,apply=False):
    candidate=Path(candidate).expanduser().resolve();target=Path(target).expanduser().resolve()
    if candidate==target or candidate in target.parents or target in candidate.parents:
        raise ValueError('Candidate and target must be separate, non-nested directories')
    if not (candidate/'DESCRIPTION').is_file() or not (candidate/'inst/python/registry.json').is_file():
        raise ValueError('Candidate must be the supplied package root')
    if Path(git(target,'rev-parse','--show-toplevel')).resolve()!=target:
        raise ValueError('Target must be the exact Git repository root')
    branch=git(target,'branch','--show-current')
    if not branch or branch in ('main','master'):
        raise ValueError('Create and select a feature branch first; main/master/detached HEAD are refused')
    if git(target,'status','--porcelain','--untracked-files=all'):
        raise ValueError('Target working tree is not clean. Preserve or commit existing work first')
    sources=[candidate/n for n in OWNED_DIRS+OWNED_FILES]
    if any(not p.exists() for p in sources):raise ValueError('Candidate is incomplete')
    for tree in (candidate,target):
        if any(p.is_symlink() for p in tree.rglob('*') if '.git' not in p.relative_to(tree).parts):
            raise ValueError('Symlinks require manual review; automatic migration refused')
    backed=[p for p in target.rglob('*') if p.is_file() and '.git' not in p.relative_to(target).parts]
    if sum(p.stat().st_size for p in backed)>100*1024**2:
        raise ValueError('Worktree exceeds 100 MB. Move data caches outside it and review manually')
    print('Branch:',branch,'\nTarget:',target)
    print('Replace:',', '.join(OWNED_DIRS+OWNED_FILES))
    print('Remove obsolete, if present:',', '.join(OBSOLETE))
    if not apply:
        print('DRY RUN ONLY. Review paths and rerun with --apply.');return None
    backup=target.parent/(target.name+'-before-candidate-'+time.strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:8]+'.zip')
    with zipfile.ZipFile(backup,'x',zipfile.ZIP_DEFLATED) as z:
        for p in backed:z.write(p,p.relative_to(target))
    with zipfile.ZipFile(backup) as z:
        if z.testzip():raise ValueError('Backup integrity check failed; no replacements performed')
    print('Backup:',backup)
    for name in OWNED_DIRS+OWNED_FILES+OBSOLETE:
        dest=target/name
        if dest.is_dir():shutil.rmtree(dest)
        elif dest.exists():dest.unlink()
    for src in sources:
        dest=target/src.name
        if src.is_dir():shutil.copytree(src,dest,ignore=shutil.ignore_patterns('__pycache__','.pytest_cache','*.pyc'))
        else:shutil.copy2(src,dest)
    print('Copied. Inspect git diff, run tests, and commit manually. No push was performed.')
    return backup

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',required=True);p.add_argument('--target',required=True)
    p.add_argument('--apply',action='store_true');a=p.parse_args()
    try:migrate(a.candidate,a.target,a.apply)
    except (ValueError,OSError) as e:p.exit(1,'ERROR: '+str(e)+'\n')

if __name__=='__main__':main()
