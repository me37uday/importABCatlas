"""Testable ABC Atlas adapter; scientific logic is independent of R and network.

The upstream SDK owns manifests, file downloads and hash validation. This module
owns explicit table joins, exact identifiers, H5AD subsetting and provenance.
No expression file is downloaded without allow_downloads=True. No network calls
or environment installation happen at import time. See docs/ARCHITECTURE.md.
"""
from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import sys
import time
import traceback
import warnings

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

VERSION = '0.2.0.9000'
DEFAULT_RELEASE = '20260711'
REGISTRY = json.loads(Path(__file__).with_name('registry.json').read_text())
KINDS = {'metadata': 'list_metadata_files',
         'expression_matrices': 'list_expression_matrix_files',
         'image_volumes': 'list_image_volume_files',
         'mapmycells': 'list_mapmycells_files'}
KEY_PATTERN = re.compile(r'(^cell_label$|label$|identifier$|alias$|barcode$|_id$|^label$|^exp_component_name$)')


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def positive_int(x, name):
    if isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x < 1:
        raise ValueError(f'{name} must be a positive integer.')
    return int(x)


def exact_directory(requested, directories):
    """Resolve case only against an existing manifest; never invent a path."""
    if requested in directories:
        return requested
    matches = [d for d in directories if d.casefold() == requested.casefold()]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(f'Directory {requested!r} absent or ambiguous. Use list_datasets(live=TRUE).')


def create_cache(download_base, release=DEFAULT_RELEASE, offline=False):
    try:
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
    except ImportError as e:
        raise ImportError('Install the Python dependencies with setup_environment(install=TRUE).') from e
    base = Path(download_base).expanduser().resolve()
    if offline and not base.exists():
        raise FileNotFoundError(f'Offline cache does not exist: {base}')
    factory = AbcProjectCache.from_local_cache if offline else AbcProjectCache.from_s3_cache
    cache = factory(base)
    if release == 'latest':
        if offline:
            raise ValueError('Offline mode requires an explicit release, not latest.')
        cache.load_latest_manifest()
    elif re.fullmatch(r'\d{8}', str(release)):
        cache.load_manifest(f'releases/{release}/manifest.json')
    elif re.fullmatch(r'releases/\d{8}/manifest\.json', str(release)):
        cache.load_manifest(release)
    else:
        raise ValueError('release must be YYYYMMDD, a releases/YYYYMMDD/manifest.json path, or latest.')
    return cache


def files(cache, directory, kind='metadata'):
    fn = getattr(cache, KINDS[kind], None)
    if fn is None:
        raise RuntimeError(f'Upstream SDK lacks {KINDS[kind]}; update abc_atlas_access.')
    # Allen's SDK raises when an existing directory has an empty optional
    # sub-directory. Treat only that documented/observed condition as empty;
    # propagate unrelated ValueError/KeyError failures so network/manifest bugs
    # are never hidden.
    try:
        return list(fn(directory))
    except (KeyError, ValueError) as error:
        message = str(error).lower()
        if ('sub-directory is empty' in message or
                ('no ' in message and ' files found in directory ' in message)):
            return []
        raise


def resolve_spec(cache, dataset, overrides=None):
    overrides = overrides or {}
    if dataset in REGISTRY:
        spec = dict(REGISTRY[dataset])
    else:
        spec = dict(directory=dataset,cell_file='cell_metadata',gene_directory=dataset,
                    gene_file='gene',taxonomy=None,spatial=False,expression_directories=[dataset])
    allowed = {'directory','cell_file','gene_directory','gene_file','taxonomy',
               'expression_directories','spatial','file_selector','imputed','mapping_directory'}
    if set(overrides)-allowed:
        raise ValueError(f'Unknown override fields: {sorted(set(overrides)-allowed)}')
    spec.update(overrides)
    directories = list(cache.list_directories)
    for name in ('directory','gene_directory','taxonomy','mapping_directory'):
        if spec.get(name):
            spec[name] = exact_directory(spec[name], directories)
    spec['expression_directories'] = [exact_directory(x,directories)
        for x in as_list(spec['expression_directories']) if any(d.casefold()==x.casefold() for d in directories)]
    if not spec['expression_directories']:
        spec['expression_directories'] = [spec['directory']]
    selector = spec.get('file_selector')
    if selector:
        for name, directory_name in [('cell_file','directory'),('gene_file','gene_directory')]:
            if name in overrides:
                continue
            available = files(cache,spec[directory_name])
            matches = [x for x in available if selector in x.lower() and
                       spec[name] in x.lower()]
            if len(matches) != 1:
                raise ValueError(f'{dataset}: cannot uniquely discover {name} for {selector}. '
                                 f'Available files: {available}. Supply explicit overrides.')
            spec[name] = matches[0]
    return spec


def get_path(cache, directory, file_name):
    path = cache.get_file_path(directory=directory,file_name=file_name)
    if isinstance(path,dict):
        path = path['local_path']
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'File unavailable in this cache: {directory}/{file_name}')
    return path


def read_csv(path, wanted=None, preview_n=None):
    """Preserve huge numeric-looking identifiers and optionally stream by key."""
    # MapMyCells CSVs can start with comment/provenance lines. Skip only a
    # leading block; a global comment='#' would corrupt unquoted hex colors.
    leading=0
    with open(path,encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('#'):
                leading+=1
            else:
                break
    header = pd.read_csv(path,nrows=0,skiprows=leading)
    types = {c:'string' for c in header.columns if KEY_PATTERN.search(c)}
    kw = dict(dtype=types,low_memory=False,keep_default_na=False,na_values=[''],skiprows=leading)
    if preview_n is not None:
        positive_int(preview_n,'preview_n')
    if not wanted:
        return pd.read_csv(path,nrows=preview_n,**kw)
    missing = set(wanted)-set(header.columns)
    if missing:
        raise ValueError(f'Missing selection columns in {path}: {missing}')
    wanted = {k:set(map(str,as_list(v))) for k,v in wanted.items()}
    kept=[]
    with pd.read_csv(path,chunksize=100000,**kw) as reader:
        for chunk in reader:
            pred=np.ones(len(chunk),dtype=bool)
            for key, vals in wanted.items():
                pred &= chunk[key].astype('string').isin(vals).to_numpy(dtype=bool,na_value=False)
            kept.append(chunk.loc[pred])
    return pd.concat(kept,ignore_index=True) if kept else header


def table(cache,directory,file_name,wanted=None,preview_n=None):
    return read_csv(get_path(cache,directory,file_name),wanted,preview_n)


def validate_ids(frame, key, context):
    if key not in frame:
        raise ValueError(f'{context}: required column {key!r} missing.')
    if frame[key].isna().any() or frame[key].astype(str).eq('').any():
        raise ValueError(f'{context}: missing {key}.')
    frame[key] = frame[key].astype('string')
    if frame[key].duplicated().any():
        examples=frame.loc[frame[key].duplicated(),key].head(3).tolist()
        raise ValueError(f'{context}: duplicate {key}: {examples}; no silent aggregation allowed.')
    return frame


def safe_join(left,right,key,suffix):
    """A left, many-to-one join; retain order, cell count and collisions."""
    if key not in left or key not in right:
        raise ValueError(f'{suffix}: join key {key} missing.')
    left=left.copy();right=right.copy()
    left[key]=left[key].astype('string');right[key]=right[key].astype('string')
    right=right.drop_duplicates()
    if right[key].isna().any() or right[key].duplicated().any():
        raise ValueError(f'{suffix}: right-hand {key} is not a non-missing unique key; join would multiply cells.')
    overlap=(set(left)&set(right))-{key}
    rename={c:f'{c}__{suffix}' for c in overlap}
    if set(rename.values()) & (set(left)|set(right)):
        raise ValueError(f'{suffix}: ambiguous annotation column suffix.')
    right=right.rename(columns=rename)
    out=left.merge(right,on=key,how='left',sort=False,validate='many_to_one')
    if len(out)!=len(left):
        raise AssertionError('Join changed the number of cells.')
    return out


def taxonomy_pivot(cache,directory):
    names=files(cache,directory)
    pivot='cluster_to_cluster_annotation_membership_pivoted'
    if pivot in names:
        result=table(cache,directory,pivot)
        return validate_ids(result,'cluster_alias',directory)
    long_name='cluster_to_cluster_annotation_membership'
    if long_name not in names:
        raise ValueError(f'{directory}: taxonomy membership file missing.')
    m=table(cache,directory,long_name)
    if 'cluster_annotation_term_name' not in m or 'cluster_annotation_term_set_name' not in m:
        if 'cluster_annotation_term' not in names:
            raise ValueError(f'{directory}: term names absent and no term table available.')
        term=table(cache,directory,'cluster_annotation_term')
        term=term.rename(columns={'label':'cluster_annotation_term_label','name':'cluster_annotation_term_name'})
        m=safe_join(m,term,'cluster_annotation_term_label','term')
        if 'cluster_annotation_term_set_name' not in m:
            sets=table(cache,directory,'cluster_annotation_term_set')
            sets=sets.rename(columns={'label':'cluster_annotation_term_set_label','name':'cluster_annotation_term_set_name'})
            m=safe_join(m,sets,'cluster_annotation_term_set_label','termset')
    keys=['cluster_alias','cluster_annotation_term_set_name']
    for key in keys+['cluster_annotation_term_name']:
        if key not in m:
            raise ValueError(f'{directory}: missing taxonomy field {key}.')
    m['cluster_alias']=m['cluster_alias'].astype('string')
    if m[keys].isna().any().any():
        raise ValueError(f'{directory}: missing taxonomy keys.')
    conflicts=m.groupby(keys)['cluster_annotation_term_name'].nunique(dropna=False)
    if (conflicts>1).any():
        raise ValueError(f'{directory}: conflicting cluster annotations.')
    result=m.pivot_table(index='cluster_alias',columns=keys[1],values='cluster_annotation_term_name',aggfunc='first')
    result.columns.name=None
    return result.reset_index()


def add_annotations(cache,spec,cell,used):
    d=spec['directory']; names=files(cache,d)
    # Preserve all donor/disease variables; no clinical interpretation is inferred.
    for fn,key in [('library','library_label'),('sample','sample_label'),('donor','donor_label'),
                   ('disease','donor_label'),('specimen_metadata','brain_section_label')]:
        if fn in names and key in cell:
            t=table(cache,d,fn,wanted={key:cell[key].dropna().astype(str).unique().tolist()})
            cell=safe_join(cell,t,key,fn);used.append([d,fn])
    if spec.get('special')=='aging':
        # Native Aging annotations are authoritative for cluster identity/age bias.
        # The separate WMB-taxonomy directory provides reference mappings.
        ad=spec.get('mapping_directory') or exact_directory('Zeng-Aging-Mouse-WMB-taxonomy',cache.list_directories)
        for td,fn,key in [(d,'cell_cluster_annotations','cell_label'),
                          (ad,'cell_cross_mapping_annotations','cell_label')]:
            if fn in files(cache,td) and key in cell:
                t=table(cache,td,fn,wanted={key:cell[key].dropna().astype(str).unique().tolist()})
                cell=safe_join(cell,t,key,fn);used.append([td,fn])
        return cell
    if spec.get('special')=='pmdbs':
        # Native PMDBS cell metadata and its own human gene table are primary.
        # SEA-AD/WHB MapMyCells tables are optional large reference mappings.
        return cell
    tax=spec.get('taxonomy')
    if tax:
        membership_dir=next((td for td in [d,tax] if 'cell_to_cluster_membership' in files(cache,td)),None)
        if membership_dir:
            t=table(cache,membership_dir,'cell_to_cluster_membership',wanted={'cell_label':cell.cell_label.tolist()})
            # Membership from the selected taxonomy is authoritative, but retain
            # any original alias under a distinct provenance-bearing name.
            if 'cluster_alias' in cell:
                cell=cell.rename(columns={'cluster_alias':'cluster_alias_original'})
            cell=safe_join(cell,t,'cell_label','membership');used.append([membership_dir,'cell_to_cluster_membership'])
        if 'cluster_alias' not in cell:
            raise ValueError(f'{d}: no cluster_alias or applicable membership table. '
                             'Use annotations=FALSE to inspect metadata; do not invent assignments.')
        t=taxonomy_pivot(cache,tax)
        cell=safe_join(cell,t,'cluster_alias','taxonomy');used.append([tax,'taxonomy annotations'])
    if spec.get('spatial'):
        for fn in ('slab_plane_coordinates','cell_anatomical_annotations'):
            if fn in names:
                t=table(cache,d,fn,wanted={'cell_label':cell.cell_label.tolist()})
                cell=safe_join(cell,t,'cell_label',fn);used.append([d,fn])
    return cell


def filter_cells(cell,filters=None,cell_ids=None,n_cells=None,seed=1):
    out=cell.copy()
    for column, values in (filters or {}).items():
        if column not in out:
            raise ValueError(f'Filter column {column!r} absent. Inspect names(cell_metadata).')
        vals=[str(v) for v in as_list(values)]
        out=out.loc[out[column].astype('string').isin(vals)]
    if cell_ids is not None:
        ids=[str(x) for x in as_list(cell_ids)]
        if len(set(ids)) != len(ids):
            raise ValueError('cell_ids contains duplicates.')
        indexed=out.set_index('cell_label',drop=False)
        missing=set(ids)-set(indexed.index)
        if missing:
            raise ValueError(f'{len(missing)} requested cells absent after filtering: {sorted(missing)[:3]}')
        out=indexed.loc[ids].reset_index(drop=True)
    if out.empty:
        raise ValueError('No cells left after filtering; inspect values and case.')
    if n_cells is not None:
        positive_int(n_cells,'n_cells')
        if n_cells>len(out):
            raise ValueError(f'Requested {n_cells} cells, but only {len(out)} qualify.')
        # Sort first for stable input-independent sampling; record IDs in provenance.
        out=out.sort_values('cell_label',kind='stable').sample(n=n_cells,random_state=int(seed))
    return out.reset_index(drop=True)


def load_data(cache,dataset,filters=None,cell_ids=None,n_cells=None,seed=1,
              preview_n=None,annotations=True,overrides=None):
    if preview_n is not None and cell_ids is not None:
        raise ValueError('preview_n and cell_ids cannot be combined.')
    spec=resolve_spec(cache,dataset,overrides)
    used=[[spec['directory'],spec['cell_file']]]
    wanted={'cell_label':as_list(cell_ids)} if cell_ids is not None else None
    cell=table(cache,spec['directory'],spec['cell_file'],wanted=wanted,preview_n=preview_n)
    cell=validate_ids(cell,'cell_label',dataset)
    if annotations:
        cell=add_annotations(cache,spec,cell,used)
    cell=filter_cells(cell,filters,cell_ids,n_cells,seed)
    gene=table(cache,spec['gene_directory'],spec['gene_file'])
    gene=validate_ids(gene,'gene_identifier',spec['gene_directory'])
    if 'gene_symbol' not in gene:
        gene['gene_symbol']=gene.gene_identifier
    used.append([spec['gene_directory'],spec['gene_file']])
    provenance={'package_version':VERSION,'dataset':dataset,'spec':spec,
                'manifest':str(cache.current_manifest),'metadata_files':used,
                'preview_only':preview_n is not None,'n_cells':len(cell),
                'selection_seed':int(seed),'filters':filters or {},
                'cell_id_sha256':hashlib.sha256('\n'.join(cell.cell_label).encode()).hexdigest()}
    if 'gene_universe_note' in spec:
        provenance['gene_universe_note']=spec['gene_universe_note']
    return cell,gene,provenance


def resolve_genes(gene,genes=None):
    validate_ids(gene,'gene_identifier','gene metadata')
    if genes is None:
        return gene.copy().reset_index(drop=True)
    query=list(map(str,as_list(genes)))
    if not query or len(set(query))!=len(query):
        raise ValueError('genes must be a nonempty vector without duplicate requests.')
    ids=set(gene.gene_identifier.astype(str));selected=[]
    for q in query:
        if q in ids:
            matches=gene.index[gene.gene_identifier.astype(str).eq(q)].tolist()
        else:
            matches=gene.index[gene.gene_symbol.astype(str).eq(q)].tolist()
        if len(matches)!=1:
            raise ValueError(f'Gene {q!r} matches {len(matches)} identifiers; use a unique gene_identifier. '
                             'Absent genes are not replaced by zeros.')
        selected.append(matches[0])
    out=gene.loc[selected].reset_index(drop=True)
    if out.gene_identifier.duplicated().any():
        raise ValueError('Multiple requests resolve to the same gene identifier.')
    return out


def expression_plan(cache,spec,cell,data_type='raw',file_map=None):
    if data_type not in ('raw','log2'):
        raise ValueError('data_type must be raw or log2.')
    if spec.get('imputed') and data_type!='log2':
        raise ValueError('Imputed MERFISH has predicted log2 expression, not measured raw counts. Use data_type="log2".')
    known=list(cache.list_directories)
    # Expression routing is dataset-local. Do not enumerate unrelated atlas
    # directories: taxonomy/metadata-only directories legitimately have no
    # expression_matrices and must not be able to break (e.g.) a MERFISH plan.
    preferred=[exact_directory(d,known) for d in spec['expression_directories']]
    inventory=[]
    for d in preferred:
        for name in files(cache,d,'expression_matrices'):
            if name.endswith('/'+data_type):
                inventory.append((d,name))
    labels = cell['feature_matrix_label'].astype('string') if 'feature_matrix_label' in cell else pd.Series([spec['directory']]*len(cell))
    if labels.isna().any():
        raise ValueError('Missing feature_matrix_label in selected cells.')
    plans=[]
    for label in labels.unique():
        ids=cell.loc[labels.to_numpy()==label,'cell_label'].tolist()
        if file_map and label in file_map:
            item=file_map[label]
            mapped_directory=exact_directory(item['directory'],known)
            pair=(mapped_directory,item['file_name'])
            mapped_files=files(cache,mapped_directory,'expression_matrices')
            if item['file_name'] not in mapped_files:
                raise ValueError(f'Explicit file_map entry is not in the manifest: {pair}')
        else:
            exact=[x for x in inventory if x[1]==f'{label}/{data_type}']
            preferred_exact=[x for x in exact if x[0] in preferred]
            matches=preferred_exact
            if spec.get('imputed'):
                matches=[x for x in inventory if x[0] in preferred]
            if not matches:
                # Safe single-file fallback only within this dataset's known expression directories.
                matches=[x for x in inventory if x[0] in preferred]
            if len(matches)!=1:
                raise ValueError(f'Cannot uniquely route feature_matrix_label={label!r}: {matches[:8]}. '
                                 'Inspect list_files(); supply file_map, rather than guessing a path.')
            pair=matches[0]
        plans.append({'directory':pair[0],'file_name':pair[1],'feature_matrix_label':str(label),
                      'n_cells':len(ids),'cell_ids':ids})
    return plans



def _manifest_tree(cache):
    """Return SDK manifest data if exposed, without network/download operations."""
    for name in ('_manifest', 'manifest', '_manifest_data', 'manifest_data'):
        value = getattr(cache, name, None)
        if isinstance(value, dict):
            return value
    return None


def _asset_file_info(cache, directory, file_name):
    """Best-effort exact-file metadata. Never use directory totals as file size.

    Allen manifests may omit ContentLength. Unknown is preferable to an
    incorrect directory-wide estimate. No get_data_path() call is made.
    """
    result = {'file_bytes': None, 'size_source': 'unknown',
              'cached': None, 'download_required_bytes': None}
    tree = _manifest_tree(cache)
    # Manifest entries normally contain a relative_path; match the complete
    # expression path, not just an ambiguous basename.
    token = f'/{directory}/'
    name_parts = file_name.split('/')
    def walk(node):
        if isinstance(node, dict):
            path = node.get('relative_path')
            if isinstance(path, str) and token in ('/' + path) and                all(part in path for part in name_parts):
                yield node
            for value in node.values():
                yield from walk(value)
        elif isinstance(node, list):
            for value in node:
                yield from walk(value)
    entries = list(walk(tree)) if tree is not None else []
    # Do not report an ambiguous entry as an exact file.
    if len(entries) == 1:
        entry = entries[0]
        rel = entry['relative_path']
        # The SDK's cache root is implementation-specific. Only use a local
        # path when its root is explicitly exposed; never trigger a download.
        for attr in ('cache_dir', '_cache_dir', 'base_dir', 'root'):
            root = getattr(cache, attr, None)
            if isinstance(root, (str, Path)):
                candidate = Path(root).expanduser() / rel
                if candidate.is_file():
                    result.update(file_bytes=candidate.stat().st_size,
                                  size_source='local_file', cached=True,
                                  download_required_bytes=0)
                    return result
        for key in ('size_bytes', 'content_length', 'ContentLength', 'file_size_bytes'):
            value = entry.get(key)
            if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
                result.update(file_bytes=value, size_source='manifest')
                break
        result['relative_path'] = rel
    return result


def make_fetch_plan(cache, spec, cell, gene, data_type='raw', file_map=None):
    plans = expression_plan(cache, spec, cell, data_type, file_map)
    public = []
    for item in plans:
        f = {key: value for key, value in item.items() if key != 'cell_ids'}
        f.update(_asset_file_info(cache, item['directory'], item['file_name']))
        public.append(f)
    # Count each source file once even if several cell groups point to it.
    unique = {(f['directory'], f['file_name']): f for f in public}
    known = all(f['file_bytes'] is not None for f in unique.values())
    total = sum(f['file_bytes'] for f in unique.values()) if known else None
    required_known = all(f['download_required_bytes'] is not None for f in unique.values())
    required = sum(f['download_required_bytes'] for f in unique.values()) if required_known else None
    return public, total, required

def h5_strings(dataset):
    """Decode H5AD strings or categorical arrays without numeric coercion."""
    if isinstance(dataset,h5py.Group):
        if 'categories' in dataset and 'codes' in dataset:
            categories=h5_strings(dataset['categories']);codes=np.asarray(dataset['codes'])
            return np.asarray([categories[i] if i>=0 else '' for i in codes],dtype=object)
        raise ValueError('Unsupported H5AD string encoding.')
    values=np.asarray(dataset)
    return np.asarray([x.decode('utf-8') if isinstance(x,(bytes,np.bytes_)) else str(x) for x in values],dtype=object)


def index_node(group):
    key=group.attrs.get('_index','_index')
    if isinstance(key,bytes):
        key=key.decode()
    if key not in group:
        if 'index' in group:
            key='index'
        else:
            raise ValueError('H5AD dataframe index is missing.')
    return group[key]


def locate_obs(obs,cell,override=None):
    names=h5_strings(index_node(obs))
    if len(set(names))!=len(names):
        raise ValueError('H5AD obs names are not unique.')
    lookup=pd.Index(names)
    columns=[override] if override else ['cell_label','exp_component_name']
    for col in columns:
        if col and col in cell and not cell[col].isna().any():
            wanted=cell[col].astype(str).tolist()
            if len(set(wanted))!=len(wanted):
                continue
            pos=lookup.get_indexer(wanted)
            if np.all(pos>=0):
                return pos,col
    # Some files store cell_label as a column rather than the dataframe index.
    if not override and 'cell_label' in obs:
        labels=h5_strings(obs['cell_label'])
        if len(set(labels))==len(labels):
            pos=pd.Index(labels).get_indexer(cell.cell_label.astype(str))
            if np.all(pos>=0):
                return pos,'obs/cell_label'
    raise ValueError('Not all selected cells can be matched exactly in H5AD. '
                     'cell_label and exp_component_name were checked; use obs_column for an explicit key. '
                     'Barcode suffixes are never stripped.')


def subset_h5ad(path,cell,gene_ids,obs_column=None,max_output_mb=1024):
    """Read only requested rows/columns of dense, CSR or CSC H5AD X.

    HDF5 datasets are read with bounded chunks; the full expression matrix is
    never materialized. Index vectors themselves must still fit in memory.
    """
    if not np.isfinite(float(max_output_mb)) or float(max_output_mb)<=0:
        raise ValueError('max_output_mb must be positive and finite.')
    rows=len(cell);cols=len(gene_ids)
    # Conservative dense-equivalent guard also bounds pathological dense input.
    if rows*cols*8 > float(max_output_mb)*1024**2:
        raise MemoryError('Requested subset exceeds max_output_mb (dense-equivalent guard). Reduce cells/genes or raise it explicitly.')
    with h5py.File(path,'r') as f:
        for key in ('obs','var','X'):
            if key not in f:
                raise ValueError(f'{path}: H5AD lacks {key}.')
        ri,matched_by=locate_obs(f['obs'],cell,obs_column)
        names=h5_strings(index_node(f['var']))
        if len(set(names))!=len(names):
            raise ValueError('H5AD gene identifiers are duplicated.')
        ci=pd.Index(names).get_indexer(list(map(str,gene_ids)))
        if np.any(ci<0):
            missing=np.asarray(gene_ids)[ci<0].tolist()
            raise ValueError(f'{len(missing)} selected genes are unmeasured/absent in this H5AD: {missing[:5]}')
        x=f['X'];rr=[];cc=[];vv=[]
        if isinstance(x,h5py.Dataset):
            if x.shape!=(len(h5_strings(index_node(f['obs']))),len(names)):
                raise ValueError('H5AD X shape does not agree with obs/var.')
            # HDF5 fancy indices must be increasing; restore requested order.
            order=np.argsort(ci);sorted_ci=ci[order];inverse=np.argsort(order)
            for out_row,r in enumerate(ri):
                val=np.asarray(x[int(r),sorted_ci])[inverse].astype(np.float64)
                nz=np.flatnonzero(val)
                rr.extend([out_row]*len(nz));cc.extend(nz);vv.extend(val[nz])
        else:
            encoding=x.attrs.get('encoding-type','')
            if isinstance(encoding,bytes):encoding=encoding.decode()
            shape=tuple(x.attrs.get('shape',()))
            if shape and shape!=(len(h5_strings(index_node(f['obs']))),len(names)):
                raise ValueError('Sparse X shape disagrees with obs/var.')
            if encoding=='csr_matrix':
                remap={int(c):i for i,c in enumerate(ci)}
                for out_row,r in enumerate(ri):
                    start,end=np.asarray(x['indptr'][int(r):int(r)+2],dtype=np.int64)
                    for at in range(int(start),int(end),65536):
                        idx=np.asarray(x['indices'][at:min(at+65536,int(end))],dtype=np.int64)
                        val=np.asarray(x['data'][at:min(at+65536,int(end))],dtype=np.float64)
                        for j,v in zip(idx,val):
                            if int(j) in remap and v!=0:
                                rr.append(out_row);cc.append(remap[int(j)]);vv.append(v)
            elif encoding=='csc_matrix':
                remap={int(r):i for i,r in enumerate(ri)}
                for out_col,c in enumerate(ci):
                    start,end=np.asarray(x['indptr'][int(c):int(c)+2],dtype=np.int64)
                    for at in range(int(start),int(end),65536):
                        idx=np.asarray(x['indices'][at:min(at+65536,int(end))],dtype=np.int64)
                        val=np.asarray(x['data'][at:min(at+65536,int(end))],dtype=np.float64)
                        for j,v in zip(idx,val):
                            if int(j) in remap and v!=0:
                                rr.append(remap[int(j)]);cc.append(out_col);vv.append(v)
            else:
                raise ValueError(f'Unsupported sparse H5AD encoding: {encoding!r}.')
        values=np.asarray(vv,dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError('Selected expression contains NaN or infinite values.')
        matrix=sparse.coo_matrix((values,(rr,cc)),shape=(rows,cols)).tocsr()
        return matrix,matched_by


def fetch_data(cache,spec,cell,gene,genes=None,data_type='raw',allow_downloads=False,
               max_output_mb=1024,file_map=None,obs_column=None,counts_semantics=None,verify_reference=False):
    cell=validate_ids(cell.copy(),'cell_label','fetch_data metadata')
    selected=resolve_genes(gene.copy(),genes)
    if not allow_downloads and not getattr(cache,'_local',False):
        raise PermissionError('Expression downloads require allow_downloads=TRUE. Use plan_fetch() first; '
                              'even a few cells can require a complete multi-GB source file. '
                              'For a populated cache, use offline=TRUE.')
    plans=expression_plan(cache,spec,cell,data_type,file_map)
    if not np.isfinite(float(max_output_mb)) or float(max_output_mb)<=0:
        raise ValueError('max_output_mb must be finite and positive.')
    if len(cell)*len(selected)*8>float(max_output_mb)*1024**2:
        raise MemoryError('Requested subset exceeds max_output_mb before downloading.')
    parts=[];ordered_ids=[];sources=[]
    if spec.get('imputed') and genes is None:
        if len(plans)!=1:
            raise ValueError('Imputed all-gene selection requires one unambiguous source file.')
        path=get_path(cache,plans[0]['directory'],plans[0]['file_name'])
        with h5py.File(path,'r') as f:
            actual=h5_strings(index_node(f['var']))
        # The imputed H5AD var axis is authoritative, not the WMB candidate list.
        selected=resolve_genes(gene.copy(),actual.tolist())
    for p in plans:
        sub=cell.set_index('cell_label',drop=False).loc[p['cell_ids']].reset_index(drop=True)
        path=get_path(cache,p['directory'],p['file_name'])
        m,matched_by=subset_h5ad(path,sub,selected.gene_identifier.tolist(),obs_column,max_output_mb)
        if verify_reference:
            verify_with_anndata(path,sub,selected.gene_identifier.tolist(),m,matched_by)
        parts.append(m);ordered_ids.extend(sub.cell_label.tolist())
        sources.append({k:v for k,v in p.items() if k!='cell_ids'} |
                       {'local_path':str(path),'file_bytes':path.stat().st_size,'matched_by':matched_by})
    joined=sparse.vstack(parts,format='csr')
    order=pd.Index(ordered_ids).get_indexer(cell.cell_label)
    if np.any(order<0) or len(set(ordered_ids))!=len(cell):
        raise AssertionError('Expression/cell alignment invariant failed.')
    matrix=joined[order,:].T.tocsc()
    # An upstream file called raw is not by itself proof of integer UMI counts.
    integer_nonnegative=bool(np.all(matrix.data>=0) and np.allclose(matrix.data,np.rint(matrix.data),rtol=0,atol=1e-8))
    if counts_semantics is None:
        is_counts=data_type=='raw' and not spec.get('imputed',False) and integer_nonnegative
    else:
        is_counts=bool(counts_semantics)
        if is_counts and (data_type!='raw' or spec.get('imputed') or not integer_nonnegative):
            raise ValueError('Cannot label transformed, imputed, negative or noninteger data as counts.')
    assay_kind='counts' if is_counts else ('imputed_log2' if spec.get('imputed') else data_type)
    provenance={'package_version':VERSION,'manifest':str(cache.current_manifest),'sources':sources,
                'data_type':data_type,'assay_kind':assay_kind,'integer_nonnegative':integer_nonnegative,
                'reference_comparison':bool(verify_reference),
                'n_cells':matrix.shape[1],'n_genes':matrix.shape[0],
                'cell_id_sha256':hashlib.sha256('\n'.join(cell.cell_label).encode()).hexdigest(),
                'gene_id_sha256':hashlib.sha256('\n'.join(selected.gene_identifier).encode()).hexdigest()}
    return matrix,cell,selected,provenance


def verify_with_anndata(path,cell,gene_ids,actual,matched_by):
    """Independent source-value comparison, intended for small live smoke tests."""
    try:
        import anndata
    except ImportError as e:
        raise ImportError('anndata is required for verify_reference=True.') from e
    ad=anndata.read_h5ad(path,backed='r')
    try:
        axis=pd.Index(ad.obs['cell_label'].astype(str)) if matched_by=='obs/cell_label' else pd.Index(ad.obs_names.astype(str))
        query=cell['cell_label' if matched_by=='obs/cell_label' else matched_by].astype(str)
        ri=axis.get_indexer(query);ci=pd.Index(ad.var_names.astype(str)).get_indexer(gene_ids)
        if np.any(ri<0) or np.any(ci<0):
            raise AssertionError('Independent AnnData axis lookup disagrees with importer.')
        for i,row_index in enumerate(ri):
            row=ad[int(row_index),:].X
            expected=(row.toarray() if sparse.issparse(row) else np.asarray(row)).reshape(-1)[ci]
            np.testing.assert_allclose(actual.getrow(i).toarray().reshape(-1),expected,rtol=0,atol=0,
                                       err_msg='Independent AnnData source-value comparison failed.')
    finally:
        ad.file.close()


def environment():
    versions={}
    for package in ('numpy','pandas','scipy','h5py','abc_atlas_access'):
        try: versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: versions[package]=None
    try:
        dist=importlib.metadata.distribution('abc_atlas_access')
        direct=dist.read_text('direct_url.json')
        versions['abc_atlas_access_source']=json.loads(direct) if direct else None
    except importlib.metadata.PackageNotFoundError:
        pass
    return {'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'packages':versions,'package_version':VERSION}


def write_frame(frame,path):
    frame.to_csv(path,index=False)
    schema={c:('string' if KEY_PATTERN.search(c) else
               'numeric' if pd.api.types.is_numeric_dtype(frame[c]) and not pd.api.types.is_bool_dtype(frame[c])
               else 'logical' if pd.api.types.is_bool_dtype(frame[c]) else 'string') for c in frame}
    Path(str(path)+'.schema.json').write_text(json.dumps(schema))


def file_catalog(cache,directory=None):
    rows=[]
    directories=[exact_directory(directory,cache.list_directories)] if directory else cache.list_directories
    for d in directories:
        for kind in KINDS:
            for name in files(cache,d,kind):
                rows.append({'directory':d,'kind':kind,'file_name':name})
    return rows


def run_request(req):
    action=req['action'];out=Path(req['out_dir']);out.mkdir(parents=True,exist_ok=True)
    args=req.get('args',{})
    if action=='diagnose':return environment()
    if action=='registry':return REGISTRY
    if action=='demo':
        root=Path(__file__).resolve().parents[1]/'extdata'
        cell=read_csv(root/'demo_cells.csv');gene=read_csv(root/'demo_genes.csv')
        matrix,_=subset_h5ad(root/'demo.h5ad',cell,gene.gene_identifier.tolist())
        matrix=matrix.T.tocsc()
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv');mmwrite(out/'matrix.mtx',matrix)
        return {'assay_kind':'counts','synthetic_fixture':True,'n_cells':len(cell),'n_genes':len(gene)}
    cache=create_cache(args.pop('download_base','abc_download_root'),
                       args.pop('release',DEFAULT_RELEASE),args.pop('offline',False))
    if action=='catalog':return {'manifest':cache.current_manifest,'files':file_catalog(cache,args.get('directory'))}
    if action=='download_file':
        path=get_path(cache,args['directory'],args['file_name'])
        return {'path':str(path),'manifest':cache.current_manifest}
    if action=='metadata':
        frame=table(cache,args['directory'],args['file_name'],preview_n=args.get('preview_n'))
        write_frame(frame,out/'table.csv');return {'manifest':cache.current_manifest}
    if action=='load':
        cell,gene,prov=load_data(cache,**args)
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv')
        return prov
    if action in ('fetch','plan'):
        cell=read_csv(args.pop('metadata_path'));gene=read_csv(args.pop('gene_path'))
        cell=filter_cells(cell,args.pop('filters',None),args.pop('cell_ids',None),args.pop('n_cells',None),args.pop('seed',1))
        spec=resolve_spec(cache,args.pop('dataset'),args.pop('overrides',None))
        if action=='plan':
            selected=resolve_genes(gene,args.get('genes'))
            plans,total_bytes,download_bytes=make_fetch_plan(cache,spec,cell,gene,args.get('data_type','raw'),args.get('file_map'))
            return {'manifest':cache.current_manifest,'files':plans,
                    'total_file_bytes':total_bytes,'total_download_bytes':download_bytes,
                    'data_type':args.get('data_type','raw'),
                    'n_cells':len(cell),'n_genes':len(selected),'dense_equivalent_mb':len(cell)*len(selected)*8/1024**2,
                    'warning':'Exact file size is unknown unless the manifest supplies byte size or the file is locally verified. No expression files were downloaded.'}
        matrix,cell,gene,prov=fetch_data(cache,spec,cell,gene,**args)
        write_frame(cell,out/'cell_metadata.csv');write_frame(gene,out/'gene_data.csv');mmwrite(out/'matrix.mtx',matrix)
        return prov
    raise ValueError(f'Unknown action: {action}')


def main():
    if len(sys.argv)!=2:
        raise SystemExit('Usage: python abc_backend.py request.json')
    req=json.loads(Path(sys.argv[1]).read_text());out=Path(req['out_dir']);out.mkdir(parents=True,exist_ok=True)
    try:
        start=time.perf_counter();result=run_request(req)
        if isinstance(result,dict) and req['action']!='registry':result['backend_elapsed_seconds']=time.perf_counter()-start
        (out/'result.json').write_text(json.dumps(result,indent=2,default=str,allow_nan=False))
    except Exception as error:
        (out/'error.json').write_text(json.dumps({'type':type(error).__name__,'message':str(error),'traceback':traceback.format_exc()},indent=2))
        traceback.print_exc();return 1
    return 0


if __name__=='__main__':
    raise SystemExit(main())
