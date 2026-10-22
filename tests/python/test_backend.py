import json
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmread
import abc_backend as b
from fixtures import CELLS, GENES, VALUES, h5ad, FakeCache

@pytest.mark.parametrize('encoding',['csr','csc','dense'])
@pytest.mark.parametrize('categorical',[False,True])
@pytest.mark.parametrize('order',[[2,0],[1,2,0]])
def test_h5ad_subset_exact(tmp_path,encoding,categorical,order):
    path=h5ad(tmp_path/'x.h5ad',encoding,categorical=categorical)
    cm=pd.DataFrame({'cell_label':[CELLS[i] for i in order]})
    actual,by=b.subset_h5ad(path,cm,[GENES[3],GENES[0]])
    assert sparse.issparse(actual)
    np.testing.assert_array_equal(actual.toarray(),VALUES[np.ix_(order,[3,0])])
    assert by=='cell_label'

@pytest.mark.parametrize('encoding',['csr','csc','dense'])
def test_zero_matrix(tmp_path,encoding):
    path=h5ad(tmp_path/'x.h5ad',encoding,values=np.zeros_like(VALUES))
    m,_=b.subset_h5ad(path,pd.DataFrame({'cell_label':CELLS}),GENES)
    assert m.nnz==0 and m.shape==VALUES.shape

def test_exact_exp_component_fallback(tmp_path):
    path=h5ad(tmp_path/'x.h5ad',cells=['rna-A','rna-B','rna-C'])
    cm=pd.DataFrame({'cell_label':CELLS,'exp_component_name':['rna-A','rna-B','rna-C']})
    m,by=b.subset_h5ad(path,cm,GENES)
    assert by=='exp_component_name';np.testing.assert_array_equal(m.toarray(),VALUES)

def test_obs_column_fallback(tmp_path):
    path=h5ad(tmp_path/'x.h5ad',cells=['a','b','c'],obs_column=CELLS)
    m,by=b.subset_h5ad(path,pd.DataFrame({'cell_label':CELLS}),GENES)
    np.testing.assert_array_equal(m.toarray(),VALUES)

@pytest.mark.parametrize('bad',['missing','9007199254740992'])
def test_missing_cells_not_silently_dropped(tmp_path,bad):
    path=h5ad(tmp_path/'x.h5ad')
    with pytest.raises(ValueError,match='match|absent|Missing|index'):
        b.subset_h5ad(path,pd.DataFrame({'cell_label':[bad]}),GENES)

def test_missing_genes_not_zeros(tmp_path):
    with pytest.raises(ValueError,match='gene|Gene'):
        b.subset_h5ad(h5ad(tmp_path/'x.h5ad'),pd.DataFrame({'cell_label':CELLS}),['human-gene'])

def test_duplicate_source_ids_fail(tmp_path):
    with pytest.raises(ValueError,match='duplicat|unique'):
        b.subset_h5ad(h5ad(tmp_path/'x.h5ad',cells=['a','a','c']),pd.DataFrame({'cell_label':['a']}),GENES)

@pytest.mark.parametrize('cap',[-1,0,float('inf'),float('nan')])
def test_invalid_memory_guard(tmp_path,cap):
    with pytest.raises(ValueError):b.subset_h5ad(h5ad(tmp_path/'x.h5ad'),pd.DataFrame({'cell_label':CELLS}),GENES,max_output_mb=cap)

def test_memory_guard(tmp_path):
    with pytest.raises(MemoryError):b.subset_h5ad(h5ad(tmp_path/'x.h5ad'),pd.DataFrame({'cell_label':CELLS}),GENES,max_output_mb=1e-7)

@pytest.mark.parametrize('encoding',['csr','csc','dense'])
def test_nonfinite_rejected(tmp_path,encoding):
    v=VALUES.copy();v[0,0]=float('nan')
    with pytest.raises(ValueError,match='NaN|finite'):
        b.subset_h5ad(h5ad(tmp_path/'x.h5ad',encoding,values=v),pd.DataFrame({'cell_label':CELLS}),GENES)

def test_csv_huge_ids_and_unknown_category(tmp_path):
    p=tmp_path/'a.csv';p.write_text('cell_label,donor_label,APOE Genotype,x\n'+CELLS[0]+',00023,Unknown,1.2\n')
    df=b.read_csv(p);assert df.cell_label[0]==CELLS[0] and df.donor_label[0]=='00023'
    assert df['APOE Genotype'][0]=='Unknown';assert df.x[0]==1.2
    b.write_frame(df,tmp_path/'b.csv')
    schema=json.loads((tmp_path/'b.csv.schema.json').read_text())
    assert schema['cell_label']=='string' and schema['x']=='numeric'

def test_stream_filter_preserves_ids(tmp_path):
    p=tmp_path/'a.csv';pd.DataFrame({'cell_label':CELLS,'cluster_alias':['001','003','004']}).to_csv(p,index=False)
    out=b.read_csv(p,wanted={'cell_label':[CELLS[2],CELLS[0]]})
    assert out.cell_label.tolist()==[CELLS[0],CELLS[2]]
    assert out.cluster_alias.tolist()==['001','004']

def test_join_order_and_collision():
    l=pd.DataFrame({'cell_label':['b','a','c'],'donor_label':['2','1','2'],'age':[3,4,5]})
    r=pd.DataFrame({'donor_label':['1','2'],'age':[10,20]})
    out=b.safe_join(l,r,'donor_label','donor')
    assert out.cell_label.tolist()==['b','a','c'];assert out.age.tolist()==[3,4,5]
    assert out.age__donor.tolist()==[20,10,20]

def test_many_to_many_join_blocked():
    with pytest.raises(ValueError,match='multiply'):
        b.safe_join(pd.DataFrame({'donor_label':['a','a']}),pd.DataFrame({'donor_label':['a','a'],'x':[1,2]}),'donor_label','donor')

def test_no_fake_cluster_alias(tmp_path):
    c=FakeCache(tmp_path)
    c.add('T','cluster_to_cluster_annotation_membership',{'cluster_alias':['991','991','42','42'],
      'cluster_annotation_term_set_name':['class','subclass']*2,'cluster_annotation_term_name':['GABA','A','Glia','B']})
    out=b.taxonomy_pivot(c,'T').set_index('cluster_alias')
    assert out.loc['991','class']=='GABA';assert out.loc['42','subclass']=='B';assert '1' not in out.index

def test_conflicting_taxonomy_rejected(tmp_path):
    c=FakeCache(tmp_path);c.add('T','cluster_to_cluster_annotation_membership',{
      'cluster_alias':['42','42'],'cluster_annotation_term_set_name':['class','class'],'cluster_annotation_term_name':['A','B']})
    with pytest.raises(ValueError,match='conflicting'):b.taxonomy_pivot(c,'T')

def test_filter_and_sample_reproducible():
    cm=pd.DataFrame({'cell_label':list(map(str,range(20))),'class':['A']*10+['B']*10})
    a=b.filter_cells(cm,{'class':'A'},n_cells=3,seed=42)
    z=b.filter_cells(cm.iloc[::-1],{'class':['A']},n_cells=3,seed=42)
    assert a.cell_label.tolist()==z.cell_label.tolist()
    assert b.filter_cells(cm,cell_ids=['3','1']).cell_label.tolist()==['3','1']
    with pytest.raises(ValueError):b.filter_cells(cm,{'not_column':'A'})
    with pytest.raises(ValueError):b.filter_cells(cm,{'class':'C'})

def test_gene_symbol_ambiguity():
    g=pd.DataFrame({'gene_identifier':['a','b'],'gene_symbol':['same','same']})
    with pytest.raises(ValueError,match='2 identifiers'):b.resolve_genes(g,['same'])
    assert b.resolve_genes(g,['b']).gene_identifier.tolist()==['b']

def test_case_resolved_from_manifest_only():
    assert b.exact_directory('HMBA-10X',''.split()+['HMBA-10x'])=='HMBA-10x'
    with pytest.raises(ValueError):b.exact_directory('HMBA',['hmba','HmBa'])

def route_fixture(tmp_path,dataset,encoding='csr',value_override=None):
    """Each route is exercised against a SYNTHETIC, not observed live, schema."""
    c=FakeCache(tmp_path);spec=b.REGISTRY[dataset];d=spec['directory'];g=spec['gene_directory'];t=spec['taxonomy']
    cellname=spec['cell_file'];genename=spec['gene_file']
    selector=spec.get('file_selector')
    if selector:
        # Native HMBA species share cell_metadata and have exact species gene
        # tables.  Include a tempting example-expression file to ensure the
        # resolver never mistakes it for gene metadata.
        cellname='cell_metadata';genename=f'{selector}_gene'
    c.add(d,cellname,{'cell_label':CELLS,'cluster_alias':['991','42','991'],
                     'feature_matrix_label':['matrix']*3,'donor_label':['001','002','001'],
                     'x':[1.,2.,3.],'y':[2.,3.,4.]})
    donor_data={'donor_label':['001','002'],'diagnosis':['Unknown','AD']}
    if selector:
        donor_data['species_genus']=[selector.capitalize(),selector.capitalize()]
    c.add(d,'donor',donor_data)
    c.add(g,genename,{'gene_identifier':GENES,'gene_symbol':['Gad1','Aqp4','X','Y']})
    if selector:
        c.add(g,f'example_{selector}_gene_expression',{'gene_identifier':['not-a-gene-table']})
    if t:c.add(t,'cluster_to_cluster_annotation_membership',{
      'cluster_alias':['991','42'],'cluster_annotation_term_set_name':['class']*2,'cluster_annotation_term_name':['Neuron','Non-neuron']})
    if spec.get('mapping_directory') and spec['mapping_directory'] not in c.list_directories:
        c.add(spec['mapping_directory'],'_placeholder',{'x':['1']})
    if spec.get('special')=='aging':
        c.add(d,'cell_cluster_annotations',{'cell_label':CELLS,'cluster_name':['native1','native2','native1'],'cluster_age_bias':['aged','adult','aged']})
        c.add(spec['mapping_directory'],'cell_cross_mapping_annotations',{'cell_label':CELLS,'wmb_cluster_alias':['990','41','990'],'wmb_class_name':['Neuron','Glia','Neuron']})
    for ed in spec['expression_directories']:
        p=h5ad(tmp_path/(dataset.replace('/','_')+str(len(c.assets))+'.h5ad'),encoding,values=value_override)
        kind=spec.get('default_data_type','log2' if spec.get('imputed') else 'raw')
        c.add(ed,'matrix/'+kind,p,'expression_matrices')
    return c

@pytest.mark.parametrize('dataset',list(b.REGISTRY))
def test_every_route_load_and_fetch_synthetic(tmp_path,dataset):
    c=route_fixture(tmp_path,dataset)
    cm,gd,prov=b.load_data(c,dataset,cell_ids=[CELLS[2],CELLS[0]])
    spec=prov['spec'];kind=spec.get('default_data_type','log2' if spec.get('imputed') else 'raw')
    # Multi-directory WMB fixtures deliberately contain the same label: select a
    # source explicitly just as the API requires for genuinely ambiguous input.
    fmap={'matrix':{'directory':spec['expression_directories'][0],'file_name':'matrix/'+kind}}
    m,cell,genes,p=b.fetch_data(c,spec,cm,gd,genes=[GENES[3],GENES[0]],data_type=kind,file_map=fmap)
    np.testing.assert_array_equal(m.toarray(),VALUES[np.ix_([2,0],[3,0])].T)
    assert cell.cell_label.tolist()==[CELLS[2],CELLS[0]]
    assert cell.diagnosis.tolist()==['Unknown','Unknown']
    expected=('imputed_log2' if spec.get('imputed') else ('log2' if kind=='log2' else 'counts'))
    assert p['assay_kind']==expected

@pytest.mark.parametrize('dataset',['AgingMouse','ConsensusMouseAIBS'])
def test_mouse_genes_not_human(dataset):assert b.REGISTRY[dataset]['gene_directory']=='WMB-10X'

def test_pmdbs_uses_own_gene_table_and_optional_reference_mappings():
    spec=b.REGISTRY['PMDBS']
    assert spec['gene_directory']=='ASAP-PMDBS-10X'
    assert spec['gene_file']=='gene'
    assert spec['mapping_directory']=='ASAP-PMDBS-taxonomy'

def test_aging_native_and_wmb_mapping_join(tmp_path):
    c=route_fixture(tmp_path,'AgingMouse')
    cm,gd,prov=b.load_data(c,'AgingMouse',cell_ids=[CELLS[0],CELLS[1]])
    assert cm.cluster_name.tolist()==['native1','native2']
    assert cm.wmb_cluster_alias.tolist()==['990','41']
    assert cm.wmb_class_name.tolist()==['Neuron','Glia']
    assert prov['spec']['mapping_directory']=='Zeng-Aging-Mouse-WMB-taxonomy'

def test_download_consent_before_file(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,p=b.load_data(c,'WHB');c._local=False;c.requested=[]
    with pytest.raises(PermissionError):b.fetch_data(c,p['spec'],cm,gd)
    assert c.requested==[]

def test_imputed_cannot_be_counts(tmp_path):
    c=route_fixture(tmp_path,'MERFISH_imputed');cm,gd,p=b.load_data(c,'MERFISH_imputed')
    with pytest.raises(ValueError,match='Imputed'):b.fetch_data(c,p['spec'],cm,gd,data_type='raw')
    with pytest.raises(ValueError,match='Cannot label'):b.fetch_data(c,p['spec'],cm,gd,data_type='log2',counts_semantics=True)

def test_imputed_all_genes_actual_var(tmp_path):
    c=route_fixture(tmp_path,'MERFISH_imputed');cm,gd,p=b.load_data(c,'MERFISH_imputed')
    gd=pd.concat([gd,pd.DataFrame({'gene_identifier':['not-imputed'],'gene_symbol':['other']})],ignore_index=True)
    m,_,selected,_=b.fetch_data(c,p['spec'],cm,gd,data_type='log2')
    assert m.shape==(4,3);assert 'not-imputed' not in selected.gene_identifier.tolist()

def test_noninteger_raw_not_counts(tmp_path):
    c=route_fixture(tmp_path,'WHB',value_override=VALUES*.5);cm,gd,p=b.load_data(c,'WHB')
    _,_,_,prov=b.fetch_data(c,p['spec'],cm,gd)
    assert prov['assay_kind']=='raw'

def test_ambiguous_plan_requires_mapping(tmp_path):
    c=route_fixture(tmp_path,'WMB');cm,gd,p=b.load_data(c,'WMB')
    with pytest.raises(ValueError,match='Cannot uniquely route'):b.expression_plan(c,p['spec'],cm)

def test_catalog_all_asset_kinds(tmp_path):
    c=route_fixture(tmp_path,'WHB');p=tmp_path/'dummy';p.write_bytes(b'fixture')
    c.add('CCF','image',p,'image_volumes');c.add('WHB-taxonomy','model',p,'mapmycells')
    assert {r['kind'] for r in b.file_catalog(c)}==set(b.KINDS)

def test_cli_demo_roundtrip(tmp_path):
    req={'action':'demo','args':{},'out_dir':str(tmp_path/'out')}
    p=tmp_path/'request.json';p.write_text(json.dumps(req))
    proc=subprocess.run([sys.executable,b.__file__,str(p)],capture_output=True,text=True)
    assert proc.returncode==0,proc.stderr
    result=json.loads((tmp_path/'out/result.json').read_text());assert result['synthetic_fixture']
    np.testing.assert_array_equal(mmread(tmp_path/'out/matrix.mtx').toarray(),VALUES.T)

def test_cli_error_is_machine_readable(tmp_path):
    p=tmp_path/'request.json';p.write_text(json.dumps({'action':'fetch','args':{'release':'invalid'},'out_dir':str(tmp_path/'out')}))
    proc=subprocess.run([sys.executable,b.__file__,str(p)],capture_output=True,text=True)
    assert proc.returncode!=0
    error=json.loads((tmp_path/'out/error.json').read_text());assert error['message'] and error['traceback']

@pytest.mark.skip(reason='R is absent in the delivery execution environment; R integration suite is supplied separately.')
def test_r_integration_not_executed():pass

@pytest.mark.skip(reason='Container DNS/network and upstream SDK unavailable; live atlas downloads were not executed.')
def test_live_atlas_not_executed():pass

def test_leading_comments_and_hex_colors(tmp_path):
    p=tmp_path/'mmc.csv';p.write_text('# provenance here\n# JSON header\ncell_label,color_hex_triplet\na,#FFFFFF\n')
    out=b.read_csv(p)
    assert out.cell_label.tolist()==['a'] and out.color_hex_triplet.tolist()==['#FFFFFF']

def test_override_bypasses_species_discovery(tmp_path):
    c=route_fixture(tmp_path,'HMBA_Human')
    s=b.resolve_spec(c,'HMBA_Human',{'cell_file':'cell_metadata','gene_file':'human_gene'})
    assert s['cell_file']=='cell_metadata' and s['gene_file']=='human_gene'


@pytest.mark.parametrize('dataset,selector,expected_gene',[
    ('HMBA_Human','human','human_gene'),
    ('HMBA_Macaque','macaque','macaque_gene'),
    ('HMBA_Marmoset','marmoset','marmoset_gene'),
])
def test_hmba_native_species_uses_shared_cells_and_exact_gene_table(tmp_path,dataset,selector,expected_gene):
    c=route_fixture(tmp_path,dataset)
    spec=b.resolve_spec(c,dataset)
    assert spec['cell_file']=='cell_metadata'
    assert spec['gene_file']==expected_gene
    assert spec['species_filter']==selector.capitalize()
    cm,gd,prov=b.load_data(c,dataset,preview_n=2)
    assert len(cm)==2 and len(gd)==len(GENES)
    assert set(cm.donor_label)=={'001','002'}
    assert [spec['gene_directory'],expected_gene] in prov['metadata_files']
    assert all('example_' not in file_name for _,file_name in c.requested)


def test_hmba_species_preview_filters_before_head(tmp_path):
    c=route_fixture(tmp_path,'HMBA_Human')
    # Put a non-human donor/cell first.  Preview must still return Human cells.
    cell_path=c.assets[('HMBA-10XMultiome-BG','metadata','cell_metadata')]
    original=pd.read_csv(cell_path,dtype=str)
    extra=original.iloc[[0]].copy();extra['cell_label']='other-species';extra['donor_label']='999'
    pd.concat([extra,original],ignore_index=True).to_csv(cell_path,index=False)
    donor_path=c.assets[('HMBA-10XMultiome-BG','metadata','donor')]
    donors=pd.read_csv(donor_path,dtype=str)
    pd.concat([pd.DataFrame({'donor_label':['999'],'diagnosis':['Unknown'],'species_genus':['Macaque']}),donors],ignore_index=True).to_csv(donor_path,index=False)
    cm,_,_=b.load_data(c,'HMBA_Human',preview_n=1,annotations=False)
    assert cm.cell_label.tolist()==[CELLS[0]]

def test_single_expression_directory_json_scalar(tmp_path):
    c=route_fixture(tmp_path,'WHB')
    s=b.resolve_spec(c,'WHB',{'expression_directories':'WHB-10Xv3'})
    assert s['expression_directories']==['WHB-10Xv3']

@pytest.mark.skipif(__import__('importlib').util.find_spec('anndata') is None,reason='AnnData absent here; independent live source comparison supplied but not executed.')
def test_independent_reference_comparison(tmp_path):
    path=h5ad(tmp_path/'x.h5ad');cm=pd.DataFrame({'cell_label':CELLS})
    m,_=b.subset_h5ad(path,cm,GENES)
    b.verify_with_anndata(path,cm,GENES,m,'cell_label')

def test_benchmark_profile_validation_and_no_fake_figures(tmp_path):
    import importlib.util
    root=Path(__file__).parents[2]/'scripts'
    for file,name in [('benchmark.py','bench'),('plot_benchmarks.py','plotter')]:
        spec=importlib.util.spec_from_file_location(name,root/file);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        if name=='bench':
            with pytest.raises(ValueError):module.validate_profile({'dataset':'WMB'})
        else:
            path=tmp_path/'fake.csv'
            pd.DataFrame([dict(dataset='WMB',cohort='neuronal',stage='load_data',cache_state='cold',replicate=1,
                 status='PASS',seconds=1,data_origin='SYNTHETIC',selection_sha256='fake',
                 return_type='sce',data_type='raw',manifest='20260711')]).to_csv(path,index=False)
            with pytest.raises(ValueError,match='Synthetic'):module.summarize(path)

def test_empty_optional_asset_category_is_not_an_error(tmp_path):
    class AllenLikeEmptyCache(FakeCache):
        def list_image_volume_files(self,d):
            raise ValueError(f'No image_volumes files found in directory {d}. image_volumes sub-directory is empty.')
    c=AllenLikeEmptyCache(tmp_path);c.add('D','cell_metadata',{'cell_label':['c']})
    assert b.files(c,'D','image_volumes') == []
    assert any(r['kind']=='metadata' for r in b.file_catalog(c))


def test_nonempty_manifest_error_is_not_swallowed(tmp_path):
    class BrokenCache(FakeCache):
        def list_image_volume_files(self,d):
            raise ValueError('manifest schema is corrupt')
    c=BrokenCache(tmp_path);c.add('D','cell_metadata',{'cell_label':['c']})
    with pytest.raises(ValueError,match='corrupt'):
        b.files(c,'D','image_volumes')


def test_expression_plan_never_queries_unrelated_directories(tmp_path):
    c=route_fixture(tmp_path,'MERFISH')
    c.add('ASAP-PMDBS-taxonomy','taxonomy_only',{'x':['y']})
    original=c.list_expression_matrix_files
    def guarded(directory):
        if directory == 'ASAP-PMDBS-taxonomy':
            raise AssertionError('unrelated directory was queried')
        return original(directory)
    c.list_expression_matrix_files=guarded
    cm,gd,p=b.load_data(c,'MERFISH')
    plans=b.expression_plan(c,p['spec'],cm)
    assert plans and {x['directory'] for x in plans} <= set(p['spec']['expression_directories'])

class _RemoteAttrs:
    def __init__(self, relative_path):
        self.relative_path=relative_path


class _FakeS3Client:
    def __init__(self, sizes=None, fail=False):
        self.sizes=sizes or {};self.fail=fail;self.head_calls=[];self.list_calls=[]
    def head_object(self,Bucket,Key):
        self.head_calls.append((Bucket,Key))
        if self.fail or Key not in self.sizes:
            raise RuntimeError('HEAD unavailable')
        return {'ContentLength':self.sizes[Key]}
    def list_objects_v2(self,Bucket,Prefix):
        self.list_calls.append((Bucket,Prefix))
        if self.fail or Prefix not in self.sizes:
            raise RuntimeError('LIST unavailable')
        return {'Contents':[{'Key':Prefix,'Size':self.sizes[Prefix]}]}


class _FakeInnerCache:
    bucket_name='synthetic-bucket'
    def __init__(self, records, sizes=None, fail_size=False):
        self.records=records
        self.s3_client=_FakeS3Client(sizes=sizes,fail=fail_size)
    def get_file_path(self,directory,file_name):
        return self.records[(directory,file_name)]
    def _check_successful_download(self,attrs):
        record=next(v for v in self.records.values() if v['file_attributes'] is attrs)
        return bool(record.get('verified',record.get('exists',False)))


def _make_remote(c,directory,file_name,size=None,cached=False,fail_size=False):
    key=f'expression_matrices/{directory}/20260711/{file_name.replace("/","-")}.h5ad'
    local=c.assets[(directory,'expression_matrices',file_name)]
    record={'local_path':local if cached else local.parent/'not-cached.h5ad',
            'exists':cached,'verified':cached,'file_attributes':_RemoteAttrs(key)}
    c._local=False
    c.cache=_FakeInnerCache({(directory,file_name):record},
                            sizes={} if size is None else {key:int(size)},
                            fail_size=fail_size)
    return key


def test_plan_uses_exact_file_size_not_directory_size(tmp_path):
    c=route_fixture(tmp_path,'PMDBS')
    cm,gd,prov=b.load_data(c,'PMDBS')
    key=_make_remote(c,'ASAP-PMDBS-10X','matrix/raw',size=123456,cached=False)
    def forbidden(*args):raise AssertionError('directory size must not be used')
    c.get_directory_expression_matrix_size=forbidden
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes']==123456
    assert result[0]['size_source']=='s3_head_object'
    assert result[0]['relative_path']==key
    assert result[0]['cached'] is False
    assert total==123456 and required==123456
    assert result[0]['file_name'].endswith('/raw')


def test_plan_cached_file_uses_exact_local_size(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    path=c.assets[('WHB-10Xv3','expression_matrices','matrix/raw')]
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=999999,cached=True)
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes']==path.stat().st_size
    assert result[0]['size_source']=='local_file'
    assert result[0]['cached'] is True and result[0]['download_required_bytes']==0
    assert total==path.stat().st_size and required==0


def test_plan_unknown_is_not_directory_size(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=None,cached=False,fail_size=True)
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes'] is None
    assert result[0]['size_source']=='unknown'
    assert total is None and required is None
    assert 'cannot obtain exact remote size' in result[0]['size_error']


def test_insufficient_disk_is_rejected_before_expression_download(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB');c.requested=[]
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=123456,cached=False)
    c._importabcatlas_download_base=tmp_path
    c._importabcatlas_disk_policy={'check_disk_space':True,'disk_reserve_bytes':10**30,
                                  'disk_safety_factor':1.05,'allow_unknown_size':False}
    with pytest.raises(OSError,match='Insufficient disk space') as error:
        b.fetch_data(c,prov['spec'],cm,gd,genes=[GENES[0]],allow_downloads=True)
    assert error.value.errno==__import__('errno').ENOSPC
    assert c.requested==[]


def test_unknown_remote_size_fails_safely_before_download(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB');c.requested=[]
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=None,cached=False,fail_size=True)
    c._importabcatlas_download_base=tmp_path
    c._importabcatlas_disk_policy={'check_disk_space':True,'disk_reserve_bytes':0,
                                  'disk_safety_factor':1.0,'allow_unknown_size':False}
    with pytest.raises(RuntimeError,match='exact size is unknown'):
        b.fetch_data(c,prov['spec'],cm,gd,genes=[GENES[0]],allow_downloads=True)
    assert c.requested==[]


def test_verified_cached_expression_needs_no_download_consent(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB');c.requested=[]
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=999999,cached=True)
    m,_,_,_=b.fetch_data(c,prov['spec'],cm,gd,genes=[GENES[0]],allow_downloads=False)
    assert m.shape==(1,len(cm))


def test_plan_log2_only_when_requested(tmp_path):
    c=route_fixture(tmp_path,'PMDBS');cm,gd,prov=b.load_data(c,'PMDBS')
    c.add('ASAP-PMDBS-10X','matrix/log2',h5ad(tmp_path/'log2.h5ad'),kind='expression_matrices')
    result,_,_=b.make_fetch_plan(c,prov['spec'],cm,gd,'log2')
    assert len(result)==1 and result[0]['file_name'].endswith('/log2')


def test_complex_multi_file_fetch_preserves_requested_axes(tmp_path):
    c=FakeCache(tmp_path)
    cells=['cell-A1','cell-B1','cell-A2','cell-B2']
    genes=['gene-1','gene-2','gene-3','gene-4','gene-5']
    symbols=['A','B','C','D','E']
    cell=pd.DataFrame({'cell_label':cells,
                       'feature_matrix_label':['group-A','group-B','group-A','group-B']})
    gene=pd.DataFrame({'gene_identifier':genes,'gene_symbol':symbols})
    # Each source contains only its own cells, as in a split 10X collection.
    values_a=np.array([[1,2,3,4,5],[11,12,13,14,15]],dtype=float)
    values_b=np.array([[21,22,23,24,25],[31,32,33,34,35]],dtype=float)
    pa=h5ad(tmp_path/'a.h5ad',values=values_a,cells=['cell-A1','cell-A2'],genes=genes)
    pb=h5ad(tmp_path/'b.h5ad',values=values_b,cells=['cell-B1','cell-B2'],genes=genes)
    c.add('D1','group-A/raw',pa,'expression_matrices')
    c.add('D2','group-B/raw',pb,'expression_matrices')
    spec={'directory':'D1','expression_directories':['D1','D2'],'imputed':False}
    requested=cell.set_index('cell_label',drop=False).loc[['cell-B2','cell-A1','cell-B1']].reset_index(drop=True)
    matrix,out_cell,out_gene,prov=b.fetch_data(
        c,spec,requested,gene,genes=['E','gene-2','A'],data_type='raw')
    expected=np.array([[35,5,25],[32,2,22],[31,1,21]],dtype=float)
    np.testing.assert_array_equal(matrix.toarray(),expected)
    assert out_cell.cell_label.tolist()==['cell-B2','cell-A1','cell-B1']
    assert out_gene.gene_identifier.tolist()==['gene-5','gene-2','gene-1']
    assert len(prov['sources'])==2


def test_log2_fetch_is_not_labeled_counts(tmp_path):
    c=FakeCache(tmp_path)
    cell=pd.DataFrame({'cell_label':CELLS,'feature_matrix_label':['matrix']*3})
    gene=pd.DataFrame({'gene_identifier':GENES,'gene_symbol':['Gad1','Aqp4','X','Y']})
    path=h5ad(tmp_path/'log2.h5ad',values=VALUES/3)
    c.add('D','matrix/log2',path,'expression_matrices')
    spec={'directory':'D','expression_directories':['D'],'imputed':False}
    matrix,_,_,prov=b.fetch_data(c,spec,cell,gene,genes=[GENES[0],GENES[1]],data_type='log2')
    assert matrix.shape==(2,3)
    assert prov['assay_kind']=='log2'
    assert prov['integer_nonnegative'] is False


def test_fetch_validation_errors_happen_before_expression_file_access(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB');c.requested=[]
    with pytest.raises(ValueError,match='data_type'):
        b.fetch_data(c,prov['spec'],cm,gd,data_type='scaled')
    assert c.requested==[]
    with pytest.raises(ValueError,match='Gene'):
        b.fetch_data(c,prov['spec'],cm,gd,genes=['definitely-not-a-gene'])
    assert c.requested==[]
    with pytest.raises(MemoryError,match='max_output_mb'):
        b.fetch_data(c,prov['spec'],cm,gd,genes=GENES,max_output_mb=1e-20)
    assert c.requested==[]


def test_duplicate_cell_query_rejected_before_plan(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    duplicate=pd.concat([cm.iloc[[0]],cm.iloc[[0]]],ignore_index=True)
    c.requested=[]
    with pytest.raises(ValueError,match='duplicate'):
        b.fetch_data(c,prov['spec'],duplicate,gd,genes=[GENES[0]])
    assert c.requested==[]


def test_s3_list_fallback_reports_exact_file_size(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    key=_make_remote(c,'WHB-10Xv3','matrix/raw',size=7654321,cached=False)
    client=c.cache.s3_client
    def head_fails(**kwargs):
        raise RuntimeError('HEAD deliberately unavailable')
    client.head_object=head_fails
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['relative_path']==key
    assert result[0]['file_bytes']==7654321
    assert result[0]['size_source']=='s3_list_objects_v2'
    assert total==required==7654321


def test_unknown_size_can_be_deliberately_allowed(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB');c.requested=[]
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=None,cached=False,fail_size=True)
    c._importabcatlas_download_base=tmp_path
    c._importabcatlas_disk_policy={'check_disk_space':True,'disk_reserve_bytes':0,
                                  'disk_safety_factor':1.0,'allow_unknown_size':True}
    m,_,_,_=b.fetch_data(c,prov['spec'],cm,gd,genes=[GENES[0]],allow_downloads=True)
    assert m.shape==(1,len(cm))
    assert ('WHB-10Xv3','matrix/raw') in c.requested


def test_metadata_download_disk_preflight_uses_exact_file_size(tmp_path):
    c=FakeCache(tmp_path)
    path=c.add('D','cell_metadata',{'cell_label':['a']})
    c._local=False
    key='metadata/D/20260711/cell_metadata.csv'
    record={'local_path':path.parent/'not-cached.csv','exists':False,
            'file_attributes':_RemoteAttrs(key)}
    c.cache=_FakeInnerCache({('D','cell_metadata'):record},sizes={key:100})
    c._importabcatlas_download_base=tmp_path
    c._importabcatlas_disk_policy={'check_disk_space':True,'disk_reserve_bytes':10**30,
                                  'disk_safety_factor':1.0,'allow_unknown_size':False}
    c.requested=[]
    with pytest.raises(OSError,match='Insufficient disk space'):
        b.get_path(c,'D','cell_metadata')
    assert c.requested==[]


def test_plan_deduplicates_same_source_file_for_totals(tmp_path):
    c=route_fixture(tmp_path,'WHB')
    cm,gd,prov=b.load_data(c,'WHB')
    # Two labels are explicitly mapped to the same physical file.
    cm=cm.copy();cm['feature_matrix_label']=['a','b','a']
    fmap={label:{'directory':'WHB-10Xv3','file_name':'matrix/raw'} for label in ['a','b']}
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=456789,cached=False)
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw',fmap)
    assert len(result)==2
    assert total==456789 and required==456789


def test_corrupt_h5ad_is_rejected_cleanly(tmp_path):
    path=tmp_path/'bad.h5ad'
    with h5py.File(path,'w') as handle:
        handle.create_group('obs')
    with pytest.raises(ValueError,match='lacks'):
        b.subset_h5ad(path,pd.DataFrame({'cell_label':['a']}),['g'])


def test_unverified_partial_local_file_does_not_define_remote_size(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    directory,file_name='WHB-10Xv3','matrix/raw'
    key='expression_matrices/WHB-10Xv3/20260711/matrix-raw.h5ad'
    partial=tmp_path/'partial.h5ad';partial.write_bytes(b'partial')
    record={'local_path':partial,'exists':True,'verified':False,
            'file_attributes':_RemoteAttrs(key)}
    c._local=False
    c.cache=_FakeInnerCache({(directory,file_name):record},sizes={key:9876543})
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['cached'] is False
    assert result[0]['file_bytes']==9876543
    assert result[0]['size_source']=='s3_head_object'
    assert total==9876543
    assert required==9876543-len(b'partial')
    assert result[0]['partial_bytes']==len(b'partial')


def test_offline_missing_expression_reports_missing_cache_not_consent(tmp_path):
    c=route_fixture(tmp_path,'WHB')
    cm,gd,prov=b.load_data(c,'WHB')
    directory,file_name='WHB-10Xv3','matrix/raw'
    # Point the manifest record at a file that is not present in an offline cache.
    missing=tmp_path/'offline-missing.h5ad'
    class LocalRecordCache:
        def get_file_path(self,directory,file_name):
            return {'local_path':missing,'exists':False,
                    'file_attributes':_RemoteAttrs('expression_matrices/missing.h5ad')}
    c._local=True
    c.cache=LocalRecordCache()
    with pytest.raises(FileNotFoundError,match='Offline cache lacks required expression files'):
        b.fetch_data(c,prov['spec'],cm,gd,genes=[GENES[0]],allow_downloads=False)


def test_offline_missing_metadata_reports_file_unavailable(tmp_path):
    c=FakeCache(tmp_path)
    c._local=True
    class LocalRecordCache:
        def get_file_path(self,directory,file_name):
            return {'local_path':tmp_path/'missing.csv','exists':False,
                    'file_attributes':_RemoteAttrs('metadata/D/missing.csv')}
    c.cache=LocalRecordCache()
    def local_get_file_path(directory,file_name,**kwargs):
        return tmp_path/'missing.csv'
    c.get_file_path=local_get_file_path
    with pytest.raises(FileNotFoundError,match='File unavailable in this cache'):
        b.get_path(c,'D','missing')


def test_plan_disk_status_reports_insufficient_reserve(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    _make_remote(c,'WHB-10Xv3','matrix/raw',size=123,cached=False)
    c._importabcatlas_download_base=tmp_path
    c._importabcatlas_disk_policy={'check_disk_space':True,'disk_reserve_bytes':10**30,
                                  'disk_safety_factor':1.05,'allow_unknown_size':False}
    _,_,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    status=b._disk_status(c,required)
    assert status['required_free_bytes']>status['free_bytes']
    assert status['enough_space'] is False


def test_stream_filter_preview_stops_at_requested_matches(tmp_path):
    path=tmp_path/'many.csv'
    rows=pd.DataFrame({
        'cell_label':[f'other-{i}' for i in range(100005)]+['target-1','target-2','target-3'],
        'donor_label':['other']*100005+['human']*3,
    })
    rows.to_csv(path,index=False)
    out=b.read_csv(path,wanted={'donor_label':['human']},preview_n=2)
    assert out.cell_label.tolist()==['target-1','target-2']


def test_hmba_species_exact_cell_selection_uses_species_and_cell_keys(tmp_path):
    c=route_fixture(tmp_path,'HMBA_Human')
    # Add one non-human donor/cell before the human cells in the shared table.
    cell_path=c.assets[('HMBA-10XMultiome-BG','metadata','cell_metadata')]
    cells=pd.read_csv(cell_path,dtype=str)
    extra=cells.iloc[[0]].copy();extra['cell_label']='non-human-cell';extra['donor_label']='non-human-donor'
    pd.concat([extra,cells],ignore_index=True).to_csv(cell_path,index=False)
    donor_path=c.assets[('HMBA-10XMultiome-BG','metadata','donor')]
    donors=pd.read_csv(donor_path,dtype=str)
    donors=pd.concat([pd.DataFrame({'donor_label':['non-human-donor'],'diagnosis':['Unknown'],
                                   'species_genus':['Macaque']}),donors],ignore_index=True)
    donors.to_csv(donor_path,index=False)
    requested=CELLS[1]
    cm,_,_=b.load_data(c,'HMBA_Human',cell_ids=[requested])
    assert cm.cell_label.tolist()==[requested]
    assert 'non-human-cell' not in cm.cell_label.tolist()


def test_manifest_attribute_file_size_is_exact_and_no_s3_request(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    directory,file_name='WHB-10Xv3','matrix/raw'
    class SizedAttrs:
        relative_path='expression_matrices/WHB-10Xv3/20260711/matrix.h5ad'
        file_size=24680
    record={'local_path':tmp_path/'not-cached.h5ad','exists':False,
            'file_attributes':SizedAttrs()}
    c._local=False
    c.cache=_FakeInnerCache({(directory,file_name):record},sizes={},fail_size=True)
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes']==24680
    assert result[0]['size_source']=='manifest_attributes'
    assert total==required==24680


def test_temporary_output_disk_check_raises_enospc(tmp_path,monkeypatch):
    c=route_fixture(tmp_path,'WHB')
    cm,gd,prov=b.load_data(c,'WHB')
    cells=tmp_path/'cells.csv';genes=tmp_path/'genes.csv'
    cm.to_csv(cells,index=False);gd.to_csv(genes,index=False)
    monkeypatch.setattr(b,'create_cache',lambda *args,**kwargs:c)
    matrix=sparse.csc_matrix(np.ones((1,1),dtype=float))
    out_cell=cm.iloc[[0]].copy();out_gene=gd.iloc[[0]].copy()
    monkeypatch.setattr(b,'fetch_data',lambda *args,**kwargs:
                        (matrix,out_cell,out_gene,{'assay_kind':'counts'}))
    usage=__import__('collections').namedtuple('usage','total used free')
    monkeypatch.setattr(b.shutil,'disk_usage',lambda path:usage(100,99,1))
    request={'action':'fetch','out_dir':str(tmp_path/'out'),'args':{
        'download_base':str(tmp_path/'cache'),'release':'20260711','offline':False,
        'dataset':'WHB','overrides':{},'metadata_path':str(cells),'gene_path':str(genes),
        'filters':{},'cell_ids':[out_cell.cell_label.iloc[0]],'n_cells':None,'seed':1,
        'genes':[out_gene.gene_identifier.iloc[0]],'data_type':'raw','allow_downloads':True,
        'max_output_mb':1024,'file_map':None,'obs_column':None,'counts_semantics':None,
        'verify_reference':False}}
    with pytest.raises(OSError,match='temporary-disk space') as error:
        b.run_request(request)
    assert error.value.errno==__import__('errno').ENOSPC


def test_empty_optional_asset_runtimeerror_is_not_an_error(tmp_path):
    class AllenLikeEmptyCache(FakeCache):
        def list_mapmycells_files(self,d):
            raise RuntimeError(f'No mapmycells files found in directory {d}. mapmycells sub-directory is empty.')
    c=AllenLikeEmptyCache(tmp_path);c.add('D','cell_metadata',{'cell_label':['c']})
    assert b.files(c,'D','mapmycells') == []


def test_patchseq_uses_declared_log2_default_and_rejects_explicit_raw(tmp_path):
    c=route_fixture(tmp_path,'HMBA_PatchSeq')
    cm,gd,prov=b.load_data(c,'HMBA_PatchSeq')
    plans,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,data_type=None)
    assert plans[0]['file_name'].endswith('/log2')
    assert b.resolve_data_type(prov['spec'],None)=='log2'
    with pytest.raises(ValueError,match='No raw expression file|Available data_type'):
        b.make_fetch_plan(c,prov['spec'],cm,gd,data_type='raw')


def test_duplicate_gene_rows_are_disambiguated_without_aggregation(tmp_path):
    source_ids=['NCBIGene:1','NCBIGene:1','NCBIGene:2']
    unique_ids=['NCBIGene:1.iabdup1','NCBIGene:1.iabdup2','NCBIGene:2']
    gene_raw=pd.DataFrame({'gene_identifier':source_ids,
                           'gene_symbol':['A-first','A-second','B']})
    gene,provenance=b.normalize_gene_metadata(gene_raw,'HMBA-Macaque')
    assert gene.gene_identifier.tolist()==unique_ids
    assert gene.source_gene_identifier.tolist()==source_ids
    assert provenance['applied'] is True
    with pytest.raises(ValueError,match='2 identifiers/rows') as error:
        b.resolve_genes(gene,['NCBIGene:1'])
    assert 'NCBIGene:1.iabdup1' in str(error.value)
    values=np.array([[1,2,3],[4,5,6]],dtype=float)
    path=h5ad(tmp_path/'dup.h5ad',values=values,cells=['c1','c2'],genes=source_ids)
    matrix,_=b.subset_h5ad(path,pd.DataFrame({'cell_label':['c2','c1']}),
                           ['NCBIGene:1.iabdup2','NCBIGene:1.iabdup1'])
    np.testing.assert_array_equal(matrix.toarray(),np.array([[5,4],[2,1]],dtype=float))


def test_transient_exact_size_lookup_is_retried(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    key=_make_remote(c,'WHB-10Xv3','matrix/raw',size=333,cached=False)
    client=c.cache.s3_client
    original=client.head_object
    calls={'n':0}
    def flaky(**kwargs):
        calls['n']+=1
        if calls['n']==1:
            raise TimeoutError('temporary read timeout')
        return original(**kwargs)
    client.head_object=flaky
    c._importabcatlas_network_policy={'network_retries':1,'network_backoff_seconds':0}
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['relative_path']==key
    assert result[0]['file_bytes']==333
    assert calls['n']==2 and total==required==333


def test_retry_call_does_not_hide_nontransient_error():
    calls={'n':0}
    def broken():
        calls['n']+=1
        raise ValueError('schema is invalid')
    with pytest.raises(ValueError,match='schema'):
        b._retry_call(broken,{'network_retries':5,'network_backoff_seconds':0},'test')
    assert calls['n']==1


def test_resilient_s3_download_resumes_partial_bytes(tmp_path,monkeypatch):
    import hashlib
    import types

    class BaseCache:
        def __init__(self,cache_dir,bucket_name,ui_class_name=None,auth_required=False):
            self._cache_dir=Path(cache_dir)
            self._bucket_name=bucket_name
            self._auth_required=auth_required

    def file_hash(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    package=types.ModuleType('abc_atlas_access')
    cache_pkg=types.ModuleType('abc_atlas_access.abc_atlas_cache')
    cloud=types.ModuleType('abc_atlas_access.abc_atlas_cache.cloud_cache')
    utils=types.ModuleType('abc_atlas_access.abc_atlas_cache.utils')
    cloud.S3CloudCache=BaseCache
    utils.file_hash_from_path=file_hash
    monkeypatch.setitem(sys.modules,'abc_atlas_access',package)
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache',cache_pkg)
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.cloud_cache',cloud)
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.utils',utils)

    payload=b'abcdefghij0123456789'
    class Body:
        def __init__(self,chunks,fail=False):self.chunks=chunks;self.fail=fail;self.closed=False
        def iter_chunks(self,chunk_size):
            for chunk in self.chunks:yield chunk
            if self.fail:raise BrokenPipeError('interrupted')
        def close(self):self.closed=True
    class Client:
        def __init__(self):self.requests=[]
        def head_object(self,Bucket,Key):return {'ContentLength':len(payload)}
        def get_object(self,**request):
            self.requests.append(dict(request))
            if len(self.requests)==1:
                return {'Body':Body([payload[:7]],fail=True),
                        'ResponseMetadata':{'HTTPStatusCode':200}}
            assert request['Range']=='bytes=7-'
            return {'Body':Body([payload[7:]]),
                    'ContentRange':f'bytes 7-{len(payload)-1}/{len(payload)}',
                    'ResponseMetadata':{'HTTPStatusCode':206}}

    policy=b._network_policy(retries=1,backoff_seconds=0,chunk_mb=0.001)
    cls=b._resilient_s3_cache_class('releases/20260711/manifest.json',policy)
    obj=cls(tmp_path,'bucket')
    client=Client();obj._s3_client=client
    destination=tmp_path/'asset.bin'
    expected=hashlib.sha256(payload).hexdigest()
    result=obj._stream_object('path/asset.bin',destination,expected_hash=expected)
    assert result==destination and destination.read_bytes()==payload
    assert not Path(str(destination)+'.part').exists()
    assert client.requests[1]['Range']=='bytes=7-'
    assert obj._list_all_manifests()==['releases/20260711/manifest.json']


def test_completed_partial_is_finalized_without_new_get_request(tmp_path,monkeypatch):
    import hashlib
    import types
    class BaseCache:
        def __init__(self,cache_dir,bucket_name,ui_class_name=None,auth_required=False):
            self._cache_dir=Path(cache_dir);self._bucket_name=bucket_name;self._auth_required=auth_required
    cloud=types.ModuleType('abc_atlas_access.abc_atlas_cache.cloud_cache');cloud.S3CloudCache=BaseCache
    utils=types.ModuleType('abc_atlas_access.abc_atlas_cache.utils')
    utils.file_hash_from_path=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    monkeypatch.setitem(sys.modules,'abc_atlas_access',types.ModuleType('abc_atlas_access'))
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache',types.ModuleType('abc_atlas_access.abc_atlas_cache'))
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.cloud_cache',cloud)
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.utils',utils)
    payload=b'complete-partial'
    class Client:
        def head_object(self,**kwargs):return {'ContentLength':len(payload)}
        def get_object(self,**kwargs):raise AssertionError('no transfer should occur')
    cls=b._resilient_s3_cache_class('releases/20260711/manifest.json',b._network_policy(retries=0))
    obj=cls(tmp_path,'bucket');obj._s3_client=Client()
    destination=tmp_path/'asset.bin';Path(str(destination)+'.part').write_bytes(payload)
    result=obj._stream_object('asset.bin',destination,hashlib.sha256(payload).hexdigest())
    assert result.read_bytes()==payload


def test_structured_s3_server_error_is_retryable_but_permanent_4xx_is_not():
    class StructuredError(Exception):
        def __init__(self,status,code):
            super().__init__(code)
            self.response={'ResponseMetadata':{'HTTPStatusCode':status},
                           'Error':{'Code':code}}
    assert b._is_transient_network_error(StructuredError(503,'ServiceUnavailable'))
    assert b._is_transient_network_error(StructuredError(429,'Throttling'))
    assert not b._is_transient_network_error(StructuredError(404,'NoSuchKey'))
    assert not b._is_transient_network_error(StructuredError(403,'AccessDenied'))


def test_wrong_range_response_is_closed_and_restarted_from_zero(tmp_path,monkeypatch):
    import hashlib
    import types

    class BaseCache:
        def __init__(self,cache_dir,bucket_name,ui_class_name=None,auth_required=False):
            self._cache_dir=Path(cache_dir)
            self._bucket_name=bucket_name
            self._auth_required=auth_required

    cloud=types.ModuleType('abc_atlas_access.abc_atlas_cache.cloud_cache')
    cloud.S3CloudCache=BaseCache
    utils=types.ModuleType('abc_atlas_access.abc_atlas_cache.utils')
    utils.file_hash_from_path=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    monkeypatch.setitem(sys.modules,'abc_atlas_access',types.ModuleType('abc_atlas_access'))
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache',types.ModuleType('abc_atlas_access.abc_atlas_cache'))
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.cloud_cache',cloud)
    monkeypatch.setitem(sys.modules,'abc_atlas_access.abc_atlas_cache.utils',utils)

    payload=b'correct-complete-payload'
    class Body:
        def __init__(self,data):self.data=data;self.closed=False;self.iterated=False
        def iter_chunks(self,chunk_size):
            self.iterated=True
            yield self.data
        def close(self):self.closed=True
    class Client:
        def __init__(self):self.requests=[];self.wrong_body=None
        def head_object(self,**kwargs):return {'ContentLength':len(payload)}
        def get_object(self,**kwargs):
            self.requests.append(dict(kwargs))
            if len(self.requests)==1:
                self.wrong_body=Body(payload)
                return {'Body':self.wrong_body,
                        'ContentRange':'bytes 0-4/23',
                        'ResponseMetadata':{'HTTPStatusCode':206}}
            assert 'Range' not in kwargs
            return {'Body':Body(payload),'ResponseMetadata':{'HTTPStatusCode':200}}

    cls=b._resilient_s3_cache_class('releases/20260711/manifest.json',
                                    b._network_policy(retries=0,chunk_mb=0.001))
    obj=cls(tmp_path,'bucket');client=Client();obj._s3_client=client
    destination=tmp_path/'asset.bin'
    Path(str(destination)+'.part').write_bytes(payload[:5])
    result=obj._stream_object('asset.bin',destination,hashlib.sha256(payload).hexdigest())
    assert result==destination and destination.read_bytes()==payload
    assert client.requests[0]['Range']=='bytes=5-'
    assert 'Range' not in client.requests[1]
    assert client.wrong_body.closed is True
    assert client.wrong_body.iterated is False


def test_registry_declares_accessible_modalities_and_representations():
    assert len(b.REGISTRY) == 23
    for dataset, spec in b.REGISTRY.items():
        available = spec.get('available_data_types')
        assert available
        assert set(available).issubset({'raw', 'log2'})
        assert spec.get('default_data_type') in available
        assert spec.get('modality') in {
            'single_cell_transcriptomics',
            'single_nucleus_transcriptomics',
            'spatial_transcriptomics',
            'imputed_spatial_transcriptomics',
            'patch_seq_transcriptomics',
        }
    assert b.REGISTRY['MERFISH_imputed']['available_data_types'] == ['log2']
    assert b.REGISTRY['HMBA_PatchSeq']['available_data_types'] == ['log2']
    assert b.REGISTRY['MERFISH']['available_data_types'] == ['raw', 'log2']


def test_unpublished_representation_rejected_before_inventory(tmp_path):
    spec = dict(b.REGISTRY['HMBA_PatchSeq'])
    with pytest.raises(ValueError, match='does not publish.*raw'):
        b.resolve_data_type(spec, 'raw')
    with pytest.raises(ValueError, match='must be one of'):
        b.resolve_data_type(spec, 'invalid')


def _add_var_column(path, name, values):
    text = h5py.string_dtype('utf-8')
    with h5py.File(path, 'r+') as handle:
        if name in handle['var']:
            del handle['var'][name]
        handle['var'].create_dataset(name, data=np.asarray(values, dtype=text))


def test_feature_matching_uses_gene_symbol_when_var_names_are_symbols(tmp_path):
    genes = ['SYMA', 'SYMB']
    path = h5ad(
        tmp_path / 'symbol-index.h5ad',
        values=np.array([[1, 2], [3, 4]], dtype=float),
        cells=['c1', 'c2'],
        genes=genes,
    )
    metadata = pd.DataFrame({
        'gene_identifier': ['NCBIGene:1', 'NCBIGene:2'],
        'gene_symbol': ['SYMA', 'SYMB'],
    })
    matrix, _, feature = b.subset_h5ad(
        path,
        pd.DataFrame({'cell_label': ['c2', 'c1']}),
        metadata.iloc[[1, 0]].reset_index(drop=True),
        return_feature_matches=True,
    )
    np.testing.assert_array_equal(matrix.toarray(), np.array([[4, 3], [2, 1]]))
    assert feature['methods'] == ['var_names/gene_symbol', 'var_names/gene_symbol']


def test_feature_matching_uses_h5ad_gene_identifier_column(tmp_path):
    path = h5ad(
        tmp_path / 'var-column.h5ad',
        values=np.array([[1, 2], [3, 4]], dtype=float),
        cells=['c1', 'c2'],
        genes=['feature-1', 'feature-2'],
    )
    _add_var_column(path, 'gene_identifier', ['NCBIGene:1', 'NCBIGene:2'])
    metadata = pd.DataFrame({
        'gene_identifier': ['NCBIGene:1', 'NCBIGene:2'],
        'gene_symbol': ['A', 'B'],
    })
    matrix, _, feature = b.subset_h5ad(
        path,
        pd.DataFrame({'cell_label': ['c1']}),
        metadata,
        return_feature_matches=True,
    )
    np.testing.assert_array_equal(matrix.toarray(), np.array([[1, 2]]))
    assert feature['methods'] == ['var/gene_identifier', 'var/gene_identifier']


def test_all_genes_uses_authoritative_h5ad_intersection_for_nonimputed_route(tmp_path):
    c = route_fixture(tmp_path, 'HMBA_Macaque')
    gene_path = c.assets[('HMBA-10XMultiome-BG', 'metadata', 'macaque_gene')]
    gene = pd.read_csv(gene_path, dtype=str)
    extra = pd.DataFrame({'gene_identifier': ['not-measured'], 'gene_symbol': ['NOPE']})
    pd.concat([gene, extra], ignore_index=True).to_csv(gene_path, index=False)
    cm, gd, prov = b.load_data(c, 'HMBA_Macaque')
    matrix, _, selected, provenance = b.fetch_data(
        c, prov['spec'], cm, gd, genes=None, data_type='raw'
    )
    assert matrix.shape == (len(GENES), len(CELLS))
    assert selected.gene_identifier.tolist() == GENES
    assert 'not-measured' not in selected.gene_identifier.tolist()
    reconciliation = provenance['gene_axis_reconciliation']
    assert reconciliation['mode'] == 'authoritative_h5ad_intersection'
    assert reconciliation['dropped_unmeasured_rows'] == 1


def test_explicit_unmeasured_gene_remains_strict_error(tmp_path):
    c = route_fixture(tmp_path, 'MERFISH_imputed')
    cm, gd, prov = b.load_data(c, 'MERFISH_imputed')
    gd = pd.concat([
        gd,
        pd.DataFrame({'gene_identifier': ['not-imputed'], 'gene_symbol': ['NOPE']})
    ], ignore_index=True)
    with pytest.raises(ValueError, match='unmeasured/absent.*genes=NULL'):
        b.fetch_data(
            c, prov['spec'], cm, gd,
            genes=['not-imputed'], data_type='log2'
        )


def test_all_gene_intersection_across_multiple_source_files(tmp_path):
    c = FakeCache(tmp_path)
    cells = pd.DataFrame({
        'cell_label': ['a', 'b'],
        'feature_matrix_label': ['one', 'two'],
    })
    gene = pd.DataFrame({
        'gene_identifier': ['g1', 'g2', 'g3'],
        'gene_symbol': ['G1', 'G2', 'G3'],
    })
    p1 = h5ad(
        tmp_path / 'one.h5ad',
        values=np.array([[1, 2, 3]], dtype=float),
        cells=['a'], genes=['g1', 'g2', 'g3']
    )
    p2 = h5ad(
        tmp_path / 'two.h5ad',
        values=np.array([[4, 5]], dtype=float),
        cells=['b'], genes=['g1', 'g3']
    )
    c.add('D', 'one/raw', p1, 'expression_matrices')
    c.add('D', 'two/raw', p2, 'expression_matrices')
    spec = {
        'directory': 'D',
        'expression_directories': ['D'],
        'default_data_type': 'raw',
        'available_data_types': ['raw'],
    }
    matrix, _, selected, provenance = b.fetch_data(
        c, spec, cells, gene, genes=None, data_type='raw'
    )
    assert selected.gene_identifier.tolist() == ['g1', 'g3']
    np.testing.assert_array_equal(matrix.toarray(), np.array([[1, 4], [3, 5]]))
    assert provenance['gene_axis_reconciliation']['dropped_unmeasured_rows'] == 1


def test_duplicate_stable_ids_align_by_occurrence_when_h5ad_index_differs(tmp_path):
    """Native HMBA Macaque duplicate NCBI IDs remain one-to-one, not aggregated."""
    values = np.array([[11, 12, 13], [21, 22, 23]], dtype=float)
    path = h5ad(
        tmp_path / 'macaque-duplicate-axis.h5ad',
        encoding='dense',
        values=values,
        cells=['c1', 'c2'],
        genes=['feature-a', 'feature-b', 'feature-c'],
    )
    _add_var_column(path, 'gene_identifier', ['NCBIGene:1', 'NCBIGene:1', 'NCBIGene:2'])
    _add_var_column(path, 'gene_symbol', ['DUP', 'DUP', 'UNIQUE'])

    raw = pd.DataFrame({
        'gene_identifier': ['NCBIGene:1', 'NCBIGene:1', 'NCBIGene:2'],
        'gene_symbol': ['DUP', 'DUP', 'UNIQUE'],
    })
    metadata, provenance = b.normalize_gene_metadata(raw, 'HMBA-Macaque')
    assert provenance['occurrence_alignment'] is True
    assert metadata.source_gene_occurrence.tolist() == [1, 2, 1]
    assert metadata.source_gene_group_size.tolist() == [2, 2, 1]

    requested = metadata.iloc[[1, 0, 2]].reset_index(drop=True)
    matrix, _, feature = b.subset_h5ad(
        path,
        pd.DataFrame({'cell_label': ['c2', 'c1']}),
        requested,
        return_feature_matches=True,
    )
    np.testing.assert_array_equal(
        matrix.toarray(),
        np.array([[22, 21, 23], [12, 11, 13]], dtype=float),
    )
    assert feature['methods'][:2] == [
        'var/gene_identifier<-source_gene_identifier/occurrence',
        'var/gene_identifier<-source_gene_identifier/occurrence',
    ]


def test_authoritative_axis_collapses_metadata_aliases_before_dense_indexing(tmp_path):
    """Imputed-style metadata supersets yield one row per measured feature."""
    values = np.array([[1, 2], [3, 4]], dtype=float)
    path = h5ad(
        tmp_path / 'imputed-authoritative-axis.h5ad',
        encoding='dense',
        values=values,
        cells=['c1', 'c2'],
        genes=['ENSG1', 'ENSG2'],
    )
    _add_var_column(path, 'gene_symbol', ['A', 'B'])
    metadata = pd.DataFrame({
        'gene_identifier': ['ENSG1', 'ALIAS-OF-A', 'ENSG2'],
        'gene_symbol': ['A', 'A', 'B'],
    })

    selected, reconciliation = b._measured_gene_intersection([path], metadata)
    assert selected.gene_identifier.tolist() == ['ENSG1', 'ENSG2']
    assert reconciliation['selected_measured_rows'] == 2
    assert reconciliation['metadata_alias_rows_collapsed'] == 1

    matrix, _, feature = b.subset_h5ad(
        path,
        pd.DataFrame({'cell_label': ['c2', 'c1']}),
        selected.iloc[[1, 0]].reset_index(drop=True),
        return_feature_matches=True,
    )
    np.testing.assert_array_equal(
        matrix.toarray(),
        np.array([[4, 3], [2, 1]], dtype=float),
    )
    assert feature['methods'] == ['var_names/unique', 'var_names/unique']


def test_explicit_aliases_to_same_feature_fail_before_h5py_indexing(tmp_path):
    path = h5ad(
        tmp_path / 'duplicate-feature-request.h5ad',
        encoding='dense',
        values=np.array([[1, 2]], dtype=float),
        cells=['c1'],
        genes=['ENSG1', 'ENSG2'],
    )
    _add_var_column(path, 'gene_symbol', ['A', 'B'])
    metadata = pd.DataFrame({
        'gene_identifier': ['ENSG1', 'ALIAS-OF-A'],
        'gene_symbol': ['A', 'A'],
    })
    with pytest.raises(ValueError, match='one-to-one.*genes=NULL'):
        b.subset_h5ad(
            path,
            pd.DataFrame({'cell_label': ['c1']}),
            metadata,
        )


def _expression_asset_path(cache):
    matches = [
        path for (directory, kind, name), path in cache.assets.items()
        if kind == 'expression_matrices'
    ]
    assert len(matches) == 1
    return matches[0]


def test_fetch_data_hmba_macaque_duplicate_source_ids_end_to_end(tmp_path):
    c = route_fixture(tmp_path, 'HMBA_Macaque', encoding='dense')
    gene_path = c.assets[('HMBA-10XMultiome-BG', 'metadata', 'macaque_gene')]
    source_ids = ['NCBIGene:1', 'NCBIGene:1', 'NCBIGene:2', 'NCBIGene:3']
    symbols = ['DUP', 'DUP', 'B', 'C']
    pd.DataFrame({
        'gene_identifier': source_ids,
        'gene_symbol': symbols,
    }).to_csv(gene_path, index=False)
    expression_path = _expression_asset_path(c)
    with h5py.File(expression_path, 'r+') as handle:
        del handle['var']['_index']
        text = h5py.string_dtype('utf-8')
        handle['var'].create_dataset(
            '_index', data=np.asarray(['f1', 'f2', 'f3', 'f4'], dtype=text)
        )
    _add_var_column(expression_path, 'gene_identifier', source_ids)
    _add_var_column(expression_path, 'gene_symbol', symbols)

    cells, genes, provenance = b.load_data(c, 'HMBA_Macaque')
    matrix, returned_cells, selected, fetch_provenance = b.fetch_data(
        c, provenance['spec'], cells, genes, genes=None, data_type='raw'
    )
    assert selected.gene_identifier.tolist() == [
        'NCBIGene:1.iabdup1', 'NCBIGene:1.iabdup2', 'NCBIGene:2', 'NCBIGene:3'
    ]
    assert matrix.shape == (4, 3)
    np.testing.assert_array_equal(matrix.toarray(), VALUES.T)
    assert returned_cells.cell_label.tolist() == CELLS
    reconciliation = fetch_provenance['gene_axis_reconciliation']
    assert reconciliation['selected_measured_rows'] == 4


def test_fetch_data_imputed_dense_collapses_metadata_alias_rows(tmp_path):
    c = route_fixture(tmp_path, 'MERFISH_imputed', encoding='dense')
    gene_path = c.assets[('WMB-10X', 'metadata', 'gene')]
    pd.DataFrame({
        'gene_identifier': [GENES[0], 'alias-of-first', *GENES[1:]],
        'gene_symbol': ['A', 'A', 'B', 'C', 'D'],
    }).to_csv(gene_path, index=False)
    expression_path = _expression_asset_path(c)
    _add_var_column(expression_path, 'gene_symbol', ['A', 'B', 'C', 'D'])

    cells, genes, provenance = b.load_data(c, 'MERFISH_imputed')
    matrix, returned_cells, selected, fetch_provenance = b.fetch_data(
        c, provenance['spec'], cells, genes, genes=None, data_type='log2'
    )
    assert selected.gene_identifier.tolist() == GENES
    assert matrix.shape == (4, 3)
    np.testing.assert_array_equal(matrix.toarray(), VALUES.T)
    assert returned_cells.cell_label.tolist() == CELLS
    reconciliation = fetch_provenance['gene_axis_reconciliation']
    assert reconciliation['selected_measured_rows'] == 4
    assert reconciliation['metadata_alias_rows_collapsed'] == 1


def test_hmba_macaque_prefers_ensembl_identifier_over_duplicate_symbol(tmp_path):
    values = np.array([[5, 6], [7, 8]], dtype=float)
    path = h5ad(
        tmp_path / 'macaque-ensembl-axis.h5ad',
        encoding='dense',
        values=values,
        cells=['c1', 'c2'],
        genes=['ENSMMUG000001', 'ENSMMUG000002'],
    )
    raw = pd.DataFrame({
        'gene_identifier': ['NCBIGene:1', 'NCBIGene:1'],
        'ensembl_gene_identifier': [
            'ENSEMBL:ENSMMUG000001', 'ENSEMBL:ENSMMUG000002'
        ],
        'gene_symbol': ['GTF2IRD1', 'GTF2IRD1'],
        'ncbi_gene_symbol': ['GTF2IRD1', 'GTF2IRD1'],
    })
    metadata, _ = b.normalize_gene_metadata(raw, 'HMBA-Macaque')
    matrix, _, feature = b.subset_h5ad(
        path,
        pd.DataFrame({'cell_label': ['c2', 'c1']}),
        metadata.iloc[[1, 0]].reset_index(drop=True),
        return_feature_matches=True,
    )
    np.testing.assert_array_equal(
        matrix.toarray(),
        np.array([[8, 7], [6, 5]], dtype=float),
    )
    assert feature['methods'] == [
        'var_names/ensembl_gene_identifier',
        'var_names/ensembl_gene_identifier',
    ]
