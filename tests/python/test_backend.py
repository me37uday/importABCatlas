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
    if spec.get('file_selector'):
        cellname=spec['file_selector']+'/'+cellname;genename=spec['file_selector']+'/'+genename
    c.add(d,cellname,{'cell_label':CELLS,'cluster_alias':['991','42','991'],
                     'feature_matrix_label':['matrix']*3,'donor_label':['001','002','001'],
                     'x':[1.,2.,3.],'y':[2.,3.,4.]})
    c.add(d,'donor',{'donor_label':['001','002'],'diagnosis':['Unknown','AD']})
    c.add(g,genename,{'gene_identifier':GENES,'gene_symbol':['Gad1','Aqp4','X','Y']})
    if t:c.add(t,'cluster_to_cluster_annotation_membership',{
      'cluster_alias':['991','42'],'cluster_annotation_term_set_name':['class']*2,'cluster_annotation_term_name':['Neuron','Non-neuron']})
    if spec.get('mapping_directory') and spec['mapping_directory'] not in c.list_directories:
        c.add(spec['mapping_directory'],'_placeholder',{'x':['1']})
    if spec.get('special')=='aging':
        c.add(d,'cell_cluster_annotations',{'cell_label':CELLS,'cluster_name':['native1','native2','native1'],'cluster_age_bias':['aged','adult','aged']})
        c.add(spec['mapping_directory'],'cell_cross_mapping_annotations',{'cell_label':CELLS,'wmb_cluster_alias':['990','41','990'],'wmb_class_name':['Neuron','Glia','Neuron']})
    for ed in spec['expression_directories']:
        p=h5ad(tmp_path/(dataset.replace('/','_')+str(len(c.assets))+'.h5ad'),encoding,values=value_override)
        c.add(ed,'matrix/log2' if spec.get('imputed') else 'matrix/raw',p,'expression_matrices')
    return c

@pytest.mark.parametrize('dataset',list(b.REGISTRY))
def test_every_route_load_and_fetch_synthetic(tmp_path,dataset):
    c=route_fixture(tmp_path,dataset)
    cm,gd,prov=b.load_data(c,dataset,cell_ids=[CELLS[2],CELLS[0]])
    spec=prov['spec'];kind='log2' if spec.get('imputed') else 'raw'
    # Multi-directory WMB fixtures deliberately contain the same label: select a
    # source explicitly just as the API requires for genuinely ambiguous input.
    fmap={'matrix':{'directory':spec['expression_directories'][0],'file_name':'matrix/'+kind}}
    m,cell,genes,p=b.fetch_data(c,spec,cm,gd,genes=[GENES[3],GENES[0]],data_type=kind,file_map=fmap)
    np.testing.assert_array_equal(m.toarray(),VALUES[np.ix_([2,0],[3,0])].T)
    assert cell.cell_label.tolist()==[CELLS[2],CELLS[0]]
    assert cell.diagnosis.tolist()==['Unknown','Unknown']
    assert p['assay_kind']==('imputed_log2' if kind=='log2' else 'counts')

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
    s=b.resolve_spec(c,'HMBA_Human',{'cell_file':'human/cell_metadata','gene_file':'human/gene'})
    assert s['cell_file']=='human/cell_metadata'

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

def test_plan_uses_exact_file_size_not_directory_size(tmp_path):
    c=route_fixture(tmp_path,'PMDBS')
    cm,gd,prov=b.load_data(c,'PMDBS')
    # A manifest per-file record is evidence; directory size is not.
    c._manifest={'expression_matrices': {'ASAP-PMDBS-10X': {
      'ASAP-PMDBS-10X': {'raw': {'files': {'h5ad': {
       'relative_path': 'expression_matrices/ASAP-PMDBS-10X/20250331/matrix-raw.h5ad',
       'size_bytes': 123456}}}}}}}
    def forbidden(*args):raise AssertionError('directory size must not be used')
    c.get_directory_expression_matrix_size=forbidden
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes']==123456
    assert result[0]['size_source']=='manifest'
    assert total==123456 and required is None
    assert result[0]['file_name'].endswith('/raw')


def test_plan_unknown_is_not_directory_size(tmp_path):
    c=route_fixture(tmp_path,'WHB');cm,gd,prov=b.load_data(c,'WHB')
    result,total,required=b.make_fetch_plan(c,prov['spec'],cm,gd,'raw')
    assert result[0]['file_bytes'] is None
    assert total is None and required is None


def test_plan_log2_only_when_requested(tmp_path):
    c=route_fixture(tmp_path,'PMDBS');cm,gd,prov=b.load_data(c,'PMDBS')
    c.add('ASAP-PMDBS-10X','matrix/log2',h5ad(tmp_path/'log2.h5ad'),kind='expression_matrices')
    result,_,_=b.make_fetch_plan(c,prov['spec'],cm,gd,'log2')
    assert len(result)==1 and result[0]['file_name'].endswith('/log2')
