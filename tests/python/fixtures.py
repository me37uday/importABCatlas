"""Synthetic H5AD writers. No fixture in this module is real Allen data."""
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy import sparse

CELLS=['182941331246012878296807398333956011710','9007199254740993','cell-3']
GENES=['ENSMUSG000000001','ENSMUSG000000002','NCBIGene:712737','gene-4']
VALUES=np.array([[1,0,3,0],[0,2,0,4],[5,6,0,0]],dtype=float)

def h5ad(path,encoding='csr',values=None,cells=None,genes=None,categorical=False,obs_column=None):
    values=VALUES if values is None else np.asarray(values)
    cells=CELLS if cells is None else cells
    genes=GENES if genes is None else genes
    text=h5py.string_dtype('utf-8')
    with h5py.File(path,'w') as f:
        f.attrs['encoding-type']='anndata';f.attrs['encoding-version']='0.1.0'
        for group,ids in [('obs',cells),('var',genes)]:
            g=f.create_group(group);g.attrs['_index']='_index'
            g.attrs['encoding-type']='dataframe';g.attrs['encoding-version']='0.2.0'
            g.attrs['column-order']=np.asarray([],dtype=text)
            if categorical:
                z=g.create_group('_index');z.attrs['encoding-type']='categorical'
                z.create_dataset('categories',data=np.array(list(reversed(ids)),dtype=text))
                z.create_dataset('codes',data=np.arange(len(ids)-1,-1,-1,dtype='i4'))
            else:
                d=g.create_dataset('_index',data=np.asarray(ids,dtype=text))
                d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0'
        if obs_column is not None:
            f['obs'].create_dataset('cell_label',data=np.array(obs_column,dtype=text))
        if encoding=='dense':
            d=f.create_dataset('X',data=values);d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
        else:
            m=sparse.csr_matrix(values) if encoding=='csr' else sparse.csc_matrix(values)
            g=f.create_group('X');g.attrs['encoding-type']=encoding+'_matrix';g.attrs['encoding-version']='0.1.0'
            g.attrs['shape']=m.shape
            for key in ['data','indices','indptr']:g.create_dataset(key,data=getattr(m,key))
    return Path(path)

class FakeCache:
    _local=True
    current_manifest='releases/20260711/manifest.json'
    def __init__(self,root):self.root=Path(root);self.assets={};self.requested=[]
    @property
    def list_directories(self):return sorted({d for d,kind,name in self.assets})
    def add(self,d,name,data,kind='metadata'):
        path=self.root/f'file-{len(self.assets)}.{"csv" if kind=="metadata" else "h5ad"}'
        if kind=='metadata':pd.DataFrame(data).to_csv(path,index=False)
        else:path=Path(data)
        self.assets[(d,kind,name)]=path;return path
    def _files(self,d,kind):
        if d not in self.list_directories:raise ValueError(d)
        return [n for dd,k,n in self.assets if dd==d and k==kind]
    def list_metadata_files(self,d):return self._files(d,'metadata')
    def list_expression_matrix_files(self,d):return self._files(d,'expression_matrices')
    def list_image_volume_files(self,d):return self._files(d,'image_volumes')
    def list_mapmycells_files(self,d):return self._files(d,'mapmycells')
    def get_file_path(self,directory,file_name,**kwargs):
        self.requested.append((directory,file_name))
        matches=[p for (d,k,n),p in self.assets.items() if d==directory and n==file_name]
        if len(matches)!=1:raise FileNotFoundError((directory,file_name))
        return matches[0]
    def get_directory_expression_matrix_size(self,d):return '0.001 MB (synthetic)'
