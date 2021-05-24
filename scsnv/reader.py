from anndata import Anndata
from scipy.io import mmread
import pandas as pd


def read_vcf(vcf_path: str) -> pd.DataFrame:
    df = pd.read_csv(vcf_path, comment="#", header=None)
    columns = []
    df.columns = columns
    return df

def read_mat(mat_path: str, vcf_path: str) -> Anndata:
    mm = mmread(mat_path)
    vcf = read_vcf(vcf_path=vcf_path)
    return Anndata(X=mm)

