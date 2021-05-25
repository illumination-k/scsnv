from anndata import AnnData
from scipy.io import mmread
import pandas as pd

from ._core import VCFDataFrame

# return dataframe with header
def read_vcf(vcf_path: str) -> pd.DataFrame:
    df = pd.read_csv(vcf_path, comment="#", header=None, sep="\t")
    columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    ]
    df.columns = columns
    return df


def read_mm(mat_path: str, vcf_path: str) -> AnnData:
    mm = mmread(mat_path)
    vcf = read_vcf(vcf_path=vcf_path)

    return AnnData(X=mm)
