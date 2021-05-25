from typing import Union
import pandas as pd

VCF_HEADERS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def read_vcf(path: str) -> tuple:
    header, content = "", ""
    with open(path) as r:
        for line in r:
            if line.startswith("##"):
                header += line
            else:
                content += line

    return header, line


class VCFDataFrame(pd.DataFrame):
    def __init__(
        self,
        data=None,
        index=None,
        columns=None,
        vcf_header=None,
        dtype=None,
        copy=None,
    ) -> None:
        if columns is not None:
            assert columns[: len(VCF_HEADERS)] == VCF_HEADERS

        super().__init__(
            data=data, index=index, columns=columns, dtype=dtype, copy=copy
        )

        self.vcf_header = vcf_header

    def read_vcf(self, path: str):
        from .reader import read_vcf

    def fetch(self, loc: str) -> pd.DataFrame:
        pass

    def mutation_uname(self) -> pd.Series:
        """
        Returns:
            pd.Series: chr:pos:ref-alt
        """
        pass

    def get_sample_info(self, tag: str) -> pd.DataFrame:
        pass

    def get_format_info(
        self, tags: Union[str, list[str]]
    ) -> Union[pd.DataFrame, pd.Series]:
        pass
