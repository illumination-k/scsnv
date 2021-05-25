from scsnv import reader


def test_read_vcf():
    vcf = reader.read_vcf("tests/data/test.vcf")
    print(vcf)

def test_read_mm():
    adata = reader.read_mm("tests/data/test_consensus.mtx", "tests/data/test.vcf")
    print(adata) 