from scsnv import reader


def test_read_vcf():
    vcf = reader.read_vcf("tests/data/test.vcf")
    print(vcf)