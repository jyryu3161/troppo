import pytest
import pandas as pd
import os

from troppo.io import detect_gene_id_type, load_expression_csv


class TestDetectGeneIdType:
    def test_entrez_ids(self):
        ids = ['2597', '3156', '5230', '5315', '1234']
        assert detect_gene_id_type(ids) == 'entrez_id'

    def test_ensembl_ids(self):
        ids = ['ENSG00000111640', 'ENSG00000130164', 'ENSG00000102144']
        assert detect_gene_id_type(ids) == 'ensembl_gene_id'

    def test_hgnc_ids(self):
        ids = ['HGNC:1', 'HGNC:2', 'HGNC:3']
        assert detect_gene_id_type(ids) == 'hgnc_id'

    def test_symbol_ids(self):
        ids = ['GAPDH', 'TP53', 'BRCA1', 'EGFR', 'MYC']
        assert detect_gene_id_type(ids) == 'symbol'

    def test_empty_list(self):
        assert detect_gene_id_type([]) == 'symbol'

    def test_mixed_ids_majority_vote(self):
        # Mostly entrez with one symbol
        ids = ['2597', '3156', '5230', '5315', 'GAPDH']
        assert detect_gene_id_type(ids) == 'entrez_id'


class TestLoadExpressionCsv:
    def test_load_csv(self, tmp_path):
        csv_file = tmp_path / "expression.csv"
        csv_file.write_text("gene,sample1,sample2\nGAPDH,10.5,8.2\nTP53,5.3,7.1\nBRCA1,3.2,4.5\n")
        df = load_expression_csv(str(csv_file))
        assert df.shape == (3, 2)
        assert 'sample1' in df.columns
        assert 'GAPDH' in df.index

    def test_load_tsv(self, tmp_path):
        tsv_file = tmp_path / "expression.tsv"
        tsv_file.write_text("gene\tsample1\tsample2\nGAPDH\t10.5\t8.2\nTP53\t5.3\t7.1\n")
        df = load_expression_csv(str(tsv_file))
        assert df.shape == (2, 2)
        assert 'sample1' in df.columns

    def test_auto_detect_gene_id(self, tmp_path):
        csv_file = tmp_path / "expression.csv"
        csv_file.write_text("gene,s1\n2597,10.5\n3156,8.2\n5230,12.1\n")
        df = load_expression_csv(str(csv_file))
        assert df.index.name == 'entrez_id'

    def test_select_sample_cols(self, tmp_path):
        csv_file = tmp_path / "expression.csv"
        csv_file.write_text("gene,s1,s2,s3\nGAPDH,10.5,8.2,3.1\nTP53,5.3,7.1,2.4\n")
        df = load_expression_csv(str(csv_file), sample_cols=['s1', 's3'])
        assert list(df.columns) == ['s1', 's3']
