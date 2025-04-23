#!/usr/bin/env python

"""Tests for `steadiercom` package."""


import unittest
import pandas as pd

from steadiercom.cli import main_run

class TestSteadiercom(unittest.TestCase):
    """Tests for `steadiercom` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        pass

    def tearDown(self):
        """Tear down test fixtures, if any."""
        pass

    def test1(self):
        """Test co-growth on media with compound identifiers"""
        df = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test1'
        )

        assert df is not None and len(df) > 0
        
        
    def test2(self):
        """Test co-growth on media with reaction identifiers"""
        df = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            output='tests/output/test2'
        )

        assert df is not None and len(df) > 0
        
    def test3(self):
        """Test that results from compound identifiers and reaction identifiers are equal"""
        df1 = pd.read_csv("tests/output/test1.tsv",sep="\t")
        df2 = pd.read_csv("tests/output/test2.tsv",sep="\t")
        
        assert df1.shape==df2.shape

    def test4(self):
        """Test co-growth on media with relative abundance ec_nh4_ko>ec_glc_ko """
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities1.tsv',
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            output='tests/output/test4'
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
        with self.subTest("Sum mass flux ec_nh4_ko>ec_glc_ko"):
            df_sum = df.groupby(["receiver"]).sum()["mass_rate"]

            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko>ec_glc_ko
        
    def test5(self):
        """Test co-growth on media with relative abundance ec_nh4_ko<ec_glc_ko """
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            output='tests/output/test5'
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
            
        with self.subTest("Sum mass flux ec_nh4_ko<ec_glc_ko"):
            df_sum = df.groupby(["receiver"]).sum()["mass_rate"]
            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko<ec_glc_ko


        
    def test6(self):
        """Test sampling for co-growth on media with relative abundance ec_nh4_ko<ec_glc_ko """
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            sample=10,
            output='tests/output/test6'
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df.columns
            
        with self.subTest("Sum mass flux ec_nh4_ko>ec_glc_ko"):
            df_copy = df.copy()
            df_copy["mass_rate*frequency"] = df_copy.mass_rate*df.frequency
            df_sum = df_copy.groupby(["receiver"]).sum()["mass_rate*frequency"]

            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko<ec_glc_ko
            
        with self.subTest("Check that results for sampling match previous results"):
            df_comparison = pd.read_csv("tests/output/test6.tsv",sep="\t")
            assert df_comparison.shape==df.shape
        
        
    def test7(self):
        """Test sampling for co-growth on media with relative abundance ec_nh4_ko<ec_glc_ko AND fixed growth rate"""
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            sample=10,
            growth=0.1,
            output='tests/output/test7'
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df.columns
            
        with self.subTest("Sum mass flux ec_nh4_ko>ec_glc_ko"):
            df_copy = df.copy()
            df_copy["mass_rate*frequency"] = df_copy.mass_rate*df.frequency
            df_sum = df_copy.groupby(["receiver"]).sum()["mass_rate*frequency"]

            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko<ec_glc_ko
            
        with self.subTest("Check that results for sampling match previous results"):
            df_comparison = pd.read_csv("tests/output/test7.tsv",sep="\t")
            assert df_comparison.shape==df.shape
            
            
        
    def test8(self):
        """Test sampling for co-growth on media with variabe relative abundance and fixed growth"""
        df = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            sample=10,
            output='tests/output/test8',
            growth=0.1,
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df.columns
            
    def test9(self):
        """Test sampling for co-growth on media with variabe relative abundance and fixed growth"""
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities_no_abundance.tsv',
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            sample=10,
            output='tests/output/test9',
            growth=0.1,
        )
        
        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df.columns
            
       