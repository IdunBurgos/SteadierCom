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

    def test0(self):
        """Test simple input"""
        df = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test00'
        )
        

        with self.subTest("Produces output"):
            assert df is not None and len(df) > 0
        
    def test1(self):
        """Test co-growth on media with compound identifiers/reaction identifiers"""
        df1 = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test01a'
        )
        
        df2 = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db2.tsv',
            output='tests/output/test01b'
        )

        with self.subTest("Produces output"):
            assert df1 is not None and len(df1) > 0
        
        with self.subTest("Produces output"):
            assert df2 is not None and len(df2) > 0
            
        with self.subTest("Check that they are equal"):
            assert df1.shape==df2.shape
        

    def test2(self):
        """Test that sum of mass rate follows protein constraints (Expected: first - ec_nh4_ko>ec_glc_ko, second- ec_nh4_ko<ec_glc_ko)"""
        df1 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities1.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test02a'
        )
        
        df2 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test02b'
        )
        
        with self.subTest("Produces output"):
            assert df1 is not None and len(df1) > 0
            
        with self.subTest("Produces output"):
            assert df2 is not None and len(df2) > 0
        
        with self.subTest("Sum mass flux ec_nh4_ko>ec_glc_ko for df1"):
            df_sum = df1.groupby(["receiver"]).sum()["mass_rate"]

            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko>ec_glc_ko
            
        with self.subTest("Sum mass flux ec_nh4_ko<ec_glc_ko for df2"):
            df_sum = df2.groupby(["receiver"]).sum()["mass_rate"]
            ec_nh4_ko = df_sum["ec_nh4_ko"]
            ec_glc_ko = df_sum["ec_glc_ko"]
            assert ec_nh4_ko<ec_glc_ko


        
    def test3(self):
        """Test sampling for co-growth on media with fixed relative abundance ec_nh4_ko<ec_glc_ko """
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            sample=10,
            output='tests/output/test03'
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
                    
        
    def test4(self):
        """Test sampling for co-growth on media with fixed relative abundance ec_nh4_ko<ec_glc_ko AND fixed growth rate"""
        df = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            sample=10,
            growth=0.1,
            output='tests/output/test04'
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
            
        
    def test5(self):
        """Test sampling for co-growth on media with variabe relative abundance and fixed growth. First case without specifying members in community, second case with specifying members"""
        df1 = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db.tsv',
            sample=10,
            output='tests/output/test05a',
            growth=0.1,
        )
        
        df2 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities_no_abundance.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            sample=10,
            output='tests/output/test05b',
            growth=0.1,
        )
        
        
        with self.subTest("Produces output"):
            assert df1 is not None and len(df1) > 0
            
        with self.subTest("Produces output"):
            assert df2 is not None and len(df2) > 0
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df1.columns
        
        with self.subTest("Frequency column is present"):
            assert "frequency" in df2.columns
            
        with self.subTest("Abundance dataframes have the same number of samples"):
            df1_abundance = pd.read_csv("tests/output/test05a_abundance.tsv")
            df2_abundance = pd.read_csv("tests/output/test05b_abundance.tsv")
            
            assert df1_abundance.shape==df2_abundance.shape 
            
            
       
    def test6(self):
        """Test co-growth on media with fixed relative abundance and growth - with some unlimited compounds"""
        
        # No prefix for unlimited compounds
        df1 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test07a',
            growth=0.1,
            unlimited='tests/data/unlimited.txt',
        )
        
        # M_ prefix for unlimited compounds
        df2 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test07b',
            growth=0.1,
            unlimited='tests/data/unlimited_M_prefix.txt',
        )
        
        # R_ prefix for unlimited compounds
        df3 = main_run(
            models=['tests/data/*.xml'],
            communities='tests/data/communities2.tsv',
            media='M9',
            mediadb='tests/data/media_db.tsv',
            output='tests/output/test07c',
            growth=0.1,
            unlimited='tests/data/unlimited_R_prefix.txt',
        )
        

        with self.subTest("No prefix for unlimited - Produces output"):
            assert df1 is not None and len(df1) > 0
            
        with self.subTest("M_ prefix for unlimited - Produces output"):
            assert df2 is not None and len(df2) > 0
        
        with self.subTest("R_ prefix for unlimited - Produces output"):
            assert df3 is not None and len(df3) > 0
            
        with self.subTest("The dataframes are equal"):
           
            assert df1.shape==df2.shape and df1.shape==df3.shape
            
            
            
