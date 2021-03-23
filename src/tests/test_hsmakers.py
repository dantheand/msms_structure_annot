"""Tests for the hsmakers module of msms_structure_annot
"""

import pytest
import pandas as pd
from msms_structure_annot.hsmakers import frag_hs, mk_charge_df
from msms_structure_annot.paths import test_data_dir

# Import the pickled dataframes with example test data to compare against
ptms_df = pd.read_pickle(test_data_dir / 'ptms_df.pkl')
hs_df = pd.read_pickle(test_data_dir / 'hs_df.pkl')
frag_df = pd.read_pickle(test_data_dir / 'frag_df.pkl')
frag_df_charged = pd.read_pickle(test_data_dir / 'frag_df_charged.pkl')

parent_seq = 'GGGG'
N_term_mod = 0
C_term_mod = 18.0027
proton_m = 1.0078
charges = [1,2,3]

def test_frag_hs():
    """Test to make sure the fragmentation function still works
    """
    expected_result = frag_df
    result = frag_hs(hs_df, ptms_df, parent_seq, N_term_mod, C_term_mod)

    assert expected_result.equals(result)

def test_mk_charge_df():
    """Test to make sure the calculation of multiply-charges species still works
    """
    expected_result = frag_df_charged
    result = mk_charge_df(frag_df, charges, proton_m)

    assert expected_result.equals(result)
