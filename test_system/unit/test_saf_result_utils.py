"""saf_result_utils: reading the published archive's dropped-E tokens."""
from pathlib import Path

import numpy as np

from saf_result_utils import loadtxt_repaired


def test_dropped_e_tokens_truncate_like_matlab(tmp_path):
    # gfortran writes 0.8384675E-101 as 0.8384675-101; MATLAB's load() reads it
    # as 0.8384675, and the published statistics were computed that way
    p = tmp_path / "t.txt"
    p.write_text("1.0 0.8384675-101\n2.0 3.5\n")
    arr, n = loadtxt_repaired(Path(p), usecols=(0, 1))
    assert n == 1
    np.testing.assert_allclose(arr, [[1.0, 0.8384675], [2.0, 3.5]])


def test_paleo_site_stats_converter_same_rule(tmp_path):
    from paleo_site_stats import _loadtxt_repaired
    p = tmp_path / "t.txt"
    p.write_text("1.0 2.0 0.5-101\n")
    arr, n = _loadtxt_repaired(Path(p), usecols=(2,))
    assert n == 1 and float(arr) == 0.5
