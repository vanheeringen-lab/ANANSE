# import numpy as np
#
# import pytest
# import ananse.distributions
#
#
# def test_distributions():
#     d = ananse.distributions.Distributions()
#     func_list = d.get()
#     assert isinstance(func_list, list)
#
#     for func in func_list:
#         d.set(func)
#
#
# scores = np.array([0, 1, 2])
#
#
# def test_scale_dist():
#     s = ananse.distributions.scale_dist(scores)
#     assert np.array_equal(s, np.array([0, 0.5, 1]))
#
#
# def test_log_scale_dist():
#     s = ananse.distributions.log_scale_dist(scores)
#     assert np.allclose(s, np.array([0.0, 0.63092975, 1.0]))
#
#
# def test_replace_infs():
#     score_w_infs = [-np.inf, 0, 1, np.inf]
#     s = ananse.distributions.replace_infs(score_w_infs)
#     assert np.array_equal(s, np.array([0, 0, 1, 1]))
#
#
# def test_scipy_dist():
#     s = ananse.distributions.scipy_dist(scores, **{"dist": "loglaplace"})
#     assert np.allclose(s, np.array([4.72219713e-05, 2.05410078e-01, 6.83221921e-01]))
#
#     s = ananse.distributions.scipy_dist(scores, **{"dist": "lognorm"})
#     assert np.allclose(s, np.array([0, 8.0793556e12, 2.8352896e-02]))
#
#     with pytest.raises(ValueError):
#         ananse.distributions.scipy_dist(scores, **{"dist": "wrongname"})
#
#
# def test_peak_rank_dist():
#     s = ananse.distributions.peak_rank_dist(scores)
#     assert np.allclose(s, np.array([0, 0.4077607, 0.4077607]))
#
#
# def test_peak_rank_file_dist():
#     s = ananse.distributions.peak_rank_file_dist(scores, **{"file": "peak_rank.txt"})
#     assert len(s) == 3
#
#     s = ananse.distributions.peak_rank_file_dist(
#         scores, **{"file": "peak_rank_hg38_h3k27ac.txt"}
#     )
#     assert len(s) == 3
#
#     # too many peaks
#     with pytest.raises(ValueError):
#         ananse.distributions.peak_rank_file_dist(
#             range(108_087), **{"file": "peak_rank.txt"}
#         )
