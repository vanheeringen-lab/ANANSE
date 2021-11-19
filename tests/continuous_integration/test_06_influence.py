# from ananse.influence import read_expression
#
#
# def test_read_expression():
#     res = read_expression("tests/data/dge.tsv")
#     assert set(res.keys()) == {"ANPEP", "CD24", "COL6A3", "DAB2", "DMKN"}
#     assert res["ANPEP"].score - 7.44242618323665 < 0.001
#     assert res["ANPEP"].realfc - 7.44242618323665 < 0.001
#     assert res["ANPEP"].absfc - 7.44242618323665 < 0.001
#     assert res["COL6A3"].score == 0
#     assert res["COL6A3"].realfc - 11.0553152937569 < 0.001
#     assert res["COL6A3"].absfc - 11.0553152937569 < 0.001
#
#
# test_read_expression()
