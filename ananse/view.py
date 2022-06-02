import pandas as pd


def get_binding_tfs(binding, all_tfs=False):
    """
    Returns a list of unique TFs found in the ananse binding output file.

    The binding file contains multiple lists of TFs:
    - set(pd.read_hdf(binding, key="_factor_activity")["factor"])
    - set(x for x in dir(hdf.root) if not x.startswith("_"))
    The former contains all TFs from the motif2factors.txt, while the latter
    contains only those for which a binding probability was calculated.

    The outersection of these sets had an error in "predict_proba()".
    """
    if all_tfs:
        # all TFs (from the motifs2factors.txt used in ananse binding)
        tfs = set(pd.read_hdf(binding, key="_factor_activity")["factor"])
    else:
        # TFs with binding probabilities
        hdf = pd.HDFStore(binding, "r")
        # TODO: This is hacky (depending on "_"), however the hdf.keys() method is
        #       much slower. Currently all TF names do *not* start with "_"
        keys = hdf.root.__members__
        tfs = set(k for k in keys if not k.startswith("_"))
        hdf.close()
    return list(tfs)


def view_h5(
    fname,
    tfs=None,
    regions=None,
    fmt="wide",
    n=None,
    list_regions=False,
    list_tfs=False,
):
    """Extract information from an ANANSE binding.h5 file.

    Parameters
    ----------
    fname : str
        File name (binding.h5).
    tfs : list, optional
        List of transcription factor names to extract. All TFs are shown
        by default.
    regions : list, optional
        List of regions to extract. All regions are shown by default.
    fmt : str, optional
        Return output in 'wide' or in 'long' format. Default is 'wide'.
    n : int, optional
        Return the first n regions and tfs. All are shown by default.
    list_regions : bool, optional
        Return a list of regions
    list_tfs : bool, optional
        Return a list of transcription factors

    Returns
    -------
    pandas.DataFrame
    """
    if list_regions:
        reg = pd.read_hdf(fname, key="_index")
        return pd.DataFrame({"region": list(set(reg.index))})

    if tfs is None:
        tfs = get_binding_tfs(fname)

    if list_tfs:
        return pd.DataFrame({"factor": list(tfs)})

    if fmt not in ["wide", "long"]:
        raise ValueError("fmt should be either 'wide' or 'long'")

    with pd.HDFStore(fname, "r") as hdf:
        if n:
            n = int(n)
            tfs = tfs[: min(len(tfs), n)]

        idx = hdf.get("_index").index
        if n:
            rows = [True for _ in range(n)] + [False for _ in range(n, len(idx))]
            idx = idx[: min(len(idx), n)]
        elif regions:
            rows = idx.isin(regions)
            idx = set(regions) & set(idx)
        else:
            rows = [True for _ in range(len(idx))]
        df = pd.DataFrame(index=idx)
        df.index.name = "loc"
        for tf in tfs:
            df[tf] = hdf.get(tf)[rows].values

    if fmt == "long":
        df = df.reset_index().melt(
            id_vars=["loc"], value_name="prob", var_name="factor"
        )

    return df
