import os.path

import numpy as np
import pandas as pd
from scipy import stats

from ananse.utils import cleanpath


class Distributions:
    def __init__(self):
        # dist_functions = [f for f in dir(ananse.distributions) if f.endswith("_dist")]
        dist_functions = [scale_dist, log_scale_dist, scipy_dist, peak_rank_dist, peak_rank_file_dist]
        self.functions = {func.__name__: func for func in dist_functions}

    def get(self):
        """list distribution methods"""
        return list(self.functions.keys())

    def set(self, dist_func):
        """return a distribution method by name"""
        dist_functions = self.get()
        if dist_func not in dist_functions:
            raise ValueError(
                f"Distribution function '{dist_func}' not recognised. Options: {', '.join(dist_functions)}"
            )
        return self.functions[dist_func]


def scale_dist(scores, **kwargs):
    """
    Scale the scores between 0 and 1
    """
    return (scores - np.min(scores)) / (np.max(scores) - np.min(scores))


def log_scale_dist(scores, **kwargs):
    """
    Scale the log of the scores between 0 and 1
    """
    scores = np.log(scores+1)
    return (scores - np.min(scores)) / (np.max(scores) - np.min(scores))


def replace_infs(dist):
    """
    Replace positive and negative infinity with the closes real value in the array
    """
    # https://stackoverflow.com/questions/12937824/lognormal-random-numbers-centered-around-a-high-value
    min_real_val = np.nanmin(dist[dist != -np.inf])
    dist[dist == -np.inf] = min_real_val
    max_real_val = np.nanmax(dist[dist != np.inf])
    dist[dist == np.inf] = max_real_val
    return dist


def scipy_dist(scores, **kwargs):
    """
    fit scores to a scipy.stats distribution.
    specified distribution name via kwargs['dist']
    """
    scores = scores + 1  # add pseudocount
    x = range(len(scores))

    dist_name = kwargs.get("dist", "lognorm")
    if dist_name not in dir(stats):
        raise ValueError(f"'{dist_name}' is not a recognized scipy.stats model.")
    distribution = getattr(stats, dist_name)  # eval(f"stats.{dist_name}")

    # fit dist to data
    params = distribution.fit(scores)

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Calculate fitted PDF
    dist = distribution.pdf(x, loc=loc, scale=scale, *arg)
    dist = replace_infs(dist)
    return dist


# def lognorm_dist(scores, **kwargs):
#     """
#     fit scores to a log normal distribution
#     """
#     scores = scores + 1  # add pseudocount
#     x = range(len(scores))
#
#     # mu = np.log(scores).mean()
#     # sigma = np.log(scores).std()
#     # dist = stats.lognorm([sigma], loc=mu).pdf(x)
#
#     s, loc, scale = stats.lognorm.fit(scores)  # floc=0
#     dist = stats.lognorm.pdf(x=x, s=s, loc=loc, scale=scale)
#     return dist


def peak_rank_dist(scores, **kwargs):
    """
    Fit scores to a distribution similar to what the p300 model was trained on
    """
    # use a lognormal distribution:
    # https://github.com/jsh58/Genrich#p-value-calculation
    # # peak_rank_file = "ananse/db/peak_rank.txt"
    # # scores = pd.read_csv(peak_rank_file, header=None)[0]
    # # mu = np.log(scores+1).mean()
    # # sigma = np.log(scores+1).std()
    # mu = 1.0500836750482117
    # sigma = 0.8000981267240566
    #
    # x = len(scores)
    # rng = np.random.default_rng(seed=None)
    # dist = rng.lognormal(mean=mu, sigma=sigma, size=x)
    #
    # print("proximity to the initial distribtion")
    # print("delta mu:", np.abs(mu - np.log(dist).mean()))
    # print("delta std:", np.abs(sigma - np.log(dist).std()))

    # best fitting distribution turns out to be this loglaplace
    x = range(len(scores))
    c = 0.92
    loc = 1.00
    scale = 1.14
    dist = stats.loglaplace.pdf(x=x, c=c, loc=loc, scale=scale)
    dist = replace_infs(dist)
    return dist


def peak_rank_file_dist(scores, **kwargs):
    """
    fit scores to the distribution in kwargs['file'].
    builtin files: "peak_rank.txt" and "peak_rank_hg38_h3k27ac.txt"
    """
    dist_filename = kwargs.get("file", "peak_rank.txt")

    # internal data or user data
    if dist_filename in ["peak_rank.txt", "peak_rank_hg38_h3k27ac.txt"]:
        package_dir = os.path.dirname(__file__)
        dist_filepath = os.path.join(package_dir, "db", dist_filename)
    else:
        dist_filepath = cleanpath(dist_filename)

    if not os.path.exists(dist_filepath):
        raise ValueError(f"Could not find file {dist_filepath}")

    dist = pd.read_csv(dist_filepath, header=None)
    n = scores.shape[0]
    max_n = dist.shape[0]
    if max_n < n:
        raise ValueError(f"Too many regions ({n}) to fit to '{dist_filename}' ({max_n})")

    dist = dist.sample(n=n, random_state=1)[0].tolist()
    return dist
