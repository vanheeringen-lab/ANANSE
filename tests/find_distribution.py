# source:
# https://stackoverflow.com/questions/6620471/fitting-empirical-distribution-to-theoretical-ones-with-scipy-python

import warnings
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as st
from tqdm import tqdm

mpl.rcParams['figure.figsize'] = (16.0, 12.0)
plt.style.use('ggplot')


# Create models from data
def best_fit_distribution(data, bins=200, ax=None):
    """Find the best fitting distribution to the data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Distributions to check
    DISTRIBUTIONS = [
        st.alpha,
        st.anglit,
        st.arcsine,
        st.argus,
        st.beta,
        st.betaprime,
        st.bradford,
        st.burr,
        st.burr12,
        st.cauchy,
        st.chi,
        st.chi2,
        st.cosine,
        st.crystalball,
        st.dgamma,
        st.dweibull,
        st.erlang,
        st.expon,
        st.exponnorm,
        st.exponweib,
        st.exponpow,
        st.f,
        st.fatiguelife,
        st.fisk,
        st.foldcauchy,
        st.foldnorm,
        st.genlogistic,
        st.gennorm,
        st.genpareto,
        st.genexpon,
        st.genextreme,
        st.gausshyper,
        st.gamma,
        st.gengamma,
        st.genhalflogistic,
        st.geninvgauss,
        st.gilbrat,
        st.gompertz,
        st.gumbel_r,
        st.gumbel_l,
        st.halfcauchy,
        st.halflogistic,
        st.halfnorm,
        st.halfgennorm,
        st.hypsecant,
        st.invgamma,
        st.invgauss,
        st.invweibull,
        st.johnsonsb,
        st.johnsonsu,
        st.kappa4,
        st.kappa3,
        st.ksone,
        st.kstwo,
        st.kstwobign,
        st.laplace,
        st.laplace_asymmetric,
        st.levy,
        st.levy_l,
        # st.levy_stable,  # unstable in v1.6.0
        st.logistic,
        st.loggamma,
        st.loglaplace,
        st.lognorm,
        st.loguniform,
        st.lomax,
        st.maxwell,
        st.mielke,
        st.moyal,
        st.nakagami,
        st.ncx2,
        st.ncf,
        st.nct,
        st.norm,
        st.norminvgauss,
        st.pareto,
        st.pearson3,
        st.powerlaw,
        st.powerlognorm,
        st.powernorm,
        st.rdist,
        st.rayleigh,
        st.rice,
        st.recipinvgauss,
        st.semicircular,
        st.skewnorm,
        st.t,
        st.trapezoid,
        st.triang,
        st.truncexpon,
        st.truncnorm,
        st.tukeylambda,
        st.uniform,
        # st.vonmises,  # does not work in v1.6.0
        st.vonmises_line,
        st.wald,
        st.weibull_min,
        st.weibull_max,
        st.wrapcauchy,
    ]

    # Best holders
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in tqdm(DISTRIBUTIONS):
        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # if ax is passed, add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(label=distribution.name, legend=True, ax=ax)
                except Exception:
                    pass

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse

        except Exception:
            pass

    return best_distribution.name, best_params


def make_pdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf


def find_best_pdf(data, outfile=None):
    # Plot for comparison
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Find best fit distribution
    best_fit_name, best_fit_params = best_fit_distribution(data, 200, ax1)
    best_dist = getattr(st, best_fit_name)

    data.plot(kind='hist', density=True, bins=50, alpha=0.5, label='Data', legend=True, ax=ax1,
              color=mpl.rcParams['axes.prop_cycle'].by_key()['color'][1])

    # Save plot limits
    dataYLim = ax1.get_ylim()
    dataXLim = ax1.get_xlim()

    # Update plots
    ax1.set_ylim(dataYLim)
    ax1.set_xlim(dataXLim)
    ax1.set_title('All Fitted Distributions\n')
    ax1.set_xlabel('Score')

    # Make PDF with best params
    pdf = make_pdf(best_dist, best_fit_params)

    # Display
    pdf.plot(lw=2, label='PDF', legend=True, ax=ax2)
    data.plot(kind='hist', density=True, bins=50, alpha=0.5, label='Data', legend=True, ax=ax2)

    param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
    param_str = ', '.join(['{}={:0.2f}'.format(k, v) for k, v in zip(param_names, best_fit_params)])
    dist_str = '{}({})'.format(best_fit_name, param_str)

    ax2.set_ylim(dataYLim)
    ax2.set_xlim(dataXLim)
    ax2.set_title('Best fit distribution \n' + dist_str)
    ax2.set_xlabel('Score')

    if outfile:
        fig.savefig(outfile)
    fig.show()


# Load data
peak_rank_file = "db/peak_rank.txt"  # "ananse/db/peak_rank.txt"
scores = pd.read_csv(peak_rank_file, header=None)[0]
data = pd.Series(scores + 1)
outfile = "../tests/output/distributions.png"  # "tests/output/distributions.png"

# run
find_best_pdf(data, outfile)
