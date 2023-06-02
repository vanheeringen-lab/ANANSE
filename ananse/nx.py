from heapq import heappop, heappush
from itertools import count

import numpy as np


# This piece of code is adapted from the networkx code licensed under a 3-clause license:
# https://networkx.org/documentation/networkx-2.7/#license
def dijkstra_prob_length(G, source, weight, cutoff=None, target=None):  # noqa
    """Uses Dijkstra's algorithm to find shortest weighted paths

    Parameters
    ----------
    G : NetworkX graph
    source : node
        Starting node for paths.
    weight: string
        Name of attribute that contains weight (as a probability between 0 and 1,
        where 0 is the minimum and 1 the maximum).
    target : node label, optional
        Ending node for path. Search is halted when target is found.
    cutoff : integer or float, optional
        Minimum combined, weighted probability.
        If cutoff is provided, only return paths with summed weight >= cutoff.

    Returns
    -------
    paths, distance : dictionaries
        Dictionary of shortest paths keyed by target, and
        a mapping from node to shortest distance to that node from one
        of the source nodes.

    Raises
    ------
    NodeNotFound
        If the `source` is not in `G`.
    """
    # cumulative weight function for values < 1
    weight_function = lambda u, v, data: -np.log10(data[weight])  # noqa
    paths = {source: [source]}

    if cutoff == 0:
        cutoff = None
    if cutoff is not None:
        if cutoff < 0 or cutoff > 1:
            raise ValueError(
                "cutoff (combined, weighted probability) should be between 0 and 1"
            )
        cutoff = -np.log10(cutoff)

    G_succ = G._succ if G.is_directed() else G._adj  # noqa

    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []

    seen[source] = 0
    push(fringe, (0, next(c), source))
    while fringe:
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if v == target:
            break
        for u, e in G_succ[v].items():
            cost = weight_function(v, u, e)
            if cost is None:
                continue

            # This is the major difference, both probability and length are taken
            # into account
            length = len(paths[v])
            vu_dist = dist[v] + cost - np.log10(1 / length)  # noqa
            if length > 1:
                vu_dist = vu_dist + np.log10(1 / (length - 1))

            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            if u in dist:
                u_dist = dist[u]
                if vu_dist < u_dist:
                    raise ValueError("Contradictory paths found:", "negative weights?")
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]  # noqa

    paths = {k: v for k, v in paths.items() if len(v) > 1}
    dist = {k: 10**-v for k, v in dist.items() if k in paths}
    return paths, dist
