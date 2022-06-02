from loguru import logger
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.preprocessing import minmax_scale
from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


def plot_influence(infile, outfile):
    """Plot TF influence score to expression."""

    df = pd.read_table(infile, sep="\t", index_col="factor")
    df = df.dropna()
    factors = list(df.sort_values("influence_score").tail(20).index)
    xcol = "factor_fc"
    plt.figure(figsize=(8, 6))
    sns.regplot(
        data=df,
        x=xcol,
        y="influence_score",
        fit_reg=False,
        scatter_kws={"s": df["direct_targets"] / 10, "alpha": 0.5},
    )
    x = df.loc[factors, xcol]
    y = df.loc[factors, "influence_score"]
    texts = []
    for s, xt, yt in zip(factors, x, y):
        texts.append(plt.text(xt, yt, s))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black"))
    plt.xlabel("Log2 fold change of TF")
    plt.ylabel("Influence score")
    plt.savefig(outfile, dpi=300)
    plt.clf()


def read_diff_network(diff_network_file, full_output):
    """read the differential network file outputed by the influence command."""
    data_columns = ["source", "target", "weight"]
    if full_output:
        data_columns = None  # all
        # data_columns = [
        #     "source",
        #     "target",
        #     "weight",
        #     "weight_target",
        #     "weight_source",
        #     "tf_expression_target",
        #     "tf_expression_source",
        #     "target_expression_target",
        #     "target_expression_source",
        #     "weighted_binding_target",
        #     "weighted_binding_source",
        #     "tf_activity_target",
        #     "tf_activity_source",
        # ]
    rnet = pd.read_csv(
        diff_network_file,
        sep="\t",
        usecols=data_columns,
        dtype="float64",
        converters={"source": str, "target": str},
    )
    grn = nx.from_pandas_edgelist(rnet, edge_attr=True, create_using=nx.DiGraph)
    return grn


def plot_grn(
    infile,
    fname,
    outfile,
    edge_info="weight",
    network_algorithm="neato",
    n_tfs=20,
    cmap="viridis",
    edge_min=0.1,
    full_output=False,
):
    """
    Plot the top differential TFs and their interactions with the highest score into a GRN network image.

    Parameters
    ----------
    infile: str
        influence output file
    fname: str
        diffnetwork text file
    outfile: str
        output file location
    edge_info: str, optional
        column to use for interaction weight, default is 'weight'. Full output diff options:
        weight, weight_target, weight_source, tf_expression_target, tf_expression_source,
        target_expression_target, target_expression_source, weighted_binding_target, weighted_binding_source,
        tf_activity_target, tf_activity_source.
    network_algorithm: str , optional
        pyviz cluster algorithm used for node placement, options include: neato, dot, fdp, twopi, sfdp, circo
    n_tfs: int, optional
        number of (significantly differential expressed) Tfs to visualize
    cmap: str, optional
        matplotlib colour pallet used
    edge_min: float, optional
        minimum value the selected edge value needs to have to be included in the GRN,
        0.1 seems like a good number for weight.
    full_output: bool, optional
        required to use edge_info other than the default 'weight'.
    """
    if full_output is False and edge_info != "weight":
        logger.info(
            "selected custom edge weight, but missing full output option, selecting 'weight' as edge weight instead."
        )
        edge_info = "weight"

    # select the top TFs:
    df = pd.read_table(
        infile,
        sep="\t",
        index_col="factor",
        usecols=["factor", "influence_score", "G_score_scaled"],
    )
    df = df[df.G_score_scaled > 0]  # plot only TFs that are differential
    top_factors = list(df.sort_values("influence_score").tail(n_tfs).index)

    if len(df) == 0:
        logger.warning(f"No differential TFs in {infile}!")
        return

    # read the diff_network file into a network
    grn = read_diff_network(fname, full_output)
    # filter the diffnetwork to only contain topTF-topTF
    tf_grn = nx.subgraph_view(grn, filter_node=lambda tf: tf in top_factors)
    # filter the diffnetwork to only contain edges above the cuttoff value
    tf_grn2 = nx.DiGraph(
        ((u, v, e) for u, v, e in tf_grn.edges(data=True) if e[edge_info] > edge_min)
    )
    # make the network directed again after filtering
    tf_grn2 = tf_grn2.to_directed()
    # remove TFs with no interactions
    tf_grn2.remove_nodes_from(list(nx.isolates(tf_grn2)))
    # load all edge info for scaling edge width
    edge_atribute = list(nx.get_edge_attributes(tf_grn2, edge_info).values())
    if len(edge_atribute) == 0:
        raise ValueError("Top factors do not share any connections.")
    edge_atribute_scaled = minmax_scale(
        edge_atribute, feature_range=(0, 1), axis=0, copy=True
    )
    # calculate edge weight quantile regions to dipslay 4 relevant numbers within the edge legend
    # normzalize the edge atributes (interaction scores) to be between 0 and 1 (where 0-1 corresponds to edge width)
    edges_norm_weight = [
        0,
        round((np.quantile(sorted(edge_atribute_scaled), 1) / 4), 3),
        round((np.quantile(sorted(edge_atribute_scaled), 1) / 2), 3),
        round(np.quantile(sorted(edge_atribute_scaled), 1), 3),
    ]
    # Also get the numbers of the unnormalized edge numbers to display within the legend next to the lines
    edges_weight = [
        0,
        round((np.quantile(sorted(edge_atribute), 1) / 4), 3),
        round((np.quantile(sorted(edge_atribute), 1) / 2), 3),
        round(np.quantile(sorted(edge_atribute), 1), 3),
    ]
    # lets calculate the nodes their outdegree (edges regulating other TFs):
    outdegree = pd.DataFrame(tf_grn2.out_degree(weight=edge_info))
    outdegree = outdegree[1]
    node_outdegree_size = 600 + outdegree * 100

    # Lets set some plotstuff
    colors = outdegree
    vmin = min(colors)
    vmax = max(colors)
    cmap = plt.get_cmap(cmap)
    plt.figure(figsize=(10, 10))
    # calculate node position of the graph
    pos = nx.drawing.nx_agraph.graphviz_layout(tf_grn2, prog=network_algorithm)
    # plot the TF nodes
    nx.draw_networkx_nodes(
        tf_grn2, pos, node_size=node_outdegree_size, node_color=outdegree, cmap=cmap
    )
    # plot the TF name labels:
    nx.draw_networkx_labels(
        tf_grn2,
        pos,
        font_color="black",
        bbox=dict(facecolor="white", alpha=0.5, pad=0),
    )
    # plot the node colour legend:
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cbar = plt.colorbar(sm)
    cbar.ax.set_ylabel("outdegree (regulation other TFs)", rotation=270, labelpad=25)
    # plot the TF-TF edges:
    nx.draw_networkx_edges(
        tf_grn2,
        pos,
        arrows=True,
        arrowstyle="->",
        arrowsize=20,
        width=edge_atribute_scaled,
        node_size=node_outdegree_size,
        connectionstyle="arc3, rad = 0.1",
    )
    # add edge width legend:
    lines = []
    for _i, width in enumerate(edges_norm_weight):
        lines.append(Line2D([], [], linewidth=width, color="black"))
    plt.legend(
        lines,
        edges_weight,
        bbox_to_anchor=(0, 0.5),
        frameon=False,
        title=f"{edge_info}",
    )
    # save the plot
    plt.savefig(outfile, dpi=300)
    plt.clf()
