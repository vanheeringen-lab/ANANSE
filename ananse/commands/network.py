import os
from tempfile import gettempdir

from dask.distributed import Client, LocalCluster, system
from dask.utils import parse_bytes
from loguru import logger

import ananse.network
from ananse.utils import check_cores, check_path, load_regions, load_tfs


@logger.catch
def network(args):
    ncore = check_cores(args.ncore)
    # actual memory usage depends on the input data
    memory_limit = parse_bytes("20GB")
    memory_limit = max(system.MEMORY_LIMIT, memory_limit)

    b = ananse.network.Network(
        genome=args.genome,  # checked in CLI
        gene_bed=check_path(args.annotation),
        include_promoter=args.include_promoter,
        include_enhancer=args.include_enhancer,
        full_output=args.full_output,
    )

    cluster = LocalCluster(
        local_directory=gettempdir(),
        scheduler_port=0,
        dashboard_address=None,  # noqa: disable dashboard
        n_workers=ncore,
        threads_per_worker=2,
        memory_limit=memory_limit,
    )
    client = Client(cluster)

    b.run_network(
        binding=check_path(args.binding),
        fin_expression=check_path(args.fin_expression),
        column=args.column,
        tfs=load_tfs(args.tfs),
        regions=load_regions(args.regions, args.genome, os.path.basename(args.outfile)),
        outfile=check_path(args.outfile, error_missing=False),
    )
    client.close()
