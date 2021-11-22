import os
import ananse.network
from ananse.utils import check_path, load_tfs, load_regions, check_cores
from dask.distributed import Client, LocalCluster
from loguru import logger


@logger.catch
def network(args):
    ncore = check_cores(args.ncore)
    memory_limit = "16GB"
    if ncore == 1:
        # With one core more memory is needed
        memory_limit = "20GB"

    b = ananse.network.Network(
        genome=args.genome,  # checked in CLI
        gene_bed=check_path(args.annotation),
        include_promoter=args.include_promoter,
        include_enhancer=args.include_enhancer,
        full_output=args.full_output,
    )

    cluster = LocalCluster(
        local_directory=os.environ.get("TMP", None),
        scheduler_port=0,
        dashboard_address=None,  # noqa
        n_workers=ncore,
        threads_per_worker=2,
        memory_limit=memory_limit,
    )
    client = Client(cluster)

    b.run_network(
        binding=check_path(args.binding),
        fin_expression=check_path(args.fin_expression),
        tfs=load_tfs(args.tfs),
        regions=load_regions(args.regions, args.genome, os.path.basename(args.outfile)),
        outfile=check_path(args.outfile, error_missing=False),
    )
    client.close()
