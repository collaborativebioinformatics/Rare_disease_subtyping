import os
import argparse

from nvflare.job_config.api import FedJob
from nvflare.recipe import SimEnv
from nvflare.app_common.workflows.scatter_and_gather import ScatterAndGather

from fed_kmeans_components import (
    FederatedKMeansExecutor,
    KMeansAggregator,
    KMeansShareableGenerator,
)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--k", type=int, default=6)
    p.add_argument("--num_rounds", type=int, default=25)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--data_dir", type=str, default="clients_data")
    p.add_argument("--init_site", type=str, default="EUR")  # which site to bootstrap centers from
    p.add_argument("--workspace", type=str, default="output_nvflare_kmeans")
    p.add_argument("--clients", type=str, default="AFR,AMR,EAS,EUR,SAS")
    p.add_argument("--threads", type=int, default=5)
    p.add_argument("--min_clients", type=int, default=5)
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    k = args.k
    num_rounds = args.num_rounds
    seed = args.seed

    data_dir = os.path.abspath(args.data_dir)
    workspace = os.path.abspath(args.workspace)
    os.makedirs(workspace, exist_ok=True)

    clients = [c.strip() for c in args.clients.split(",") if c.strip()]
    init_csv = os.path.abspath(os.path.join(args.data_dir, f"client_{args.init_site}.csv"))

    job = FedJob(name=f"fed_kmeans_k{k}_r{num_rounds}")

    workflow = ScatterAndGather(
        min_clients=args.min_clients,
        num_rounds=num_rounds,
        train_task_name="kmeans_step",
        aggregator_id="kmeans_agg",
        shareable_generator_id="kmeans_sg",
        allow_empty_global_weights=True,
    )

    job.to_server(workflow, id="kmeans_wf")
    job.to_server(KMeansAggregator(k=k, seed=seed), id="kmeans_agg")
    job.to_server(KMeansShareableGenerator(k=k, init_csv=init_csv, seed=seed), id="kmeans_sg")

    job.to_clients(
        FederatedKMeansExecutor(data_dir=data_dir, k=k, task_name="kmeans_step"),
        id="kmeans_exec",
        tasks=["kmeans_step"],
    )

    env = SimEnv(num_clients=len(clients))

    run = job.simulator_run(
        workspace=workspace,
        clients=clients,
        threads=args.threads,
    )

    print("workspace used:", workspace)