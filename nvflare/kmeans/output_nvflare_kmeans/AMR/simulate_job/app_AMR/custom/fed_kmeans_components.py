import os
import csv
import numpy as np
import pandas as pd

from nvflare.apis.executor import Executor
from nvflare.apis.fl_constant import ReturnCode
from nvflare.apis.shareable import Shareable
from nvflare.apis.dxo import DXO, DataKind, from_shareable

from nvflare.app_common.abstract.aggregator import Aggregator
from nvflare.app_common.abstract.learnable import Learnable
from nvflare.app_common.abstract.shareable_generator import ShareableGenerator


class FederatedKMeansExecutor(Executor):
    """
    Client side:
    receives centers, computes local sums, counts, inertia, returns them.
    It loads the client csv based on the client site name.
    """
    def __init__(self, data_dir: str, k: int, task_name: str = "kmeans_step"):
        super().__init__()
        self.data_dir = data_dir
        self.k = int(k)
        self.task_name = task_name
        self._cache = {}  # site_name -> X

    @staticmethod
    def _to_matrix(df: pd.DataFrame) -> np.ndarray:
        drop_cols = {"patient_id", "superpopulation", "n_variants_carried", "primary_gene", "primary_condition", "primary_subtype", "onset_type","base_impact","allelic_ratio","ancestral_modifier","modifier_type","consequence_modifier","interaction_score","predicted_progression", "all_genes","all_subtypes"}
        feat_cols = [c for c in df.columns if c not in drop_cols]
        return df[feat_cols].to_numpy(dtype=np.float32)

    def _load_X_for_site(self, site_name: str) -> np.ndarray:
        if site_name in self._cache:
            return self._cache[site_name]

        path = os.path.join(self.data_dir, f"client_{site_name}.csv")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing client csv for site={site_name}: {path}")

        df = pd.read_csv(path)
        X = self._to_matrix(df)
        self._cache[site_name] = X
        return X

    def execute(self, task_name: str, shareable: Shareable, fl_ctx, abort_signal):
        if task_name != self.task_name:
            s = Shareable()
            s.set_return_code(ReturnCode.TASK_UNKNOWN)
            return s

        # client/site name
        site = None
        if hasattr(fl_ctx, "get_identity_name"):
            site = fl_ctx.get_identity_name()
        site = site or os.environ.get("NVFLARE_SITE_NAME") or os.environ.get("NVFLARE_CLIENT_NAME")
        if not site:
            s = Shareable()
            s.set_return_code(ReturnCode.BAD_TASK_DATA)
            return s

        X = self._load_X_for_site(site)

        dxo = from_shareable(shareable)
        centers = np.asarray(dxo.data.get("centers", None), dtype=np.float32)

        if centers.ndim != 2 or centers.shape[0] != self.k:
            s = Shareable()
            s.set_return_code(ReturnCode.BAD_TASK_DATA)
            return s

        # squared euclidean
        d2 = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        labels = d2.argmin(axis=1)

        sums = np.zeros_like(centers, dtype=np.float64)
        counts = np.zeros((self.k,), dtype=np.int64)
        for j in range(self.k):
            mask = labels == j
            counts[j] = int(mask.sum())
            if counts[j] > 0:
                sums[j] = X[mask].sum(axis=0)

        inertia = float(d2[np.arange(X.shape[0]), labels].sum())

        out = DXO(
            data_kind=DataKind.WEIGHTS,
            data={"sums": sums, "counts": counts, "inertia": inertia},
        )
        return out.to_shareable()


class KMeansAggregator(Aggregator):
    """
    Server side:
    accumulates sums and counts from clients, updates centers, returns new centers.
    """
    def __init__(self, k: int, seed: int = 0):
        super().__init__()
        self.k = int(k)
        self.rng = np.random.default_rng(seed)
        self.centers = None
        self._sum = None
        self._count = None
        self._inertia = 0.0
        self._round = 0

    def reset(self):
        self._sum = None
        self._count = None
        self._inertia = 0.0

    def accept(self, shareable, fl_ctx) -> bool:
        dxo = from_shareable(shareable)
        sums = np.asarray(dxo.data["sums"], dtype=np.float64)
        counts = np.asarray(dxo.data["counts"], dtype=np.int64)
        inertia = float(dxo.data.get("inertia", 0.0))

        if self._sum is None:
            self._sum = np.zeros_like(sums, dtype=np.float64)
            self._count = np.zeros_like(counts, dtype=np.int64)

        self._sum += sums
        self._count += counts
        self._inertia += inertia

        # init centers if first time
        if self.centers is None:
            feat_dim = sums.shape[1]
            self.centers = self.rng.normal(size=(self.k, feat_dim)).astype(np.float32)

        return True

    def aggregate(self, fl_ctx):
        if self.centers is None or self._sum is None:
            # nothing received yet, return current (or empty) centers
            out = DXO(data_kind=DataKind.WEIGHTS, data={"centers": self.centers})
            return out.to_shareable()

        new_centers = self.centers.copy()
        for j in range(self.k):
            if self._count[j] > 0:
                new_centers[j] = (self._sum[j] / self._count[j]).astype(np.float32)

        moved = float(np.linalg.norm(new_centers - self.centers))
        self.centers = new_centers
    
        # prints show up in server log (good enough for demo)
        print(f"round {self._round:02d} inertia {self._inertia:.2f} moved {moved:.4f} counts {self._count.tolist()}")
        workspace = os.environ.get("NVFLARE_WORKSPACE") or os.path.abspath("nvflare_workspace_fed_kmeans")
        os.makedirs(workspace, exist_ok=True)

        centers_path = os.path.join(workspace, "final_centers.npy")
        np.save(centers_path, self.centers)
        metrics_path = os.path.join(workspace, "kmeans_metrics.csv")
        
        new_file = not os.path.exists(metrics_path)
        with open(metrics_path, "a", newline="") as f:
            w = csv.writer(f)
            if new_file:
                w.writerow(["round", "inertia", "moved"] + [f"count_{i}" for i in range(self.k)])
            w.writerow([self._round, self._inertia, moved] + self._count.tolist())

        self._round += 1

        self.reset()
        out = DXO(data_kind=DataKind.WEIGHTS, data={"centers": self.centers})
        return out.to_shareable()


class KMeansShareableGenerator(ShareableGenerator):
    def __init__(self, k: int, init_csv: str, seed: int = 0):
        super().__init__()
        self.k = int(k)
        self.init_csv = init_csv
        self.rng = np.random.default_rng(seed)
        self._centers = None

    def _bootstrap_centers(self) -> np.ndarray:
        df = pd.read_csv(self.init_csv)
        drop_cols = {"patient_id", "superpopulation", "n_variants_carried", "primary_gene", "primary_condition", "primary_subtype", "onset_type","base_impact","allelic_ratio","ancestral_modifier","modifier_type","consequence_modifier","interaction_score","predicted_progression", "all_genes","all_subtypes"}
        feat_cols = [c for c in df.columns if c not in drop_cols]
        X = df[feat_cols].to_numpy(dtype=np.float32)

        if X.shape[0] < self.k:
            raise ValueError(f"init_csv has only {X.shape[0]} rows but k={self.k}")

        idx = self.rng.choice(X.shape[0], size=self.k, replace=False)
        return X[idx].copy()

    def learnable_to_shareable(self, learnable: Learnable, fl_ctx) -> Shareable:
        centers = None
        if learnable is not None:
            centers = learnable.get("centers", None)

        if centers is None:
            if self._centers is None:
                self._centers = self._bootstrap_centers()
            centers = self._centers
        else:
            self._centers = np.asarray(centers, dtype=np.float32)

        out = DXO(data_kind=DataKind.WEIGHTS, data={"centers": self._centers})
        return out.to_shareable()

    def shareable_to_learnable(self, shareable: Shareable, fl_ctx) -> Learnable:
        dxo = from_shareable(shareable)
        centers = dxo.data.get("centers", None)

        l = Learnable()
        if centers is not None:
            self._centers = np.asarray(centers, dtype=np.float32)
            l["centers"] = self._centers
        return l
