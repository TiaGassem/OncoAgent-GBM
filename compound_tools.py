"""OncoAgent-GBM: compound triage tools (RDKit): structural alerts (PAINS/Brenk), QED,
similarity and clustering. Alerts are WARNINGS for follow-up, never proof of assay interference."""
from __future__ import annotations
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")
_CATS = {"PAINS": FilterCatalogParams.FilterCatalogs.PAINS,
         "Brenk": FilterCatalogParams.FilterCatalogs.BRENK}
_cache: dict = {}


def _catalog(name):
    if name not in _cache:
        p = FilterCatalogParams(); p.AddCatalog(_CATS[name]); _cache[name] = FilterCatalog(p)
    return _cache[name]


def _mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    return m


def structural_alerts(smiles: str) -> list:
    """List of (catalog, alert description) matched."""
    m = _mol(smiles)
    out = []
    for cat in _CATS:
        for e in _catalog(cat).GetMatches(m):
            out.append((cat, e.GetDescription()))
    return out


def qed_score(smiles: str) -> float:
    """Bickerton QED (0-1), a desirability heuristic, not a probability of success."""
    return float(QED.qed(_mol(smiles)))


def murcko_scaffold(smiles: str) -> str:
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(_mol(smiles)))


def _fps(smiles_list, radius=2, nbits=2048):
    return [AllChem.GetMorganFingerprintAsBitVect(_mol(s), radius, nbits) for s in smiles_list]


def tanimoto_matrix(smiles_list):
    fps = _fps(smiles_list)
    return [[DataStructs.TanimotoSimilarity(a, b) for b in fps] for a in fps]


def butina_clusters(smiles_list, cutoff: float = 0.4):
    """Butina clustering on Tanimoto distance (1 - similarity); returns list of index tuples."""
    fps = _fps(smiles_list)
    dists = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1 - s for s in sims)
    return [tuple(c) for c in Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)]
