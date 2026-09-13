import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.metrics import silhouette_score, davies_bouldin_score

RDLogger.DisableLog("rdApp.*")

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

df = pd.read_csv("data/prob_veriseti.csv", sep=";", encoding="utf-8")


def safe_mol_from_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        return mol
    except Exception:
        return None


def mol_features(smiles):
    mol = safe_mol_from_smiles(smiles)
    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None


def get_murcko_scaffold(smiles):
    mol = safe_mol_from_smiles(smiles)
    if mol is None:
        return None
    return MurckoScaffold.GetScaffoldForMol(mol)


valid_rows = []
X_list = []
y_em_list = []
y_abs_list = []

for idx, s in enumerate(df["SMILES"]):
    f = mol_features(s)
    if f is not None:
        valid_rows.append(idx)
        X_list.append(f)
        y_em_list.append(df.loc[idx, "Lambda_em_nm"])
        y_abs_list.append(df.loc[idx, "Lambda_abs_nm"])

X = np.array(X_list)
y_em = np.array(y_em_list)
y_abs = np.array(y_abs_list)
df_valid = df.loc[valid_rows].reset_index(drop=True)

print(f"X shape: {X.shape}, y_em shape: {y_em.shape}, y_abs shape: {y_abs.shape}")

n = len(df_valid)
tanimoto_matrix = np.ones((n, n))
fps = [morgan_gen.GetFingerprint(safe_mol_from_smiles(s)) for s in df_valid["SMILES"]]

for i in range(n):
    for j in range(i + 1, n):
        sim = TanimotoSimilarity(fps[i], fps[j])
        tanimoto_matrix[i, j] = sim
        tanimoto_matrix[j, i] = sim

off_diag = tanimoto_matrix[np.triu_indices(n, k=1)]
print(f"Tanimoto similarity range: {off_diag.min():.2f} - {off_diag.max():.2f}")
print(f"Pairs above 0.95 similarity: {int((off_diag > 0.95).sum())}")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2, random_state=42)
pca_result = pca.fit_transform(X_scaled)
df_valid["pca1"] = pca_result[:, 0]
df_valid["pca2"] = pca_result[:, 1]

print(f"PC1 explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%")
print(f"PC2 explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%")

best_k = 4
kmeans = KMeans(n_clusters=best_k, init="k-means++", n_init=50, random_state=42)
cluster_labels = kmeans.fit_predict(X_scaled)
df_valid["Cluster"] = cluster_labels

sil_score = silhouette_score(X_scaled, cluster_labels)
db_score = davies_bouldin_score(X_scaled, cluster_labels)
print(f"Silhouette score: {sil_score:.3f}")
print(f"Davies-Bouldin score: {db_score:.3f}")

for c in sorted(df_valid["Cluster"].unique()):
    subset = df_valid[df_valid["Cluster"] == c]
    print(f"Cluster {c}: n={len(subset)}, emission range={subset['Lambda_em_nm'].min():.0f}-{subset['Lambda_em_nm'].max():.0f} nm")

stokes_shift = y_em - y_abs
print(f"Stokes shift range: {stokes_shift.min():.1f} - {stokes_shift.max():.1f} nm")
print(f"Positive Stokes shift for all molecules: {bool((stokes_shift > 0).all())}")

X_train_em, X_test_em, y_train_em, y_test_em = train_test_split(
    X, y_em, test_size=0.2, random_state=42
)
em_model = RandomForestRegressor(n_estimators=200, random_state=42)
em_model.fit(X_train_em, y_train_em)
em_pred_test = em_model.predict(X_test_em)
em_r2 = r2_score(y_test_em, em_pred_test)
em_mae = mean_absolute_error(y_test_em, em_pred_test)
em_cv_scores = cross_val_score(em_model, X, y_em, cv=5, scoring="r2")

print(f"Emission R2 (test): {em_r2:.2f}")
print(f"Emission MAE (test): {em_mae:.1f} nm")
print(f"Emission R2 (5-fold CV): {em_cv_scores.mean():.2f} +/- {em_cv_scores.std():.2f}")

X_train_abs, X_test_abs, y_train_abs, y_test_abs = train_test_split(
    X, y_abs, test_size=0.2, random_state=42
)
abs_model = RandomForestRegressor(n_estimators=200, random_state=42)
abs_model.fit(X_train_abs, y_train_abs)
abs_pred_test = abs_model.predict(X_test_abs)
abs_r2 = r2_score(y_test_abs, abs_pred_test)
abs_mae = mean_absolute_error(y_test_abs, abs_pred_test)
abs_cv_scores = cross_val_score(abs_model, X, y_abs, cv=5, scoring="r2")

print(f"Absorption R2 (test): {abs_r2:.2f}")
print(f"Absorption MAE (test): {abs_mae:.1f} nm")
print(f"Absorption R2 (5-fold CV): {abs_cv_scores.mean():.2f} +/- {abs_cv_scores.std():.2f}")

em_importances = em_model.feature_importances_
abs_importances = abs_model.feature_importances_
top_em_bits = np.argsort(em_importances)[::-1][:10]
top_abs_bits = np.argsort(abs_importances)[::-1][:10]

print(f"Top 10 fingerprint bits for emission prediction: {top_em_bits.tolist()}")
print(f"Top 10 fingerprint bits for absorption prediction: {top_abs_bits.tolist()}")


def find_most_similar(smiles, dataset_smiles):
    mol1 = safe_mol_from_smiles(smiles)
    if mol1 is None:
        return None, 0
    fp1 = morgan_gen.GetFingerprint(mol1)
    max_sim = -1
    best_match = None
    for idx, s in enumerate(dataset_smiles):
        mol2 = safe_mol_from_smiles(s)
        if mol2 is None:
            continue
        fp2 = morgan_gen.GetFingerprint(mol2)
        sim = TanimotoSimilarity(fp1, fp2)
        if sim > max_sim:
            max_sim = sim
            best_match = idx
    return best_match, max_sim


def predict_new_molecule(smiles):
    features = mol_features(smiles)
    if features is None:
        print("Geçersiz SMILES")
        return

    features = features.reshape(1, -1)
    em_pred = em_model.predict(features)[0]
    abs_pred = abs_model.predict(features)[0]

    idx, sim = find_most_similar(smiles, df["SMILES"])
    if idx is None:
        print("Benzer molekül bulunamadı")
        return

    most_similar = df.iloc[idx]

    print("Benzer molekülün Prob IDsi:", most_similar["Prob_ID"])
    print(f"En benzer molekül: {most_similar['SMILES']} (Similarity: {sim:.2f})")
    print("Tahmini analit:", most_similar["Hedef_Analit"])
    print("Tahmini fotofiziksel mekanizma:", most_similar["Fotofiziksel_Mekanizma"])
    print(f"Tahmini Emisyon (nm): {em_pred}")
    print(f"Tahmini Absorbans (nm): {abs_pred}")
    print("Emisyon (gerçek):", most_similar["Lambda_em_nm"])
    print("Absorbans (gerçek):", most_similar["Lambda_abs_nm"])
    print("Biyogörüntüleme sonucu:", most_similar["Ana_Biyogoruntuleme_Sonuc"])
    print("Görüntüleme tipi:", most_similar["Goruntuleme_Tipi"])
    print("Subselüler lokalizasyon:", most_similar["Subseluler_Lokalizasyon"])
    print("Yapısal özellik notları:", most_similar["Yapı_ozellik_Notları"])


predict_new_molecule("c1ccc2c(c1)nc3ccccc23")
predict_new_molecule("c1ccc2c(c1)nc3ccccc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccncc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccccn3c2")
predict_new_molecule("c1ccc2c(c1)n(c3ccccc23)C")
predict_new_molecule("c1ccc2c(c1)nc3ccc(Cl)cc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccc(F)cc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccc(Br)cc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccc(C#N)cc3c2")
predict_new_molecule("c1ccc2c(c1)nc3ccc(C(C)C)cc3c2")
