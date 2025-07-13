# ✅ 日本語フォント設定（Colab用）-----------------------------
!apt-get -y install fonts-noto-cjk
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()

# ✅ ライブラリ
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

# ✅ スペクトル特徴読み込み
features = np.load("local_clustered_spectral_features.npy")  # shape = (N, freq)

# ✅ クラスタ数の設定（2〜4で選べるように）
n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
labels = kmeans.fit_predict(features)

# ✅ ラベル保存
np.save("local_clustered_spectral_labels.npy", labels)
print(f"✅ スペクトルクラスタリング完了: クラスタ数 = {n_clusters}")
print(f"✅ ラベル保存ファイル: local_clustered_spectral_labels.npy")
