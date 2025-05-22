#日本語フォント設定（Colab用）-----------------------------
!apt-get -y install fonts-noto-cjk
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()

#ライブラリ
import numpy as np

#スペクトル特徴の読み込み（各: shape = (N,) の1次元配列かも）
fft_info = np.load("fft_info_potential.npy")
fft_dc_dt = np.load("fft_dc_dt.npy")
fft_kappa = np.load("fft_kappa.npy")

#次元が1Dなら明示的に2Dへ変換（reshape）
if fft_info.ndim == 1:
    fft_info = fft_info.reshape(-1, 1)
if fft_dc_dt.ndim == 1:
    fft_dc_dt = fft_dc_dt.reshape(-1, 1)
if fft_kappa.ndim == 1:
    fft_kappa = fft_kappa.reshape(-1, 1)

#結合して統合ベクトルを作成（axis=1で特徴ベクトル化）
features = np.concatenate([fft_info, fft_dc_dt, fft_kappa], axis=1)

#ローカルに保存
np.save("local_clustered_spectral_features.npy", features)

print("スペクトル特徴ベクトル（local_clustered_spectral_features.npy）を生成・保存しました。")
print("shape:", features.shape)
