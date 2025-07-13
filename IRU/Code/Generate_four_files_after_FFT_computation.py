#日本語フォント設定（Colab用）-----------------------------
!apt-get -y install fonts-noto-cjk
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()

#必要ライブラリ
import numpy as np
from scipy.fft import fft, fftfreq

#データ読み込み
info_potential = np.load("info_potential.npy")
dc_dt = np.load("dc_dt.npy")
kappa = np.load("kappa_list.npy")
spike_indices = np.load("spike_indices.npy")
y_data = np.load("y_data.npy")

#スパイク点に対応するデータ抽出
ip_spike = info_potential[spike_indices]
dc_spike = dc_dt[spike_indices]
kp_spike = kappa[spike_indices]

#補間：等間隔でないため線形補間してFFT（共通のx軸を用意）
n_interp = 300
y_interp = np.linspace(y_data[spike_indices].min(), y_data[spike_indices].max(), n_interp)

from scipy.interpolate import interp1d
f_ip = interp1d(y_data[spike_indices], ip_spike, kind='linear', fill_value="extrapolate")
f_dc = interp1d(y_data[spike_indices], dc_spike, kind='linear', fill_value="extrapolate")
f_kp = interp1d(y_data[spike_indices], kp_spike, kind='linear', fill_value="extrapolate")

ip_uniform = f_ip(y_interp)
dc_uniform = f_dc(y_interp)
kp_uniform = f_kp(y_interp)

#FFT計算
fft_ip = np.abs(fft(ip_uniform))[:n_interp // 2]
fft_dc = np.abs(fft(dc_uniform))[:n_interp // 2]
fft_kp = np.abs(fft(kp_uniform))[:n_interp // 2]
freqs = fftfreq(n_interp, d=(y_interp[1] - y_interp[0]))[:n_interp // 2]

#保存
np.save("fft_info_potential.npy", fft_ip)
np.save("fft_dc_dt.npy", fft_dc)
np.save("fft_kappa.npy", fft_kp)
np.save("fft_freqs.npy", freqs)

print("4ファイルをローカルに保存しました")
