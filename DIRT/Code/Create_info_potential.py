import numpy as np

# ===== 既存ファイルの読み込み =====
x_data = np.load("x_data.npy")[:9999]          # Δt（ゼロ点間隔）
primes = np.load("primes.npy")[:9999]          # 素数列

# ===== prime_density の生成（1 / log(p)）=====
prime_density = 1 / np.log(primes)

# ===== 情報ポテンシャル = Δt × prime_density =====
info_potential = x_data * prime_density

# ===== 保存 =====
np.save("info_potential.npy", info_potential)

print("✅ info_potential.npy および prime_density.npy を保存しました。")