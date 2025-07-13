#事前にロードしておく
scores = (fft_info + fft_dc_dt + fft_kappa).ravel()
labels = np.load("local_clustered_spectral_labels.npy")   # shape = (N,)

#実データの観測値
idx_max   = np.argmax(scores)
label_obs = labels[idx_max]          # 極大スコアが属するクラスタ
print(f"実データ: 最大スコア点のクラスタ = {label_obs}")

#permutation test
n_trials   = 5000
match_cnt  = 0
for _ in range(n_trials):
    shuffled_labels = np.random.permutation(labels)
    if shuffled_labels[idx_max] == label_obs:
        match_cnt += 1

p_value = match_cnt / n_trials       # “同じクラスタになる” 片側確率
print(f"p-value (クラスタ偏り) = {p_value:.4f}")