# IRU

## Overview

This directory contains supplementary materials for the paper:

**Integrated Analysis of Information Structure Preservation and Dimensional Reconstruction: Verification of a Recursive Universe Based on DIRT Theory**

Authored by **Hiroto Iwasaki**, this study investigates the hierarchical robustness and structural significance of the "#1 point" within spectral-preserving data, as proposed in DIRT theory.

---

## Structure

- `/Code/` – Contains `p-value.py`, a permutation test script evaluating the statistical uniqueness of the #1 point's cluster assignment.
- `/Data/` – Contains `local_clustered_spectral_labels.npy`, a NumPy array representing spectral cluster labels.

---

## Key Supplementary Analysis

The file `p-value.py` performs a **permutation test** to assess the robustness of the "#1 point" classification:

```python
scores = (fft_info + fft_dc_dt + fft_kappa).ravel()
labels = np.load("local_clustered_spectral_labels.npy")   # shape = (N,)
idx_max = np.argmax(scores)
label_obs = labels[idx_max]

n_trials = 5000
match_cnt = sum(np.random.permutation(labels)[idx_max] == label_obs for _ in range(n_trials))
p_value = match_cnt / n_trials

print(f"p-value (クラスタ偏り) = {p_value:.4f}")

