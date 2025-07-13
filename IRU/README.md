IRU

Overview

This directory contains supplementary materials for the paper:

Integrated Analysis of Information Structure Preservation and Dimensional Reconstruction: Verification of a Recursive Universe Based on DIRT Theory

Authored by Hiroto Iwasaki, this study investigates the hierarchical robustness and structural significance of the “#1 point” within spectral-preserving data, as proposed in DIRT theory.

⸻

Structure

/Code/

Four executable Python scripts are provided in the order they should be run:

Step	Script	Purpose
1️⃣	Generate_four_files_after_FFT_computation.py	Computes FFT for the three key indicators (info_potential, dc_dt, kappa) and generates four files: fft_info_potential.npy, fft_dc_dt.npy, fft_kappa.npy, fft_freqs.npy (saved to Data/after_fft/).
2️⃣	Construct_spectral_feature_vectors.py	Builds spectral feature vectors from the FFT outputs.
3️⃣	Spectral_features_are_clustered_into_three_groups.py	Clusters the spectral feature vectors into three groups and stores the label array as local_clustered_spectral_labels.npy in /Data/.
4️⃣	p-value.py	Performs a permutation test to evaluate the statistical uniqueness of the “#1 point” cluster assignment.

/Data/
	•	All .npy data files—including the FFT results, spectral cluster labels, and other analysis matrices—are stored here.
	•	For detailed descriptions of each file, please refer to the parent-level README.md in the repository root.

Directory example:

Data/
├── after_fft/
│   ├── fft_dc_dt.npy
│   ├── fft_freqs.npy
│   ├── fft_info_potential.npy
│   └── fft_kappa.npy
├── dc_dt.npy
├── info_potential.npy
├── kappa_list.npy
├── local_clustered_spectral_labels.npy
├── spike_indices.npy
└── y_data.npy


⸻

Key Supplementary Analysis

The final step (p-value.py) assesses the robustness of the “#1 point” classification by permutation testing:

scores = (fft_info + fft_dc_dt + fft_kappa).ravel()
labels = np.load("local_clustered_spectral_labels.npy")  # shape = (N,)
idx_max = np.argmax(scores)
label_obs = labels[idx_max]

n_trials = 5000
match_cnt = sum(np.random.permutation(labels)[idx_max] == label_obs for _ in range(n_trials))
print(f"p-value (クラスタ偏り) = {match_cnt / n_trials:.4f}")

A p-value < 0.01 indicates that the occurrence of the #1 point in its cluster is statistically significant and not due to random chance, supporting its theoretical role in DIRT.

⸻

License

All contents are released under CC0 1.0 Universal.

⸻

Contact
	•	Name: Hiroto Iwasaki
	•	Email: hirotoiwasaki25@gmail.com
	•	GitHub: Hiroto-Iwasaki