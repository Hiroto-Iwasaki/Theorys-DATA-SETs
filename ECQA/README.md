# Observational Evidence for a Cosmic Quadrupole Anisotropy (ECQA) - Supplementary Materials

## Overview / 概要

This directory contains supplementary materials, including data and analysis code, for the research paper "Observational Evidence for a Cosmic Quadrupole Anisotropy: A Test of the Information Lock Coefficient in DIRT Cosmology."

This paper challenges the cosmological principle by testing a framework derived from Dimensional Infinite Regression Theory (DIRT), which predicts an intrinsic anisotropy in the universe's expansion. Through a joint analysis of Type Ia supernovae (Pantheon+) and Baryon Acoustic Oscillation (BAO) data, we report a definitive, multi-sigma detection of a cosmic quadrupole anisotropy. This result offers a new physical perspective on long-standing cosmological anomalies and tensions.

本ディレクトリは、研究論文「DIRT宇宙論における非等方性の観測的証拠と理論的帰結」の補足資料（データおよび解析コード）を含みます。

本研究は、宇宙の膨張が内因的に非等方的であると予測する次元無限回帰理論（DIRT）の枠組みを、Ia型超新星（Pantheon+）とバリオン音響振動（BAO）データの統合解析によって検証し、宇宙の四重極子非等方性の決定的証拠を報告するものです。この結果は、長年の宇宙論的アノマリーやテンション問題に対し、新たな物理的視点を提供します。

---

## Code / コード内容

The `/Code` directory contains the main Jupyter Notebook (`ECQA_analysis.ipynb`) used for the analysis presented in the paper and its supplementary information.

* **Description**: The notebook covers all stages of the analysis, including:
    * Robustness evaluation of the \#1 point structure using t-SNE and FFT.
    * Formulation and calculation of the DIRT Quadrupole model's predictions (anisotropic Hubble parameter, distances, BAO scales).
    * The final joint MCMC analysis using Pantheon+ and BAO data.
* **Dependencies**: The code relies on standard Python libraries such as `numpy`, `pandas`, `matplotlib`, `scipy`, `astropy`, `emcee`, and `corner`.
* **How to Run**: Execute the cells in the notebook sequentially. Note that the MCMC analysis cells are computationally intensive and may take a significant amount of time to run.

`/Code`ディレクトリには、論文および補足資料で提示された解析を実行するための主要なJupyter Notebook (`ECQA_analysis.ipynb`) が含まれています。ノートブックには、頑健性評価（t-SNE, FFT）、DIRT Quadrupoleモデルの理論予測計算、そして最終的なSNe+BAO統合MCMC解析の全工程が含まれています。実行には `numpy`, `pandas`, `astropy`, `emcee`, `corner` などのライブラリが必要です。

---

## Data / データ概要

The `/Data` directory contains the input data required to run the analysis scripts and the key output files used to generate the figures in the paper.

* **Input Data**:
    * `x_data.npy`, `y_data.npy`, `dc_dt.npy`, `kappa_list.npy`, etc.: The foundational time-series data derived from DIRT simulations (from Paper C).
    * `Pantheon+SH0ES.dat`: The public data file for the Pantheon+ Type Ia supernovae compilation.
* **Key Output Files**:
    * `Mcmc_not_test.png`: The main result of this work, the corner plot from the joint SNe+BAO MCMC analysis.
    * `S1_*.png`, `S2_*.png`: Figures used in the Supplementary Information document.

`/Data`ディレクトリには、解析の実行に必要な入力データ（DIRTシミュレーションから得られた基礎データやPantheon+カタログ）と、論文中の図の作成に使用された主要な出力ファイルが格納されています。

---

## Citation / 引用情報

If you use any material from this study, please cite the main paper:

研究のいずれかの資料を使用する場合は、主要な論文を引用してください：

> Iwasaki, H. (2025). "Observational Evidence for a Cosmic Quadrupole Anisotropy: A Test of the Information Lock Coefficient in DIRT Cosmology." *(Journal/Preprint details, e.g., DOI: 10.13140/RG.2.2.12237.14565)*

---