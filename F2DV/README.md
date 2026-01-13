⸻

GAV — Geometric Architecture of Vacuum

Dual Integration, Dimensional Transmutation, and the Origin of Flavor
（真空の幾何学的建築：二重の双対性・次元転換・フレーバーの起源）
This represents the geometric architecture of vacuum: double duality, dimensional transmutation, and the origin of flavor.

⸻

Overview / 概要

このディレクトリは、研究論文
「The Geometric Architecture of Vacuum: Dual Integration, Dimensional Transmutation, and the Origin of Flavor」
（真空の幾何学的建築：二重の双対性、次元転換、そしてフレーバーの起源）
の解析コード・幾何学図・バイナリデータをまとめたものです。
This directory contains the analysis scripts, geometric diagrams, and binary data associated with the research paper “The Geometric Architecture of Vacuum: Dual Integration, Dimensional Transmutation, and the Origin of Flavor.”

⸻

Research Context / 研究背景

本論文では、標準模型におけるフレーバー階層が、
真空の離散的幾何学構造（Platonic dual pairs / Lie group dimensions）
によって決定されるという
GAV 質量生成スキーム（GAV Mass-Generation Scheme）
を提示します。
This work proposes that the flavor hierarchy of the Standard Model is determined by the discrete geometric structure of the vacuum—specifically Platonic dual pairs and the dimensional structure of Lie groups—through the GAV Mass-Generation Scheme.

⸻

1. 質量 = 幾何学的指数 K による指数法則

1. Mass follows an exponential law determined by the geometric index K.

m \propto (\sqrt{2})^{K}

質量の倍化は振幅の √2 倍に対応するという場の幾何学的構造から導かれる。
Doubling of mass corresponds to a √2 scaling of field amplitude, derived from the geometric structure of quantum fields.

⸻

2. K は真空の対称性（プラトン立体 + リー群）の開放履歴

2. K encodes the unfolding history of vacuum symmetry (Platonic solids + Lie groups).
	•	レプトン：SU(5) の K ≈ 24 で停止
Leptons terminate at the SU(5) boundary with K \approx 24.
	•	クォーク：G₂ を経由し K ≈ 32 (SO(32) 的深層) まで到達
Quarks pass through G₂ and reach K \approx 32, corresponding to an SO(32)-like depth.
	•	これは カシミール選択則 によって数学的に決定される
This structure is determined by a Casimir-selection rule.

⸻

3. トップクォークは「二重の双対性」を持つ唯一の粒子

3. The top quark is the only particle exhibiting “double duality.”
	•	増分：K ≈ 14 = Cube(6) + Octa(8)（G₂）
Increment: K \approx 14 = 6 + 8 from cube and octahedral duality (G₂).
	•	累積：K ≈ 32 = Dodeca(12) + Icosa(20)（SO(32) 的）
Cumulative: K \approx 32 = 12 + 20 from dodeca–icosa duality (SO(32)-like).

→ “Double Dual Integration”
→ This structure is termed “Double Dual Integration.”

⸻

4. ニュートリノは正四面体の自己双対性 → K = −48

4. Neutrinos inherit the self-duality of the tetrahedron → K = -48.
	•	幾何学的シーソー機構
A geometric seesaw mechanism emerges.
	•	予言値：
m_{\nu 3} \approx 0.030\ \mathrm{eV}
Predicted mass:
m_{\nu 3} \approx 0.030\ \mathrm{eV}

⸻

Directory Structure / ディレクトリ構造

GAV/
 ├── Code/
 │    ├── G2_root_system_diagram.py
 │    ├── gav_k_index_identification.py
 │    ├── gav_unified_fermion_spectrum.py
 │    └── top_quark_mass_projection.py
 │
 └── Data/
       ├── G2_root_system_diagram.png
       └── GAV_Cumulative_Geometric_Indices.json

上記は GAV のコード群とデータファイルの構成を示します。
The above represents the directory structure containing all GAV-related code and data files.

⸻

Code / コード内容と役割

● G2_root_system_diagram.py

役割
	•	G₂のロングルート／ショートルート（3:1 異方性）の図示
Visualizes the long/short root structure of G₂ (3:1 anisotropy).
	•	フレーバー階層（アップ型 vs ダウン型）に対応させる
Connects the G₂ structure to the up/down flavor hierarchy.

論文対応箇所：Sec.3, Sec.6（Fig.1）
Corresponds to Sections 3 and 6 (Fig.1) of the paper.

⸻

● gav_k_index_identification.py

役割
	•	観測質量データから指数 K を逆算
Back-computes the geometric index K from observed mass data.
	•	K ≈ 8, 11, 14, 15, 24, 32 の整数構造を抽出
Extracts the integer-like structure of K (≈8, 11, 14, 15, 24, 32).

論文対応箇所：Sec.2.3（Table 1）
Corresponds to Section 2.3 (Table 1).

⸻

● gav_unified_fermion_spectrum.py

役割
	•	SU(5) → G₂ → SO(32) への階層的質量生成を可視化
Visualizes hierarchical mass generation from SU(5) to G₂ to SO(32).
	•	レプトン／クォーク分岐の カシミール選択則 を反映
Implements the Casimir-selection rule for lepton/quark branching.

論文対応箇所：Sec.4
Corresponds to Section 4.

⸻

● top_quark_mass_projection.py

役割
	•	トップクォークの「増分 14」「累積 32」の幾何学的投影の再現
Reproduces the geometric projection of the top quark (increment 14, cumulative 32).
	•	K に基づくトップ質量の幾何学的スケーリング確認
Verifies geometric mass scaling using the index K.

論文対応箇所：Sec.3
Corresponds to Section 3.

⸻

Data / データ内容

● G2_root_system_diagram.png
	•	G₂ root system の図
Diagram of the G₂ root system.
	•	論文 Fig.1 の生成に使用
Used to produce Fig.1 of the paper.

⸻

● GAV_Cumulative_Geometric_Indices (JSON/CSV)
	•	各粒子の K（increment/cumulative）の一覧
List of geometric indices K (increment/cumulative) for all fermions.
	•	SU(5) → G₂ → SO(32) への階層構造の計算結果
Computed hierarchical structure: SU(5) → G₂ → SO(32).

⸻

How to Use / 使用方法

1. 環境構築 / Environment Setup

pip install numpy scipy sympy matplotlib pandas

Installs all required dependencies.

⸻

2. 指数 K の抽出 / Extract the K indices

python gav_k_index_identification.py

Computes geometric indices from observed fermion masses.

⸻

3. 幾何学スペクトルの生成 / Generate the geometric spectrum

python gav_unified_fermion_spectrum.py

Creates the SU(5) → G₂ → SO(32) cumulative-K visualization.

⸻

4. G₂ root system の可視化 / Visualize the G₂ root system

python G2_root_system_diagram.py

Outputs the G₂ geometric diagram.

⸻

5. トップクォーク幾何学予言の検証 / Validate top-quark geometric predictions

python top_quark_mass_projection.py

Checks the K=14 and K=32 dual-projection consistency.

⸻

Citation / 引用

このディレクトリの内容を使用する場合、以下を引用してください：
If you use any material from this directory, please cite the following:

Iwasaki, H. (2025).
“The Geometric Architecture of Vacuum: Dual Integration, Dimensional Transmutation, and the Origin of Flavor.”
Preprint.

⸻