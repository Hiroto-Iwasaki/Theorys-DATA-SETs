GAV — Geometric Architecture of Vacuum

Dual Integration, Dimensional Transmutation, and the Origin of Flavor

（真空の幾何学的建築：二重の双対性・次元転換・フレーバーの起源）

⸻

Overview / 概要

このディレクトリは、研究論文
「The Geometric Architecture of Vacuum: Dual Integration, Dimensional Transmutation, and the Origin of Flavor」
（真空の幾何学的建築：二重の双対性、次元転換、そしてフレーバーの起源）
￼
の解析コード・幾何学図・バイナリデータをまとめたものです。

⸻

Research Context / 研究背景

本論文では、標準模型におけるフレーバー階層が、
真空の離散的幾何学構造（Platonic dual pairs / Lie group dimensions）
によって決定されるという
GAV 質量生成スキーム（GAV Mass-Generation Scheme）
を提示します。

中心となる主張：

1. 質量 = 幾何学的指数 K による指数法則

m \propto (\sqrt{2})^{K}

質量の倍化は振幅の √2 倍に対応するという場の幾何学的構造から導かれる。

⸻

2. K は真空の対称性（プラトン立体 + リー群）の開放履歴
	•	レプトン：SU(5) の K ≈ 24 で停止
	•	クォーク：G₂ を経由し K ≈ 32 (SO(32) 的深層) まで到達
	•	これは カシミール選択則 によって数学的に決定される

⸻

3. トップクォークは「二重の双対性」を持つ唯一の粒子
	•	増分：K ≈ 14 = Cube(6) + Octa(8)（G₂）
	•	累積：K ≈ 32 = Dodeca(12) + Icosa(20)（SO(32) 的）
→ “Double Dual Integration”

⸻

4. ニュートリノは正四面体の自己双対性 → K = −48
	•	幾何学的シーソー機構
	•	予言値：
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
       └── GAV_Cumulative_Geometric_Indices.json（あるいは類似ファイル）


⸻

Code / コード内容と役割

● G2_root_system_diagram.py

役割
	•	G₂のロングルート／ショートルート（3:1 異方性）の図示
	•	フレーバー階層（アップ型 vs ダウン型）に対応させる
論文対応箇所：Sec.3, Sec.6（Fig.1） ￼

⸻

● gav_k_index_identification.py

役割
	•	観測質量データから指数 K を逆算
	•	K ≈ 8, 11, 14, 15, 24, 32 の整数構造を抽出
論文対応箇所：Sec.2.3（Table 1） ￼

⸻

● gav_unified_fermion_spectrum.py

役割
	•	SU(5) → G₂ → SO(32) への階層的質量生成を可視化
	•	レプトン／クォーク分岐の カシミール選択則 を反映した計算
論文対応箇所：Sec.4（Casimir operator formalism） ￼

⸻

● top_quark_mass_projection.py

役割
	•	トップクォークの「増分 14」「累積 32」の幾何学的投影の再現
	•	K に基づくトップ質量の幾何学的スケーリング確認
論文対応箇所：Sec.3（Double Dual Integration） ￼

⸻

Data / データ内容

● G2_root_system_diagram.png
	•	G₂ root system の図
	•	論文 Fig.1 の生成に使用  ￼

⸻

● GAV_Cumulative_Geometric_Indices (JSON/CSV)
	•	各粒子の K（increment/cumulative）の一覧
	•	SU(5) → G₂ → SO(32) への階層構造の計算結果

⸻

How to Use / 使用方法
	1.	環境構築

pip install numpy scipy sympy matplotlib pandas


	2.	指数 K の抽出

python gav_k_index_identification.py


	3.	幾何学スペクトルの生成

python gav_unified_fermion_spectrum.py


	4.	G₂ root system の可視化

python G2_root_system_diagram.py


	5.	トップクォーク幾何学予言の検証

python top_quark_mass_projection.py



⸻

Citation / 引用

このディレクトリの内容を使用する場合、以下を引用してください：

Iwasaki, H. (2025).
“The Geometric Architecture of Vacuum: Dual Integration, Dimensional Transmutation, and the Origin of Flavor.”
Preprint.  ￼

⸻