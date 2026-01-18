README（提案：簡潔版・日英併記）

F2DV — √2 離散真空からのフレーバー（F2DV）

F2DV — Flavor from a √2-Discrete Vacuum (F2DV)

TL;DR / 要点

【日本語】 本リポジトリは、論文「Flavor from a √2-Discrete Vacuum (F2DV)」の全結果（数値計算・図・表・主要チェック）を、固定乱数シードの下で端から端まで再現するための最小セットです。
【English】 This repository is the minimal package to reproduce all paper results end-to-end (numerics, figures, tables, and key checks) under fixed random seeds.

⸻

1. 最優先：Colab ノートブック（これだけで全再現）

【日本語】 まずは Code/F2DV_paper_full_repro_end_to_end.ipynb を Colab で開き、**「セッションを再起動してすべて実行」**でワンクリック再現してください（日本語フォントの !apt 問題は解消済み）。
【English】 Open Code/F2DV_paper_full_repro_end_to_end.ipynb in Colab and run “Runtime → Restart and run all” for one-click reproduction (the Japanese font !apt issue is fixed).

⸻

2. 生成物は Data/ に集約

【日本語】 Data/ には、ノートブック／スクリプトが生成した図（png）と要約CSVを格納します（論文に貼る最終成果物置き場）。
【English】 Data/ stores all generated figures (png) and summary CSVs (final artifacts used in the paper).

例（あなたのスクショの構成と一致）
	•	gaussian_mc_NO_dist_epsilon_l_13.png, gaussian_mc_NO_epsilon_l_13_vs_sum_mnu.png
	•	gaussian_mc_IO_dist_epsilon_l_13.png, gaussian_mc_IO_epsilon_l_13_vs_sum_mnu.png
	•	mbb_phase_scan_summary_NO.csv, mbb_phase_scan_summary_IO.csv
	•	pmns_epsilon13_gaussian_mc_*_with_mbeta.csv など

⸻

3. 再現性：乱数設定と代表値の定義（論文と同一）

3.1 Gaussian MC（ε13^(ℓ), Σmν, mβ）
	•	【日本語】 N=200000, seed=12345 を固定しています。
	•	【English】 We fix N=200000 and seed=12345.

3.2 0νββ 位相走査（mββ）
	•	【日本語】 N_phase=300000, seed_phase=20260111（※採用値）を固定しています。
	•	【English】 We fix N_phase=300000 and seed_phase=20260111 (adopted value).

3.3 代表質量の定義（重要）
	•	【日本語】 m_i^{\rm rep} は 各オクタントの MC サンプルの mean（平均）で定義します（実装通り）。
	•	【English】 We define m_i^{\rm rep} as the mean of the MC samples for each octant (as implemented).

⸻

4. Code/ の役割（ノートブック優先、.py は補助）

【日本語】 Code/ の .py は、再現実行の補助・切り出し用です。優先順位は ipynb が最上位です（今後 .py が増えても「再現は ipynb」で統一）。
【English】 The .py files are helper/standalone scripts. The notebook is the single source of truth for reproduction (even if more .py files are added later).

主なスクリプト例：
	•	Gaussian_NO_and_IO.py：NO/IO の Gaussian MC と図・CSV生成
	•	delta-angle_check_for_m_bb.py：δ固定の妥当性チェック（同一の Majorana 位相サンプルを固定して比較）
	•	NO-main_seesaw_validation_Casas–Ibarra.py：Casas–Ibarra / 直交性検証（論文の検証パート）

⸻

5. 主要式（最小） / Core equations (minimal)

5.1 混合則（GAV型）

\log_2(\sin\theta_{ij})\simeq-(p\Delta K_{ij}+r\Delta C_{ij})
	•	意味（各項）：\theta_{ij} 混合角、\Delta K_{ij},\Delta C_{ij} 離散差分、p,r フィット係数
	•	役割：(12),(23) から p,r を決め、(13) の基準値に接続する

5.2 残差保存（レプトン側）

\epsilon_{13}^{(\ell)}=\sqrt2-\log_2\!\left(\frac{s_{13,\rm base}}{s_{13,\rm obs}}\right),\qquad
x_{13,\rm eff}^{(\ell)}=\sqrt2-\epsilon_{13}^{(\ell)}
	•	意味（各項）：s_{13,\rm base} は(12,23)由来の基準、s_{13,\rm obs} は観測入力、\epsilon_{13}^{(\ell)} は残差、x_{13,\rm eff}^{(\ell)} は実効変数
	•	役割：\sqrt2 原理からのズレを \epsilon_{13}^{(\ell)} として分離し、以降の予言（EFT/シーソー）に渡す

5.3 観測量

\sum m_\nu=m_1+m_2+m_3,\quad
m_\beta^2=\sum_i |U_{ei}|^2 m_i^2,\quad
m_{\beta\beta}=\left|\sum_i U_{ei}^2 m_i e^{i\alpha_i}\right|
	•	意味（各項）：m_i 質量固有値、|U_{ei}| PMNS e行、\alpha_i Majorana位相
	•	役割：\sum m_\nu, m_\beta（位相非依存）と m_{\beta\beta}（位相依存の帯）を、論文の比較パネルとして出力する

⸻

6. 注意（定義のズレ防止）

【日本語】 \epsilon_{13}^{(\ell)} は「log2 比そのもの」ではなく、論文の x=\sqrt2-\epsilon 形式に整合するよう 残差として定義しています。
【English】 \epsilon_{13}^{(\ell)} is defined as a residual to match the paper convention x=\sqrt2-\epsilon, not merely as the raw log2-ratio itself.

⸻

Citation / 引用

【日本語】 本リポジトリの内容を用いる場合は論文（プレプリント）を引用してください。
【English】 If you use materials from this repository, please cite the paper/preprint.

⸻