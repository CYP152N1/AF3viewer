# AF3viewer.py 
# AlphaFold Server (AF3用)データ解析ツール

iPTM, PAE, pLDDT viewer via local html for Alphafold server (AF3)

このツールは、AlphaFold Serverによって生成されたJSONおよびCIFデータを処理し、タンパク質モデルの評価や可視化を行うための解析ツールです。主にpLDDTスコアの評価、アミノ酸配列の解析、コンタクト確率とPAE（Pairwise Aligned Error）の可視化、Cα構造の3Dプロットなどを行います。
このプログラムとREADMEはほとんどCHatGPT4oさんに作ってもらいました。自己責任でご利用ください。

## 機能詳細

### 1. pLDDTスコアのプロット
- 各モデルの残基ごとにpLDDTスコアをプロットし、タンパク質構造の局所的な信頼性を評価します。
- プロットはPNG形式で保存され、対応するCSVファイルも出力されます。

### 2. アミノ酸頻度のプロット
- 各タンパク質チェーンにおけるアミノ酸の出現頻度を計算し、ヒストグラムとして可視化します。
- 結果はPNG形式で保存されます。

### 3. コンタクト確率とPAEのプロット
- モデル間のコンタクト確率マップとPAEマップを作成し、タンパク質構造間の相互作用を評価します。
- 結果はPNGおよびCSV形式で保存されます。

### 4. Cα構造の3Dプロット
- CIFファイルから抽出したCα原子の座標を使用して、3Dプロットを作成します。
- 各チェーンを異なる色で表示し、構造の可視化を行います。
- プロットはPNG形式で保存されます。

### 5. チェーンペア行列のプロット
- チェーン間のiPTMスコアおよび最小PAE値を行列としてプロットします。
- 行列データはCSV形式でも出力されます。

### 6. モデルの整列
- 指定したチェーンを基準にして、複数のモデルを整列させます。
- 整列後の構造はCIF形式で保存されます。

### 7. GIFの生成
- 各プロットを連続画像としてGIF形式で保存し、構造の変化や特徴を動的に表示します。

### 8. HTMLレポートの生成
- 全ての解析結果をまとめたHTMLレポートを生成し、ブラウザで閲覧可能な形で提供します。

## 使用方法

### 前提条件

以下のパッケージが必要です。

- Python 3.x
- Biopython
- Matplotlib
- Pillow
- NumPy

これらのパッケージは以下のコマンドでインストールできます。

```bash
pip install biopython matplotlib pillow numpy
```

### 実行手順

コマンドラインから以下のように実行します。

```bash
python3 AF3viewer.py -i <入力フォルダ> [-r <ランキング基準>] [-a <整列チェーン>]
```

- `-i`: 入力フォルダの名前（必須）alphafoldサーバーからダウンロードしたzipファイルを解凍してできたフォルダを指定してください。フォルダ名を変えると動きません。まzip解凍後のフォルダの上の真階層でしか動きません。
- `-r`: ランキングの基準（`fraction_disordered`, `iptm`, `ptm`, `ranking_score`から選択、デフォルトは`ranking_score`）
- `-a`: 整列に使用するチェーン（例: `"A"` または `"AB"`、オプション、選択しない場合はすべてのchainでアライメント）

### 例

```bash
python3 AF3viewer -i example_data -r ranking_score -a A
```

## 出力ファイル

- **PNG画像**: 各種プロット（pLDDTスコア、アミノ酸頻度、コンタクト確率、PAE、Cα構造、チェーンペア行列など）
- **CSVファイル**: 各種データのCSVファイル（pLDDTスコア、コンタクト確率、PAE、チェーンペア行列など）
- **GIFファイル**: 各プロットの連続画像GIF
- **HTMLレポート**: 全ての解析結果をまとめたHTML形式のレポート

## ライセンス

このプロジェクトはMITライセンスの下で提供されています。

## 貢献

バグの報告や新機能のリクエストは、Issueを作成してください。また、プルリクエストも歓迎します。プロジェクトの改善にご協力いただけることを楽しみにしています。


## 出力されるhtmlのイメージ

![image](https://github.com/user-attachments/assets/8c4f5302-4e10-4ffb-a6a5-98e77f65bcc3)

![image](https://github.com/user-attachments/assets/b880350b-fa2f-4c01-b231-b87b936feb4d)

https://github.com/user-attachments/assets/ce3479ed-293e-4b80-9f85-05427f391e3c

https://github.com/user-attachments/assets/2e75a488-0d38-41ec-ab9e-857032b51760

![image](https://github.com/user-attachments/assets/fa156379-be3e-4d07-91cd-74618ec12eed)

![image](https://github.com/user-attachments/assets/bc7e108d-add1-484a-9aa0-2d72436738f0)

![image](https://github.com/user-attachments/assets/4c9143a9-a9e2-4bd6-baa2-77d025bf5d20)
