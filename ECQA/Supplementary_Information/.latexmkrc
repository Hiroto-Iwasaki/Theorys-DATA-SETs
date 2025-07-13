# 使用エンジンを xelatex に設定
$latex = 'xelatex -interaction=nonstopmode -synctex=1 %O %S';
$pdf_mode = 1;

# 出力ディレクトリ
$out_dir = './out';

# デバッグ表示
print "🚀 .latexmkrc is being used: Now using xelatex\n";
