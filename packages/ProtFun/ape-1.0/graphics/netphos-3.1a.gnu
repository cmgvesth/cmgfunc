set terminal pbm small color
set output "| cat ";
set title "NetPhos 3.1a: predicted phosphorylation sites in _NAME_"
set size 1.0,0.6
set xrange [0:_LEN_]
set yrange [0:1.3]
set tics out
set ytics 0,1
set xlabel "Sequence position" 0,0
set ylabel "Phosphorylation potential" 0,0
plot 'S.dat' t"Serine" w i 3, \
     'T.dat' t"Threonine" w i 8, \
     'Y.dat' t"Tyrosine" w i 1, \
     0.5 t"Threshold"w l 0
