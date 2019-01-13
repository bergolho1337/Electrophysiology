set title "Noble Model"
set grid
set xlabel "t (s)"
set terminal png
set output "noble.png"
plot "data.dat" using 1:2 title "v" w l