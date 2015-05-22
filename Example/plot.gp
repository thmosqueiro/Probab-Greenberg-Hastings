set border lw 4
set tics scale 1.

set log xy

set terminal postscript eps enhanced "Helvetica" 20
set output '| epstopdf --filter --outfile=plot_PDFs.pdf'

set xrange[0.3:100000]

plot 'pdf_experiment_s08.dat' w points pt 5 ps 1.5 lc rgb "#00AA55",\
     'pdf_experiment_s10.dat' w points pt 7 ps 1.5 lc rgb "#2255DD"



set terminal postscript eps enhanced "Helvetica" 20
set output '| epstopdf --filter --outfile=plot_CDFS.pdf'

set xrange[1:100000]

plot 'cdf_experiment_s08.dat' w points pt 5 ps 1.5 lc rgb "#00AA55",\
     'cdf_experiment_s10.dat' w points pt 7 ps 1.5 lc rgb "#2255DD"
