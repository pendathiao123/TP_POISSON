# Set the output file and the terminal type
set terminal png
set output 'GS_plot.png'

# Set the title, labels, and key
set title "Convergence gauss siedl"
set xlabel "Iteration"
set ylabel "Error"
set key outside

# Use logarithmic scale for both x and y axes
#set logscale x
set logscale y

# Plot data from RESVEC1.dat and RESVEC2.dat
plot 'ALPHA_RESVEC.dat' using 1 with lines title 'gauss-siedl method', \
     #'ALPHA_RESVEC.dat' using 1 with lines title 'richardson_alpha method', \
     #'GS_RESVEC.dat' using 1 with lines title 'gauss-siedl method', \

# Export to png
set output