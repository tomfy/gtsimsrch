width = ARG1
height = ARG2
input_file = ARG3
nsnps = ARG4
xmax = ARG5
prsst = ARG6
trmnl = ARG7
png_filename = ARG8
ttl = ARG9


if(prsst == 1){
	 set terminal ARG7 persist size width,height
}else{
	set terminal ARG7 size width,height
}

f(x) = -log(x)/log(10)
bonf_threshold = -1*f(20*nsnps)

# input file format is
# snpid  chrom  position_on_chrom  cume_position  p-value  causal
# causal is 0 for non-causal snps, 1 for causal snps.

set xlabel 'marker position'
set ylabel '-log(p)'
set grid
set title ttl

plot [0:xmax] \
     input_file u 4:(f($5)):(($6==1)? 1 : 1):((int($2)%2==0)? 3 : 4) pt 6  ps variable lc variable t'', \
      '' u 4:(f($5)):(($6==1)? 1 : 0):(($6==1)? 8 : 0) pt 7  ps variable lc variable t'', \
      bonf_threshold

set terminal png
set out png_filename

plot [0:xmax] \
     input_file u 4:(f($5)):(($6==1)? 1 : 1):((int($2)%2==0)? 3 : 4) pt 6  ps variable lc variable t'', \
      '' u 4:(f($5)):(($6==1)? 1 : 0):(($6==1)? 8 : 0) pt 7  ps variable lc variable t'', \
      bonf_threshold

# pause -1