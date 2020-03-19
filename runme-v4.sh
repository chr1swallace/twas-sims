## output in ~/share/Projects/twas/sims

## version 4 uses eqtl genotype matrices and extends to 1 or 2 cv simulations
for e1 in A AB; do
    for t in A AB AC C; do
	for e2 in same weaker diff t -; do
	    for iN in 500; do
		qR.rb -r simq-v4.R --args NSIM=500 N=$iN NSNP=500 e1=$e1 e2=$e2 t=$t
	    done
	done
    done
done

# 2 x 4 x 5 x 1
