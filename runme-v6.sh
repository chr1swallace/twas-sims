## output in ~/share/Projects/twas/sims

## version 4 uses eqtl genotype matrices and extends to 1 or 2 cv simulations
iN=500
for t in A ; do # AB
    for e1 in A B ; do # - BA AB AC 

        ## A-AB equiv to A-AC
        if [[ "$e1" = "A"  && "$t" = "AC" ]]; then
            continue
        fi
        ## A-AB, A-AC, AB-B equiv to AB-A
        if [[ "$t" = "AB"  && ("$e1" = "B" || "$e1" = "A" || "$e1" = "-") ]]; then
            continue
        fi
        if [[ "$t" = "A"  && "$e1" = "AC" ]]; then
            continue
        fi


	for e2 in same diff t -; do
            ## must be some eqtl
            if [[ "$e1" = "-" && ("$e2" = "-" || "$e2" = "same") ]]; then
                continue
            fi
            ## if e1==t, then same and t are equivalent
            if [[ "$e1" = "$t"  && "$e2" = "t" ]]; then
                continue
            fi
            echo "$t :: $e1 :: $e2"
	    qR.rb -r simq-v6.R --args rho=0.25 NSIM=400 N=$iN NSNP=500 e1=$e1 e2=$e2 t=$t
	    qR.rb -r simq-v6.R --args rho=0 NSIM=400 N=$iN NSNP=500 e1=$e1 e2=$e2 t=$t
	done
    done
done

# 4 x 4 x 4 x 1

## should be 20000
## block size of 20 means 1000 runs
echo simulations
ruby -e 'print Dir.glob("/home/cew54/share/Projects/twas/sims/simv6*").count'
echo results
ruby -e 'print Dir.glob("/home/cew54/share/Projects/twas/sims/results/*").count'

## running time for one dataset = ~ 7 minutes with 2 runs of fuser 
date
Rscript analyse-v6.R --args BATCH.SIZE=1 --taskid=0
date

qR.rb -y 0-800 -r -c 1 analyse-v6.R --args BATCH.SIZE=10

## put all the runs together
qR.rb -r collate-v6.R

## plot
Rscript plot.R
