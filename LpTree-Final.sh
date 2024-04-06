cat Order_LpTree.txt | xargs -n 1 -I % awk '{if ( $1 == "%") print $0}' final.results100.txt > HMMER_HitsLp.txt
