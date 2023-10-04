perl -ane '$h{$F[0]} = $F[1] if (!defined $h{$F[0]} || $h{$F[0]} > $F[1]);                       
           END {foreach (keys %h) {print "$_ $h{$_}\n"}}' HMMER-1-3.order.txt > out_HMMER-1-3.order.txt'.txt
