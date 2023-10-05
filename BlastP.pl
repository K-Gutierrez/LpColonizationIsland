open OUT,  ">>entry_batch";

#Blasting a query against all .faa dbs
if (defined $ARGV[0]){
if (defined $ARGV[1]){
if (defined $ARGV[2]){
if (defined $ARGV[3]){
$query= "$ARGV[0]";
@fastas=qx/ls|grep .faa|grep -v .faa./;
foreach (@fastas){
    chomp $_;
    $dbname= "$_";
    @blastouts=qx/blastp -query $query -db $_ -evalue $ARGV[2] -outfmt 6/;
        foreach (@blastouts){
          @columns=split(/\s+/,$_);
            if ($columns[11]>=$ARGV[3]){
            print "$dbname\t$_";   
            $entry="$columns[1]";
            print OUT  "$entry\n";
            
            }
        }
}



print "This search was performed with -evalue $ARGV[2] and score cutoff of $ARGV[3]\n";


}}}}

open FASTAS,  ">$ARGV[1]";

system "sort entry_batch|uniq>entry_batch_uniq";

foreach(@fastas){
    chomp $_;
    $dbname= "$_";
system "blastdbcmd -db $dbname -dbtype prot -entry_batch entry_batch_uniq -out temp";

system "cat temp >>$ARGV[1]";

}



system "rm entry_batch* temp";

print "#####USAGE######\n";
print "blastin.pl [query] [fastahits filename] [-evalue 0.000001 is recomended] [score cutoff i.e. 100]\n";
    exit;


#to get the number of hits per genome type perl blastin.pl query.fasta |cut -f 1|sort|uniq -c\n";
