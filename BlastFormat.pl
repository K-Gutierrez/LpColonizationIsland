#making a blast database of all .faa files

@fastas=qx/ls|grep .faa/;
foreach (@fastas){
    chomp $_;
    system "makeblastdb -in $_ -dbtype prot -parse_seqids -logfile makeblastdb.log";
    $cont++;
    print "$cont $_ has been added to the blast database\n";

}
