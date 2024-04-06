open FILE,  "$ARGV[0]" or die "I can not open the input FILE\n";

while ($line=<FILE>){
    if ($line=~/\>lcl\|.+peg\./){
        $line=~/(\>lcl\|)(.+)(.peg.\d+)( unnamed protein product)/g;
        $rast_id="$2";
        $peg_number="$3";
        $peg_number=~s/\./_/g;
        
        open RASTIDS,  "RAST.ids" or die "I can not open the RASTIDS\n";
        while ($ids=<RASTIDS>){
            if ($ids=~/$rast_id/){
                $ids=~/(.+)\t(.+)\t(.+)/;
                $specie_gi="$3";
                $specie_gi=~/(.+)\s+(.+)$/;
                $specie="$1";
                $specie=~s/ sp. /_/g;
                $specie=~s/ \= /_/g;
        $specie=~s/\=/_/g;
                $specie=~s/\-/_/g;
                $specie=~s/ cf. /_/g;
                $specie=~s/ subsp. /_/g;
                $specie=~s/ str. /_/g;
                $specie=~s/\./_/g;
                $specie=~s/\s+/_/g;
                $specie=~s/\//_/g;
                $specie=~s/\)//g;
                $specie=~s/\(//g;
                $rast_id=~s/6666666.//g;
                print ">$rast_id$peg_number\_$specie\n";
            }
            
        }


    }
    else {
      print $line;
    }
}
