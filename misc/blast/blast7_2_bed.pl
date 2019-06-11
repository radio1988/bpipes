use strict; use warnings;

while(<>){
    chomp;
    next if /^#/;
    my@F=split"\t", $_;
    if (int($F[8]) > int($F[9])) {
        print join"\t", $F[1],$F[9],$F[8],$F[0],$F[11],"-\n"
    }else{
        print join"\t", $F[1],$F[8],$F[9],$F[0],$F[11],"+\n"
    }
}
