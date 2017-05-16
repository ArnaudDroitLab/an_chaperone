$F[2]="chr" . $F[2] unless ($F[2] eq "*");

if($F[6] ne "*" && $F[6] ne "=") { 
    $F[6]="chr".$F[6]
}

print join("\t", @F);