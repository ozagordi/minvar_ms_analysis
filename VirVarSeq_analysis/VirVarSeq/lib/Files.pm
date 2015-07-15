package Files;

sub openFile{
   my ($in) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my @tmp = ();
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      push(@tmp, $line);
   }
   close(IN);
   return @tmp;
}

sub openFileColumn{
   my ($in, $colnr) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my @tmp = ();
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      my @values = split("\t",$line);
      push(@tmp, $values[$colnr-1]);
   }
   close(IN);
   return @tmp;
}


sub openFileAssociative{
   my ($in) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my %tmp;
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      my @tmparray = split("\t",$line);
      my $key = shift(@tmparray);
      $tmp{$key} = join("\t",@tmparray);
      # print $key."\t".$value."\n";
   }
   close(IN);
   return %tmp;
}

sub openFileAssociativeWhiteSpace{
   my ($in) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my %tmp;
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      my @tmparray = split("\s+",$line);
      my $key = shift(@tmparray);
      $tmp{$key} = join("\t",@tmparray);
      # print $key."\t".$value."\n";
   }
   close(IN);
   return %tmp;
}

sub openFileAssocColumn{
   my ($in, $col1, $col2) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my %tmp;
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      my @tmparray = split("\t",$line);
      my $key = $tmparray[$col1-1];
      my $value = $tmparray[$col2-1];
      $tmp{$key} = $value;
   }
   close(IN);
   return %tmp;
}

sub openFileToHash{
   my ($in) = @_;
   
   open(IN, $in) or die "File $in does not exists\n";
   my %tmp;
   foreach my $line(<IN>){
      $line =~ s/\r|\n//sg;
      $tmp{$line} =1;
   }
   close(IN);
   return %tmp;
}

1;
