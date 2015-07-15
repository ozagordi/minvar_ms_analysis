package Log;

sub screen {
   my ($script, $txt) = @_;
   my $date = `date +%F" "%T`;
   chomp($date);
   print $date . "\t" . $script . "\t" . $txt . "\n";
   return();
}

1;
