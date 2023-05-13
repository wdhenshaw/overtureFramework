eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to exttract the performance results from Ogmg output
#  usage: 
#         getPerfResults
# 

@fileNames = @ARGV;

foreach $fileName ( @fileNames )  # process all files
{
  
  open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  while( <FILE> )
  {
    $line = $_;
    if( $line =~ /\% PerfInfo/ )
    {
      # printf("%s",$line);
      # Clean up the line : 
      $line =~ s/ \% PerfInfo//;
      $line =~ s/\.order2.*\.hdf//;
      $line =~ s/\.order4.*\.hdf//;
      $line =~ s/(Ogmg[^&]*)/\1    /;                    # add some spaces to solver name to line up with others
      $line =~ s/PETSc.*hypre[^&]*/AMG               /;  # shorten name to AMG
      $line =~ s/PETSc.*bi-conjugate gradient stabilized.*(ILU\(.\))[^&]*/Bi-CG-Stab \1 /; # shorten name to Bi-CG-Stab
      $line =~ s/PETSc.*generalized minimal residual.*(ILU\(.\))[^&]*/GMRES \1      /;

      printf("%s",$line);

    }

  }

  # close(OUTFILE);
  close(FILE);


}