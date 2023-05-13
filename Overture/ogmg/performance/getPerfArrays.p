eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to extra performance data and put into Matlab arrays
#  usage: 
#         getPerfArrays.p perf.dat > perfData.m
# 

@fileNames = @ARGV;

foreach $fileName ( @fileNames )  # process all files
{
  
  open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  @gridNames = ();
  $gridNum=-1; 
  @solverNames = ();
  $cpuTotal = 6; # column location of the total cpu
  $cpuSolve = 8; # solve time  
  $storage  = 9; # storage, reals/pt
  printf("cpuTotal=1; %% location in array with cpu total\n");
  printf("cpuSolve=2; %% location in array with cpu total\n");
  printf("storage =3; %% location in array with storage in reals/pt\n");
  while( <FILE> )
  {
    $line = $_;
    chop($line);
    if( $line =~ /&/ )
    {
      @cols = split('&', $line);

      # printf("%s",$line);
      if( 1==0 )
      {
        for( $i=0; $i< @cols; $i++ )
        {
          printf("cols[$i]=$cols[$i]\n");
        }
      }
      $grid=$cols[0];
      $grid =~ s/^[ ]*//;  # remove leading blanks
      $grid =~ s/[ ]*$//;  # ... and trailing    
      if( $gridNum<0 || $grid ne $gridNames[$gridNum] )
      {
        $gridNum=$gridNum+1;
        $gridNames[$gridNum]=$grid;
        printf("$grid = %d; %% grid name enumerator\n",$gridNum+1);
      }
      $order = $cols[1]; 
      $solver = $cols[2];
      $solver =~ s/^[ ]*//;  # remove leading blanks
      $solver =~ s/[ ]*$//;
      $solverNum=-1;
      for( $i=0; $i<@solverNames; $i++ )
      {
        if( $solver eq $solverNames[$i] )
        {
          $solverNum=$i; break;
        }
      }
      if( $solverNum==-1 )
      { $solverNames[@solverNames]=$solver; $solverNum=@solverNames-1; 
        printf("solverName{%d}=\"$solver\";\n",$solverNum+1); 
      }

      $storageValue = $cols[$storage];
      $storageValue =~ s/\\\\//; # remove ending //
      $sn=$solverNum+1;
      printf("data(cpuTotal,$order,$sn,$gridNames[$gridNum])=$cols[$cpuTotal]; %% data(:,order,solver,grid)= value\n");
      printf("data(cpuSolve,$order,$sn,$gridNames[$gridNum])=$cols[$cpuSolve]; \n");
      printf("data(storage ,$order,$sn,$gridNames[$gridNum])=$storageValue; \n");


    }

  }

  # close(OUTFILE);
  close(FILE);


}