#**************************************************************************
#
#  Overset grid for two domains to illustrate FSI and CHT with a moving heated annulus
#     OUTER square
#     INNER solid annulus 
#
# *  Usage:
#    ogen noplot io -factor=[1|2...] -interp=[e,i] -order=[2,4,6,8] -rgd=[fixed|var] -outerBC=[d|yPeriodic]
# 
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance* 
#  -outerBC : sets BC on outer boundary
# Examples:
#   ogen -noplot solidAnnulusInASquareGrid -factor=2 -order=2 -interp=e
#   ogen -noplot solidAnnulusInASquareGrid -factor=4 -order=2 -interp=e
# 
# 
#**************************************************************************
$prefix = "solidAnnulusInASquareGrid"; 
$order=2; $factor=1; $interp="i"; $name=""; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$outerBC="dirichlet"; # sets BC on outer boundary
$rgd="var";    $ml=0;
$xa=-1.;  $xb=1.0;  $ya=-1.;  $yb=1.0; 
$deltaRadius0=-1.;  # If set >0 then use this value instead of the default
$numGhost=-1;  # if this value is set, then use this number of ghost points
GetOptions( "order=i"=>\$order,"factor=f"=>\$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "angle=f"=> \$pipeAngle,"startAngle=f"=> \$pipeAngleStart,"interp=s"=> \$interp,\
            "option=s"=> \$option,"rgd=s"=> \$rgd,"name=s"=> \$name,"deltaRadius0=f"=> \$deltaRadius0,\
            "numGhost=i"=>\$numGhost,"outerBC=s"=> \$outerBC,"prefix=s"=> \$prefix,"ml=i"=>\$ml );
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $option ne "" ){ $prefix = $option . "CurvedPipe"; }
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
if( $outerBC eq "yPeriodic" ){ $prefix = $prefix . "Yp"; }
if( $name eq "" ){$name = "$prefix" . "$interp$factor" . $suffix . ".hdf";}
# 
if( $deltaRadius0 < 0 ){ $deltaRadius0=.175 +($order-2)*.1; }
$pi=4.*atan2(1.,1.);
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
# printf("ml=$ml\n");
# pause
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
#
# domain parameters:  
$ds = .05/$factor; # target grid spacing
#
#
$bcInterface=100;  # bc for interfaces
$ishare=100;
# 
create mappings 
#
  annulus 
    mappingName 
      innerAnnulus 
    boundary conditions 
      -1 -1 5 $bcInterface 
    share
       # material interfaces are marked by share>=100
      0 0 0 $ishare
    $outerRadius=.4; 
    $innerRadius=.3;
    inner and outer radii 
      $innerRadius $outerRadius
    lines 
      $nx0= int((2.*$pi*$outerRadius)/$ds+1.5 );
      $nx = intmg($nx0);
      # printf(" nx0=$nx0, nx=$nx, ml=$ml\n");
      $ny0=int( ($outerRadius-$innerRadius)/$ds+2.5 );
      $ny=intmg( $ny0 );
      # printf(" ny0=$ny0, ny=$ny, ml=$ml\n");
      # pause
      $nTheta=$nx;
      $nx $ny 
    exit 
#
  # rectangle 
  #   mappingName
  #     innerSquare
  #   $delta = $ds + ($order-2)*$ds;
  #   $xai=-$innerRadius-$delta;  $xbi=$innerRadius+$delta; 
  #   $yai=-$innerRadius-$delta;  $ybi=$innerRadius+$delta; 
  #   set corners
  #    $xai $xbi $yai $ybi 
  #   lines
  #     $nx=int( ($xbi-$xai)/$ds+1.5 );
  #     $ny=int( ($ybi-$yai)/$ds+1.5 );
  #     $nx $ny
  #   boundary conditions
  #     0 0 0 0
  #   exit 
#
  $deltaRadius=0.1; 
  annulus 
    mappingName 
      outerAnnulus 
    $innerRadius=$outerRadius; 
    $outerRadius=$innerRadius+$deltaRadius;
    inner and outer radii 
      $innerRadius $outerRadius
    lines 
      $nx=$nTheta; 
      $ny=intmg( ($deltaRadius)/$ds+2.5 );
      $nx $ny
    boundary conditions 
      -1 -1 $bcInterface 0 
    share
 # material interfaces are marked by share>=100
      0 0 $ishare 0   
    exit 
#
  rectangle 
    mappingName
      outerSquare
    set corners
     $xa $xb $ya $yb 
    lines
      $nx=intmg( ($xb-$xa)/$ds+1.5 );
      $ny=intmg( ($yb-$ya)/$ds+1.5 );
      $nx $ny
    boundary conditions
      $cmd="1 2 3 4 5 6";
      if( $outerBC eq "yPeriodic" ){ $cmd="1 2 -1 -1"; }
      $cmd
    exit 
  exit this menu 
#
generate an overlapping grid 
  outerSquare
  outerAnnulus 
  # innerSquare
  innerAnnulus 
  done 
#
  change parameters 
 # define the domains -- these will behave like independent overlapping grids
    specify a domain
 # domain name:
      outerDomain 
 # grids in the domain:
      outerSquare
      outerAnnulus
      done
    specify a domain
 # domain name:
      innerDomain 
 # grids in the domain:
      # innerSquare
      innerAnnulus
      done
    ghost points 
      all 
      $ng $ng $ng $ng $ng $ng
    order of accuracy
     $orderOfAccuracy
    interpolation type
      $interpType
    exit 
    # open graphics
    compute overlap
# pause
  exit
save a grid
  $name
  solidAnnulusInABox
exit
