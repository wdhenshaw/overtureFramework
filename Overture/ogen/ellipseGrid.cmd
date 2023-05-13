#
#   Grid for the interior of an ellipse
#
# usage: ogen [-noplot] ellipseGrid -radX=<f> -radY=<f> -factor=<num> -order=[2/4/6/8] -interp=[e/i] -numGhost=<i> ...
#
# Example:
#    ogen -noplot ellipseGrid -interp=e -factor=2
#    ogen -noplot ellipseGrid -interp=e -order=4 -factor=2
#    ogen -noplot ellipseGrid -interp=e -order=4 -factor=4
#    ogen -noplot ellipseGrid -interp=e -order=4 -numGhost=3 -factor=2
# 
$prefix="ellipseGrid";
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
# Fixed radial distance
$nDistInterface=.2;
$ml=0; 
#
$cx=0.; $cy=0;          # ellipse center
$radX=1.2; $radY=1.0; 
$theta=0.;    # angle in degrees
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"radX=f"=> \$radX,"radY=f"=> \$radY,"cx=f"=> \$cx,"cy=f"=> \$cy,\
            "interp=s"=> \$interp,"name=s"=> \$name,"numGhost=i"=> \$numGhost,"prefix=s"=> \$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
$pi=4.*atan2(1.,1.);
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
$theta=$theta*$pi/180.;
$numberOfVolumeSmooths=0; 
#
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
$minRad=min($radX,$radY);
#
#
# --- make the ellipse using a spline
#    NOTE: go in a counter clockwise direction
#
$scmd="#";
$ns=50*$factor; $arcLength=0.;
for( $i=1; $i<=$ns; $i++ ){ $s=2.*$pi*($i-1.)/($ns-1.); \
  $x=$radX*cos($s)*cos($theta)-$radY*sin($s)*sin($theta)+$cx;           \
  $y=$radY*sin($s)*cos($theta)+$radX*cos($s)*sin($theta)+$cy;           \
   if( $i > 1 ){ $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 );} $x0=$x; $y0=$y; \
   $scmd .= "\n $x $y"; }
# 
create mappings
#
spline
  #
  enter spline points
    $ns
    $scmd
  lines
    $ns
    periodicity
      2
  mappingName
    splineBoundary
  # pause
 exit
# 
# -- Make a hyperbolic grid --
#
  # $nDistInterface=.2 + ($order-2)*$ds; 
  # $nr = intmg( 5 + $order/2 );
  hyperbolic
    backward
    $nDistExtraFactor=1.2; # make grid a bit finer in the normal direction
    $nrm = intmg(($nDistInterface*$nDistExtraFactor)/$ds+.5);
    # $nDistInterface=($nr-3)*$ds;
    # $nrm=$nr-1; 
    distance to march $nDistInterface
    lines to march $nrm
    $nThetaInterface = intmg($arcLength/$ds+.5);
    points on initial curve $nThetaInterface
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    # pause
    generate
    # use fourth order interpolant to define the mapping: *wdh* May 8, 2018
    fourth order
    # evaluate as nurbs 1
    boundary conditions
      -1 -1 1 0 0 0
    share 
       0 0 1 0 0 0
    name ellipseGrid
    # pause
  exit 
#
#  --- inner background grid
#
 rectangle
  set corners
    # FIX FOR ROTATED GRID: 
    $xas=$cx-$radX; $xbs=$cx+$radX; $yas=$cy-$radY; $ybs=$cy+$radY; 
    $xas $xbs $yas $ybs
    $refineFactor=1.;
  lines
    $nx = int( ($xbs-$xas)/$ds*$refineFactor + .5 ); 
    $ny = int( ($ybs-$yas)/$ds*$refineFactor + .5 ); 
    $nx $ny
  boundary conditions
    0 0 0 0 
  share 
    0 0 0 0 
  mappingName
    backGround
  exit
#
exit this menu
#
generate an overlapping grid
  backGround
  ellipseGrid
  done
  #
  change parameters
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
  # open graphics
  compute overlap
exit
#
# save an overlapping grid
save a grid (compressed)
$name
ellipseGrid
exit



#
# --- make the ellipse using a spline
#    NOTE: go in a counter clockwise direction
#
$scmd="#";
$ns=50*$factor; $arcLength=0.;
for( $i=1; $i<=$ns; $i++ ){ $s=2.*$pi*($i-1.)/($ns-1.); \
  $x=$radX*cos($s)*cos($theta)-$radY*sin($s)*sin($theta)+$cx;           \
  $y=$radY*sin($s)*cos($theta)+$radX*cos($s)*sin($theta)+$cy;           \
   if( $i > 1 ){ $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 );} $x0=$x; $y0=$y; \
   $scmd .= "\n $x $y"; }
# 

#
#  ----- Given arrays of points xc,yc, build curves and grids -----
#
# Update lists 
$count +=1; $shareInterface += 1; $bcInterface+=1;
$fluidGridName="fluidInterface$count";
$solidGridName="solidInterface$count";
$solidBackGroundName="solidBackGround$count";
#
$solidCoreName="solidCore$count";
## $solidGridNames .= "\n $solidBackGroundName\n $solidCoreName\n $solidGridName"; 
# -- remove the core grid 
$solidGridNames .= "\n $solidBackGroundName\n$solidGridName"; 
#
$fluidGridNames .= "\n $fluidGridName";
## $solidDomains .= "\n specify a domain\n solidDomain$count\n $solidBackGroundName\n $solidCoreName\n $solidGridName\n done";
$solidDomains .= "\n specify a domain\n solidDomain$count\n $solidBackGroundName\n $solidGridName\n done";
#
if( $degree eq "" ) { $degree=3; } # degree of Nurbs 
if( $numberOfVolumeSmooths eq "" ){ $numberOfVolumeSmooths=20; }  #
# rotation angle (degress)
if( $angle eq "" ){ $angle=0.; }  
if( $pi eq "" ){ $pi=4.*atan2(1.,1.); }  
# Core offset from center
if( $xCore eq "" ){ $xCore=0.; }
if( $yCore eq "" ){ $yCore=0.; }
if( $coreRadius eq "" ){ $coreRadius=.05; }
#
$xcMin=1e10; $xcMax=-1e10; $ycMin=1e10; $ycMax=-1e10;  # for inner background grid 
$cmd="#";  $ns=0;  $arcLength=0.; $x0=$xc[0]; $y0=$yc[0];
for( $ic=0; $ic<$nc-1; $ic++ ){ $numpt=$numPerSegment; if( $ic == $nc-2 ){ $numpt=$numPerSegment+1;} \
for( $i=0; $i<$numpt; $i++ ){ $s=($i)/($numPerSegment); $ns=$ns+1;  \
   $xs=(1.-$s)*$xc[$ic]+$s*$xc[$ic+1]; \
   $ys=(1.-$s)*$yc[$ic]+$s*$yc[$ic+1];   \
   $ct=cos($angle*$pi/180.); $st=sin($angle*$pi/180.);  \
   $x=$cx + $ct*($xs-$cx)-$st*($ys-$cy);   \
   $y=$cy + $st*($xs-$cx)+$ct*($ys-$cy);   \
   $xcMin = min($xcMin,$x); $xcMax = max($xcMax,$x);  \
   $ycMin = min($ycMin,$y); $ycMax = max($ycMax,$y);  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); printf(" ns=$ns: x0=$x0 x=$x y0=$y0 y=$y arcLength=$arcLength\n");  $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; } }\
  $knots="#"; for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; } \
#
nurbs (curve)
  periodicity
   2
  enter control points
    $degree
    $ns
    $knots
    $cmd 
 #
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
  $curveName="boundaryInitial$count";
  $curveName
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    $curveName
  mappingName
   $curveName2="boundaryCurve$count";
   $curveName2
exit 
# 
# -- Make EXTERIOR hyperbolic grid --
#
  $nr = intmg( 6 + 3*($order-2)/2 );
  hyperbolic
    forward
    $nDist=($nr-3.25)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    # pause
    boundary conditions
      -1 -1 $bcInterface 0 0 0
    share 
       0 0 $shareInterface 0 0 0
    name $fluidGridName
  exit
# 
# -- Make INTERIOR hyperbolic grid --
#
   $dsi=$ds*.75;  # target grid spacing for interior
# 
  $nr = intmg( 6 + 3*($order-2)/2 );
  hyperbolic
    backward
    $nDist=($nr-3.25)*$dsi;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths 20 
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    # pause
    boundary conditions
      -1 -1 $bcInterface 0 0 0
    share 
       0 0 $shareInterface 0 0 0
    name $solidGridName
  exit
#
#   Inner background -- make a bit finer than target
#
   rectangle
   set corners
    $xai=$xcMin; $xbi=$xcMax;
    $yai=$ycMin; $ybi=$ycMax;
    $xai $xbi $yai $ybi
   lines
    $nxi = int( ($xbi-$xai)/$dsi +1.5 ); 
    $nyi = int( ($ybi-$yai)/$dsi +1.5 ); 
    $nxi $nyi
    boundary conditions
      0 0 0 0 
    mappingName
      $solidBackGroundName
  exit 

