#
#   Objects in a big domain
#
#
# usage: ogen [noplot] objectsGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot objectsGrid -interp=e -order=2 -factor=2 
#     ogen -noplot objectsGrid -interp=e -order=2 -factor=4 
#     ogen -noplot objectsGrid -interp=e -order=4 -factor=4 
#     ogen -noplot objectsGrid -interp=e -order=4 -factor=8
#
#     ogen -noplot objectsGrid -prefix=multiObjectsGrid -addBackGround=0 -interp=e -order=2 -factor=2 
#     ogen -noplot objectsGrid -prefix=multiObjectsGrid -addBackGround=0 -interp=e -order=2 -factor=4 
#     ogen -noplot objectsGrid -prefix=multiObjectsGrid -addBackGround=0 -interp=e -order=4 -factor=4 
#
# *wdh* Oct 3, 2024 -- change the boundary conditions to be distinct
#
$prefix="objectsGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-3.; $xb=5.5; $ya=-3.5; $yb=3.5; $ml=0; 
$addBackGround=1; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix, "addBackGround=i"=>\$addBackGround );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.05/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
create mappings
#
#  Coarse background 
#
  rectangle
    mappingName
      backGround
   set corners
    $dsc=2.*$ds; # coarse spacing 
    $width=4;
    $xac=$xa-$width; $xbc=$xb+$width; $yac=$ya-$width; $ybc=$yb+$width; 
    $xac $xbc $yac $ybc
   lines
    $nx = int( ($xbc-$xac)/$dsc +1.5 ); 
    $ny = int( ($ybc-$yac)/$dsc +1.5 ); 
    $nx $ny
  exit
  # Inner refinement 
  rectangle
    mappingName
      innerBackGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds +1.5 ); 
    $ny = int( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      # if( $addBackGround eq 1 ){ $bcCmd="0 0 0 0"; }else{ $bcCmd="1 1 1 1"; }
      if( $addBackGround eq 1 ){ $bcCmd="0 0 0 0"; }else{ $bcCmd="1 2 3 4"; }
      # 0 0 0 0
      $bcCmd
  exit
#
$gridNames="#"; 
#
# --- build I-beam grid 
#
$count=0;
# For IBeam: 
$centerHeight=.75; $centerWidth=.25;
$edgeHeight=.25;   $edgeWidth=1.;
include $ENV{CG}/mx/runs/scattering/iBeam.h
# include iBeam.h
$mapName="iBeamGridBase";
# 
$xr=.0; $yr=0; # center for rotation
$xShift=-1; $yShift=0; $angle=0; 
$count=$count+1; $gridName="iBeam$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
$xShift=-1; $yShift=1.5; $angle=90; 
$count=$count+1; $gridName="iBeam$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
$xShift=-1; $yShift=-1.5; $angle=90; 
$count=$count+1; $gridName="iBeam$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
# --- build base split-ring
#
$count=0;
include $ENV{CG}/mx/runs/scattering/splitRing.h
$mapName="splitRingGridBase";
# 
$xr=.0; $yr=0; # center for rotation
$xShift=3.5; $yShift=0; $angle=0; 
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
$xShift=3.5; $yShift=1.5; $angle=90; 
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
$xShift=3.5; $yShift=-1.5; $angle=-90; 
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
#
#  -- build the baseTriangle
include $ENV{CG}/mx/runs/scattering/triangle.h
#
# build shifted/rotated/scaled versions
$count=0;
$mapName="triangleBase";
# 
$xr=.5; $yr=0; # center for rotation
$xShift=0; $yShift=0; $angle=45; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
# 
$xShift=0; $yShift=1.5; $angle=20; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
# 
$xShift=0; $yShift=-1.5; $angle=-30; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
# 
$xShift=1.5; $yShift=0; $angle=-45; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
# 
$xShift=1.5; $yShift=1.5; $angle=80; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
# 
$xShift=1.5; $yShift=-1.5; $angle=0; 
$count=$count+1; $gridName="triangle$count"; $gridNames .= "\n" . $gridName;
include $ENV{CG}/mx/runs/scattering/transform.h
#
exit 
#
generate an overlapping grid
  if( $addBackGround eq 1 ){ $cmd="backGround\n innerBackGround"; }else{ $cmd="innerBackGround"; }
  $cmd 
  # backGround
  # innerBackGround
  $gridNames
  done
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
   # pause
  # open graphics
  compute overlap
exit
#
save an overlapping grid
$name
triangle
exit
