#
# Usage:
#    plotStuff plotCircleDeform.cmd -show=<s> -name=<s> -field=[vorticity|streamLines] -res=<i> -vorMin=<f> -vorMax=<f> -numStreamLines=<i> -sol=<list of integers>
#
# Examples: 
# Vorticity:
#  plotStuff plotCircleDeform.cmd -show=iceDeformG8.show -name=iceDeformG8 -field=vorticity -vorMin=-50 -vorMax=50 -sol=100 200
#
# Streamlines:
#  plotStuff plotCircleDeform.cmd -show=iceDeformG8.show -name=iceDeformG8 -field=streamLines -numStreamLines=200 -arrowSize=.1 -res=2048 -sol=100 200
# 
#
$show="iceDeform4.show"; $field="vorticity"; $name="iceDeform"; 
$vorMin=100; $vorMax=-100.; # only use if $vorMin < $vorMax
$numStreamLines=100; $arrowSize=0.25; 
$res=1024; 
@sol= ();  # list of solutions to plot, this must be null for GetOptions to work, defaults are given below
# 
* ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "field=s"=>\$field,"vorMin=f"=>\$vorMin,"vorMax=f"=>\$vorMax,\
            "numStreamLines=i"=>\$numStreamLines,"arrowSize=f"=>\$arrowSize,"res=i"=>\$res,"sol=i{1,}"=>\@sol  );
#
$show
# 
derived types
  vorticity
exit
#
DISPLAY AXES:0 0
solution:$sol[0]
#
contour
  plot:vorticity
  coarsening factor 1 (<0 : adaptive)
  if( $vorMin < $vormMax ){ $cmd = "min max $vorMin $vorMax"; }else{ $cmd="#"; }
  $cmd
  vertical scale factor 0.
  plot contour lines (toggle)
  exit
# zoom: 
set view:0 0.245334 -0.00773827 0 2.79288 1 0 0 0 1 0 0 0 1
# plot streamlines if chosen
if( $field eq "streamLines" ){ $cmd="erase\n stream lines\n  streamline density $numStreamLines\n arrow size $arrowSize\n exit"; }else{ $cmd="#"; }
$cmd
#
pause
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd 
#
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


DISPLAY AXES:0 0
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd 
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=0; $i<$numFreq; $i++ ){ $cmd .= "\n plot:$comp$i\n \$plotName = $name . \"$comp$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd
exit


# -- movie
DISPLAY AXES:0 0
DISPLAY COLOUR BAR:0 0
set view:0 -0.00504559 0.0647128 0 1.45693 1 0 0 0 1 0 0 0 1
movie file name: flowPastDeformingCylinder

save movie files 1

show movie


stream lines
  streamline density $numStreamLines
 exit





derived types
vorticity
exit
#
contour
  plot:vorticity
  coarsening factor 1 (<0 : adaptive)
  min max $vorMin $vorMax
  vertical scale factor 0.
  plot contour lines (toggle)
  exit
# 


# -- movie
DISPLAY AXES:0 0
DISPLAY COLOUR BAR:0 0
set view:0 -0.00504559 0.0647128 0 1.45693 1 0 0 0 1 0 0 0 1
movie file name: twoBeamsInAChannel
save movie files 1
show movie


# -- Hardcopy:
set view:0 0.103601 0.038383 0 1.94067 1 0 0 0 1 0 0 0 1
DISPLAY AXES:0 0
DISPLAY COLOUR BAR:0 0
previous
  hardcopy file name:0 twoBeamsInAChannelVor_t10.ps
  hardcopy save:0