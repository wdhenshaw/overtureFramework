# 
# Plot results from heated bodies
# *new way:
#   plotStuff plotSolution -show=movingHeatedAnnulusG4.show -name=movingHeatedAnnulusG4 -sol=51 61 71
#
#
#   plotStuff plotSolution -show=fourHeatedDisks.show -name=fourHeatedDiskTempt0p7 -solution=15 -Tmin=0 -Tmax=10
#   plotStuff plotSolution -show=16HeatedDisks.show -name=16HeatedDiskTempt0p7 -solution=15 -Tmin=0 -Tmax=10
#
#   plotStuff plotSolution -show=16HeatedDisksG4.show -name=16HotDisksINS -solution=1 -Tmin=0 -Tmax=10  [movie]
#   plotStuff plotSolution -show=16HeatedDisksG4.show -name=16HeatedDisk -Tmin=0 -Tmax=10 -numToPlot=6
# 
# RPI:
#   plotStuff plotSolution -show=rpiCHT.show -name=rpiCHTt0p5 -solution=11 -Tmin=0 -Tmax=10 -domain1=fluidDomain -domain2=solidDomain
#   plotStuff plotSolution -show=rpiCHT.show -name=rpiCHTt0p75 -solution=16 -Tmin=0 -Tmax=10 -domain1=fluidDomain -domain2=solidDomain
#   plotStuff plotSolution -show=rpiCHT.show -name=rpiCHTt1p0 -solution=21 -Tmin=0 -Tmax=10 -domain1=fluidDomain -domain2=solidDomain
# 
$show="fourHeatedDisks.show"; $name="fourHeatedDisks"; 
$Tmin=0; $Tmax=-1; $field="T"; $res=1024;
$domain1="outerDomain"; $domain2="innerDomain"; 
@sol= ();  # list of solutions to plot, this must be null for GetOptions to work, defaults are given below
# ----------------------------- get command line arguments --------------------------------------- 
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"solution=i"=>\$solution,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax ,\
            "numToPlot=i"=>\$numToPlot,"domain1=s"=>\$domain1,"domain2=s"=>\$domain2,"sol=i{1,}"=>\@sol, "field=s"=>\$field ); 
# 
$show 
# 
# -- set displacement components for Godunov 
frame series:$domain1
# derived types 
# speed 
# exit 
# 
# previous 
# plot:stressNorm 
DISPLAY AXES:0 0
solution:$sol[0]
#
plot:T
contour 
  plot contour lines (toggle)
  vertical scale factor 0.0
  if( $Tmax>$Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
# plot streamlines if chosen
if( $field eq "streamLines" ){ $cmd="erase\n stream lines\n  streamline density $numStreamLines\n arrow size $arrowSize\n exit"; }else{ $cmd="#"; }
$cmd
#
frame series:$domain2
contour
  plot contour lines (toggle)
  vertical scale factor 0.0
  if( $Tmax>$Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
#
pause
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd 
#
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd



pause
# --- movie:
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
bigger:0
DISPLAY COLOUR BAR:0 0
movie file name: $name
save movie files 1
show movie

# -- save one hardcopy
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . $label . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0




# Foreground colour:0 white
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0



# --plot solutions at different times:
#
@sol=(  31,   51,   71,  101,  151, 201);
@time=("0p3","0p5","0p7","1p0","1p5","2p0"); 
$cmd="#"; 
for( $i=0; $i<$numToPlot; $i++ ){ $plotName = $name . "Tt$time[$i].ps"; $cmd .="\n solution: $sol[$i]\n hardcopy file name:0 $plotName\n hardcopy save:0"; }
printf("cmd=$cmd\n");
$cmd



