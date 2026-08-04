#
# cgmp: Solve the CHT problem for the RPI solids with INS 
#     : flow around a heated cylinder (2d or 3d)
#
echo to terminal 0
# 
# Usage:
#    cgmp [-noplot] rpiCHT.cmd -g=<name> -method=[ins|cns] -nu=<> -mu=<> -kappa=<num> -tf=<tFinal> -tp=<tPlot> ...
#           -solver=<yale/best> -ktcFluid=<> -ktcFluid=<> -tz=[poly/trig/none] -bg=<backGroundGrid> ...
#           -nc=<num> degreex=<num> -degreet=<num> -ts=[fe|be|im|pc] -coupled=[0|1]  -mixedInterface=[0|1] ...
#           -go=[run/halt/og]
# 
#  -ktcFluid -ktcSolid : thermal conductivities 
#  -ts = time-stepping-method, be=backward-Euler, fe=forward-Euler, im=implicit-multistep
#  -nc : number of correction steps for implicit predictor-corrector
# 
# Examples:  
#  ogen -noplot rpiGrid -interp=e -order=2 -factor=4 -xa=.25 -yb=3
#  $cgmp -noplot rpiCHT.cmd -g=rpiGride4.order2.hdf -nu=.05 -kappa=.01 -tf=2 -tp=.05 -solver=yale -show=rpiCHT.show -go=go
# 
# --- set default values for parameters ---
# 
$grid="innerOuter2d.hdf"; $method="ins"; $go="halt"; 
$ts="pc";  $numberOfCorrections=1; $coupled=1; $iOmega=1.; $iTol=1.e-3; $dtMax=.1;
$degreeSpace=2; $degreeTime=2;  $u0=0.; $T0=10.; $Twall=10.; 
$nu=.025; $kThermal=-1.; $ktcFluid=.05; $kappa=.04; $ktcSolid=.5; $ad2=0; 
$gravity = "0 -10. 0."; $solver="best"; $dtMax=.05; 
$rtolp=1.e-5; $atolp=1.e-7; $rtoli=1.e-7; $atoli=1.e-9; 
$fx1=1.; $fx2=2.; 
# 
$ksp="bcgs"; $pc="bjacobi"; $subksp="preonly"; $subpc="ilu"; $iluLevels=3;
#
# ----------------------------- get command line arguments ---------------------------------------
#  -- first get any commonly used options: (requires the env variable CG to be set)
$getCommonOptions = "$ENV{'CG'}/mp/cmd/getCommonOptions.h";
include $getCommonOptions
#  -- now get additional options: 
GetOptions( "degreex=i"=>\$degreex, "degreet=i"=>\$degreet,"fx1=f"=>\$fx1,"fx2=f"=>\$fx2,"T0=f"=>\$T0,"Twall=f"=>\$Twall,\
            "ad2=i"=>\$ad2 );
# -------------------------------------------------------------------------------------------------
if( $kThermal < 0 && $method eq "ins" ){ $kThermal=$nu/$prandtl; }; 
if( $kThermal < 0 && $method eq "cns" ){ $kThermal=$mu/$prandtl; }; 
# 
$grid
# 
# ------- specify fluid domain ----------
$domainName=fluidDomain; $solverName="fluid"; 
$bc = "all=noSlipWall";
$bc .= "\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$bc .= "\n bcNumber101=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber101=heatFluxInterface";
$bc .= "\n bcNumber102=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber102=heatFluxInterface";
# $bc = "all=slipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
if( $tz eq "turn off twilight zone" ){$ic = "uniform flow\n p=1., u=$u0";}
$ktc=$ktcFluid; 
$fx=$fx1; $fy=$fx1; $fz=$fx1; $ft=$fx1;
if( $method eq "ins" ){ $cmd = "include $ENV{CG}/mp/cmd/insDomain.h"; }else{ $cmd ="*"; };
$cmd
#
#  Cgcns:
$bc = "all=noSlipWall uniform(u=.0,T=$T0)\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
# gravitationally stratified : rho(y) = rho0*exp( gravity[1]/(Rg*T0) ( y - y0 ))
$rho0=1.; $y0=0.; 
if( $tz eq "turn off twilight zone" ){ $ic = "OBIC:user defined...\n gravitationally stratified\n $rho0 $y0\n r=$rho0 u=0. v=0. T=$T0\n exit";}
if( $method eq "cns" ){ $cmd = "include cnsDomain.h"; }else{ $cmd ="*"; };
$cmd
# 
# ------- specify solid domain ----------
$domainName=solidDomain; $solverName="solid"; 
# $bc = "all=dirichletBoundaryCondition, uniform(T=$Twall)\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$bc = "all=neumannBoundaryCondition";
$bc .= "\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$bc .= "\n bcNumber101=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber101=heatFluxInterface";
$bc .= "\n bcNumber102=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber102=heatFluxInterface";
if( $tz eq "turn off twilight zone" ){ $ic = "uniform flow\n" . "T=$Twall\n";}
$ktc = $ktcSolid;
$fx=$fx2; $fy=$fx2; $fz=$fx2; $ft=$fx2;
include $ENV{CG}/mp/cmd/adDomain.h
# 
continue
#
# -- set parameters for cgmp ---
# 
  final time $tFinal
  times to plot $tPlot
  cfl $cfl
  $ts
  $tz
  debug $debug
  number of PC corrections $numberOfCorrections
  OBPDE:interface tolerance $iTol
  OBPDE:interface omega $iOmega
  OBPDE:solve coupled interface equations $coupled
  OBPDE:use mixed interface conditions $mixedInterface
# 
  show file options
    compressed
      open
       $show
    frequency to flush
      2
    exit
  continue
continue
#
plot:fluid : T
echo to terminal 1
$go



 plot:solid : T


      contour
 #min max 0 $Twall
        exit
 #min max 0 $Twall
        exit


  movie mode
  finish
