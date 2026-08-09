#
# cgmp: solve for INS/CNS and AD in two or more domains
#     : flow around a heated cylinder (2d or 3d)
#
echo to terminal 0
# 
# Usage:
#    cgmp [-noplot] io -g=<name> -method=[ins|cns] -nu=<> -mu=<> -kappa=<num> -tf=<tFinal> -tp=<tPlot> ...
#           -solver=<yale/best> -ktcFluid=<> -ktcFluid=<> -tz=[poly/trig/none] -bg=<backGroundGrid> ...
#           -nc=<num> degreex=<num> -degreet=<num> -ts=[fe|be|im|pc] -coupled=[0|1]  -mixedInterface=[0|1] ...
#           -go=[run/halt/og]
# 
#  -ktcFluid -ktcSolid : thermal conductivities 
#  -ts = time-stepping-method, be=backward-Euler, fe=forward-Euler, im=implicit-multistep
#  -nc : number of correction steps for implicit predictor-corrector
# 
# Examples:  
# -- incompressible NS examaples
#  $cgmp movingHeatedBody.cmd -g=solidAnnulusInASquareGride2.order2.hdf -method=ins -nu=.1 -kappa=.05 -tf=1.0 -tp=0.01 -solver=yale -show=movingHeatedAnnulus.show
#
# -- compressible NS examples
#  $cgmp movingHeatedBody.cmd -g=solidAnnulusInASquareGride2.order2.hdf -method=cns -mu=.05 -T0=300 -Twall=400 -tp=.01
#
# --- set default values for parameters ---
# 
$grid="innerOuter2d.hdf"; $method="ins"; $go="halt"; 
$ts="pc";  $numberOfCorrections=1; $coupled=1; $iOmega=1.; $iTol=1.e-3; $dtMax=.1;
$degreeSpace=2; $degreeTime=2; $T0=10.; $Twall=10.; 
$u0=0.;  $Tf0=0; # initial velocity and temp in the fluid
$nu=.025; $kThermal=-1.; $ktcFluid=.05; $kappa=.04; $ktcSolid=.5; $ad2=0; 
$bcTypeFluid="wall"; # [wall|inflow]
$uInflow=0.;             # inflow velocity
$gravity = "0 -10. 0."; $solver="best"; $dtMax=.05; $psolver="best"; $pc="ilu"; $ogmgCoarseGridSolver="yale"; $ogmgIlucgFill=20; $ogmgAutoChoose=0;
$rtolp=1.e-5; $atolp=1.e-7; $rtoli=1.e-7; $atoli=1.e-9; 
$fx1=1.; $fx2=2.; 
$flushFrequency=20;
#
# Translate along a line
$tx=0.707107; $ty=0.707107;
# Sinsusoid motion:
$amps=0.5; 
# 
$ksp="bcgs"; $pc="bjacobi"; $subksp="preonly"; $subpc="ilu"; $iluLevels=3;
#
# ----------------------------- get command line arguments ---------------------------------------
#  -- first get any commonly used options: (requires the env variable CG to be set)
$getCommonOptions = "$ENV{'CG'}/mp/cmd/getCommonOptions.h";
include $getCommonOptions
#  -- now get additional options: 
GetOptions( "degreex=i"=>\$degreex, "degreet=i"=>\$degreet,"fx1=f"=>\$fx1,"fx2=f"=>\$fx2,"T0=f"=>\$T0,"Twall=f"=>\$Twall,\
            "ad2=i"=>\$ad2, "method=s"=>\$method , "flushFrequency=s"=>\$flushFrequency, "solver=s"=>\$solver, "psolver=s"=>\$psolver,\
            "ogmgAutoChoose=i"=>\$ogmgAutoChoose,  "bcTypeFluid=s"=>\$bcTypeFluid,"uInflow=f"=>\$uInflow,\
            "tx=f"=>\$tx,"ty=f"=>\$ty ,"amps=f"=>\$amps   );
# -------------------------------------------------------------------------------------------------
if( $kThermal < 0 && $method eq "ins" ){ $kThermal=$nu/$prandtl; }; 
if( $kThermal < 0 && $method eq "cns" ){ $kThermal=$mu/$prandtl; }; 
# 
$grid
#
# NOTE: these MOVE COMMANDS ARE USED IN adDomain.h and insDomain.h below
# Sinusoid time function:  f(t)=b0*sin(2.*Pi*f0*(t-t0));
$moveCmds = \
  "turn on moving grids\n" . \
  "specify grids to move\n" . \
  "    matrix motion\n" . \
  "      translate along a line\n" . \
  "        point on line: 0 0 0\n" . \
  "        tangent to line: $tx,$ty,0  \n" . \
  "        edit time function\n" . \
  "         sinusoidal function\n" . \
  "         sinusoid parameters: $amps,1,0 (b0,f0,t0)\n" . \
  "        done\n" . \
  "      done\n" . \
  "      choose grids by share flag\n" . \
  "         100 \n" . \
  "   done\n" . \
  "done";
# 
# ------- specify fluid domain ----------
$domainName=outerDomain; $solverName="fluid"; 
$bc = "all=noSlipWall\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
# try inflow outflow
$halfH=0.2; $cpn=1; 
$bcInflow = "all=noSlipWall\n  bcNumber1=inflowWithVelocityGiven, parabolic(d=$halfH, p=1.,u=$uInflow,T=$Tf0)\n bcNumber2=outflow, pressure(10.*p+$cpn*p.n=0.)\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
#
if( $bcTypeFluid eq "inflow" ){ $bc=$bcInflow; }
#
# $bc = "all=slipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
if( $tz eq "turn off twilight zone" ){$ic = "uniform flow\n p=1., u=$uInflow, T=$Tf0";}
$ktc=$ktcFluid; 
$fx=$fx1; $fy=$fx1; $fz=$fx1; $ft=$fx1;
if( $method eq "ins" ){ $cmd = "include $ENV{CG}/mp/cmd/insDomain.h"; }else{ $cmd ="*"; };
$cmd
#
#
#  Cgcns:
$bc = "all=noSlipWall uniform(u=.0,T=$T0)\n bcNumber100=noSlipWall, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
# gravitationally stratified : rho(y) = rho0*exp( gravity[1]/(Rg*T0) ( y - y0 ))
$rho0=1.; $y0=0.; 
if( $tz eq "turn off twilight zone" ){ $ic = "OBIC:user defined...\n gravitationally stratified\n $rho0 $y0\n r=$rho0 u=$uInflow v=0. T=$T0\n exit";}
if( $method eq "cns" ){ $cmd = "include $ENV{CG}/mp/cmd/cnsDomain.h"; }else{ $cmd ="*"; };
$cmd
#
# ------- specify solid domain ----------
$domainName=innerDomain; $solverName="solid"; 
# $bc = "all=dirichletBoundaryCondition, uniform(T=$Twall)\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$bc = "all=neumannBoundaryCondition\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
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
      $flushFrequency
    exit
  continue
continue
#
plot:fluid : T
contour
  plot contour lines (toggle)
  exit
  plot contour lines (toggle)
  exit
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
