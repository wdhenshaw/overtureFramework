# Run ogmgt and extract timing and memory usage for ogmg/doc/performanceComparison.tex
# Examples: 
# make -s square
# make -s square2
# make -s square4

ogmgt = $(OvertureCheckout)/ogmg/ogmgt

all: square cic shapes

allSerial: square2s square4s shapes2s shapes4s	

square: square2 square4
cic: cic2 cic4
shapes: shapes2 shapes4

# order=2
square2:
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out;
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=4e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# serial: (no Hypre)
square2s:
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -cgSolver=yale -omega=1.15 -bs=rb -nsm="2 1" -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out;
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# order=4
square4:
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=6e-15 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# serial 
square4s:
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -maxIt=9 -cgSolver=yale -omega=1.15 -bs=rbj -nsm="2 1" -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square1024.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"


# order=2 finer grid:
square2b:
	$(ogmgt) tz -noplot -g=square2048.order2 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out;
	$(ogmgt) tz -noplot -g=square2048.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=4e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square2048.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square2048.order2 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=4e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# order=4 finer grid
square4b:
	$(ogmgt) tz -noplot -g=square2048.order4 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square2048.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=6e-15 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square2048.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=square2048.order4 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=1e-12 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# order=2
cic2:
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=4e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=3e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

#serial
cic2s:
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=3e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=cice32.order2.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# order=4
cic4:
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=1e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=2e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=1e-13 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# serial	
cic4s:
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=0 -debug=0 -levels=1                 -rtol=2e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=cice32.order4.ml3 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=1e-13 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"
# order=2 - note: set omega to match r2p paper
shapes2:
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=2e-15 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# serial
shapes2s:
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=1 -maxIt=9 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -g=shapesBige16.order2.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=5e-14 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"


# order=4 note: set omega to match r2p paper
shapes4: 
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=1 -maxIt=11 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=AMG   -rtol=8e-15 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=2e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=4e-13 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"

# serial
shapes4s: 
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=1 -maxIt=10 -orderCoarse=2 -debug=0 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=bcgs  -rtol=2e-13 -opave=0 -opavCoarseGrid=0 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	$(ogmgt) tz -noplot -numberOfParallelGhost=4 -g=shapesBige16.order4.ml3 -omega=1.15 -autoChoose=0 -debug=0 -levels=1 -cgSolver=gmres -rtol=4e-13 -opave=0 -opavCoarseGrid=0 -ilucg=3 -solveTwice=1 > ogmg.out; getPerfResults.p ogmg.out
	echo "\hline"
