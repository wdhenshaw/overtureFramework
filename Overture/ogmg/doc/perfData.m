cpuTotal=1; % location in array with cpu total
cpuSolve=2; % location in array with cpu total
storage =3; % location in array with storage in reals/pt
square1024 = 1; % grid name enumerator
solverName{1}="Ogmg V[2,1]";
data(cpuTotal,  2  ,1,square1024)=    0.28  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,1,square1024)=    0.27  ; 
data(storage ,  2  ,1,square1024)=     8.7 ; 
solverName{2}="AMG";
data(cpuTotal,  2  ,2,square1024)=    6.26  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,2,square1024)=    3.19  ; 
data(storage ,  2  ,2,square1024)=    99.8 ; 
solverName{3}="Bi-CG-Stab ILU(1)";
data(cpuTotal,  2  ,3,square1024)=   23.27  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,3,square1024)=   22.30  ; 
data(storage ,  2  ,3,square1024)=    64.3 ; 
solverName{4}="GMRES ILU(3)";
data(cpuTotal,  2  ,4,square1024)=   68.64  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,4,square1024)=   67.50  ; 
data(storage ,  2  ,4,square1024)=    86.4 ; 
data(cpuTotal,  4  ,1,square1024)=    0.49  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,1,square1024)=    0.48  ; 
data(storage ,  4  ,1,square1024)=     8.4 ; 
data(cpuTotal,  4  ,2,square1024)=   15.43  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,2,square1024)=   10.66  ; 
data(storage ,  4  ,2,square1024)=   178.0 ; 
data(cpuTotal,  4  ,3,square1024)=   27.42  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,3,square1024)=   25.63  ; 
data(storage ,  4  ,3,square1024)=   113.1 ; 
data(cpuTotal,  4  ,4,square1024)=   61.82  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,4,square1024)=   58.97  ; 
data(storage ,  4  ,4,square1024)=   155.8 ; 
cice32 = 2; % grid name enumerator
data(cpuTotal,  2  ,1,cice32)=    0.71  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,1,cice32)=    0.67  ; 
data(storage ,  2  ,1,cice32)=     9.7 ; 
data(cpuTotal,  2  ,2,cice32)=   20.53  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,2,cice32)=   15.58  ; 
data(storage ,  2  ,2,cice32)=   115.2 ; 
data(cpuTotal,  2  ,3,cice32)=   60.80  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,3,cice32)=   59.27  ; 
data(storage ,  2  ,3,cice32)=    62.0 ; 
data(cpuTotal,  2  ,4,cice32)=  128.30  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,4,cice32)=  126.50  ; 
data(storage ,  2  ,4,cice32)=    83.7 ; 
data(cpuTotal,  4  ,1,cice32)=    1.37  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,1,cice32)=    1.31  ; 
data(storage ,  4  ,1,cice32)=    10.8 ; 
data(cpuTotal,  4  ,2,cice32)=  100.31  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,2,cice32)=   92.72  ; 
data(storage ,  4  ,2,cice32)=   172.9 ; 
data(cpuTotal,  4  ,3,cice32)=   59.14  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,3,cice32)=   56.36  ; 
data(storage ,  4  ,3,cice32)=   109.4 ; 
data(cpuTotal,  4  ,4,cice32)=   90.67  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,4,cice32)=   86.09  ; 
data(storage ,  4  ,4,cice32)=   151.4  ; 
shapesBige16 = 3; % grid name enumerator
data(cpuTotal,  2  ,1,shapesBige16)=    0.76  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,1,shapesBige16)=    0.62  ; 
data(storage ,  2  ,1,shapesBige16)=    14.9 ; 
data(cpuTotal,  2  ,2,shapesBige16)=    8.22  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,2,shapesBige16)=    6.27  ; 
data(storage ,  2  ,2,shapesBige16)=   112.6 ; 
data(cpuTotal,  2  ,3,shapesBige16)=   11.67  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,3,shapesBige16)=   11.10  ; 
data(storage ,  2  ,3,shapesBige16)=    61.9 ; 
data(cpuTotal,  2  ,4,shapesBige16)=   22.20  ; % data(:,order,solver,grid)= value
data(cpuSolve,  2  ,4,shapesBige16)=   21.56  ; 
data(storage ,  2  ,4,shapesBige16)=    83.6 ; 
data(cpuTotal,  4  ,1,shapesBige16)=    1.51  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,1,shapesBige16)=    1.27  ; 
data(storage ,  4  ,1,shapesBige16)=    20.3 ; 
data(cpuTotal,  4  ,2,shapesBige16)=   69.68  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,2,shapesBige16)=   66.39  ; 
data(storage ,  4  ,2,shapesBige16)=   168.1  ; 
data(cpuTotal,  4  ,3,shapesBige16)=   11.74  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,3,shapesBige16)=   10.69  ; 
data(storage ,  4  ,3,shapesBige16)=   107.9 ; 
data(cpuTotal,  4  ,4,shapesBige16)=   16.17  ; % data(:,order,solver,grid)= value
data(cpuSolve,  4  ,4,shapesBige16)=   14.28  ; 
data(storage ,  4  ,4,shapesBige16)=   149.0 ; 
