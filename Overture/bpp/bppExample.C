// This file automatically generated from bppExample.bC with bpp.
// =====================================================================
//
//  C++ example using bpp macros 
//    
// =====================================================================

#include <stdio.h>
#include <math.h>

typedef double Real;

// Inline bpp macro:

// Shows one inline macro using another:

// ---------------------------------------------------
// Macro takeStep: 
//   ORDER  : order of accuracy, 2,4
//   OPTION : STENCIL,FAME 
// ---------------------------------------------------

///Example showing how to save output to a file 
// Macro to build separate files 

// create two files 


int 
main(int argc, char *argv[])
{

    int nx=100;

    const int numGhost=2;

    int n1a=0, n1b=nx;
    int nd1a=n1a-numGhost;
    int nd1b=n1b+numGhost;
    int nd1= nd1b-nd1a+1;

    Real *up_p = new Real[nd1];  // preivous
    Real *uc_p = new Real[nd1];  // current
    Real *un_p = new Real[nd1];  // next 

    #define up(i) up_p[i-nd1a]
    #define uc(i) uc_p[i-nd1a]
    #define un(i) un_p[i-nd1a]

    for( int i=nd1a; i<=nd1b; i++ )
    {
          	up(i) = 0.;
          	uc(i) = 1.;
    }

    Real c=1.;
    Real dx = 1./nx;
    Real dt = dx/c;
    Real cdtdxSq = c*c*dt*dt/(dx*dx);

          	printf("Run order=2, option=STENCIL\n");
        for( int i=n1a; i<=n1b; i++ )// start loop i 
        {
                    		    un(i) = 2.*uc(i) - up(i) + cdtdxSq*(uc(i+1)-2.*uc(i)+uc(i-1));
        } // end loop i 
          	printf("Run order=4, option=STENCIL\n");
        for( int i=n1a; i<=n1b; i++ )// start loop i 
        {
                    		    un(i) = 2.*uc(i) - up(i) + cdtdxSq*((uc(i-1+1)-2.*uc(i-1)+uc(i-1-1))-2.*(uc(i+1)-2.*uc(i)+uc(i-1))+(uc(i+1+1)-2.*uc(i+1)+uc(i+1-1)));
        } // end loop i 
    
          	printf("Run order=2, option=FAME\n");
        for( int i=n1a; i<=n1b; i++ )// start loop i 
        {
                    		    un(i) = 2.*uc(i) - up(i) + cdtdxSq*( uc(i-1) -2.*uc(i) + uc(i+1) );
        } // end loop i 
          	printf("Run order=4, option=FAME\n");
        for( int i=n1a; i<=n1b; i++ )// start loop i 
        {
                    		    un(i) = 2.*uc(i) - up(i) + cdtdxSq*( -uc(i-2) + 4.*uc(i-1) -6.*uc(i) + 4.*uc(i+1) -uc(i+2) );
        } // end loop i 
    

    delete [] up_p;
    delete [] uc_p;
    delete [] un_p;

    printf("bppExample: done.\n");
    return 0;
}

