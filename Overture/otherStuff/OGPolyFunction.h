#ifndef OGPolyFunction_h
#define OGPolyFunction_h "OGPolyFunction.h"

#include "OGFunction.h"

//===========================================================================================
// Define a Function and it derivatives
//   This function can be used to define the "exact solution" for
//   an Overlapping Grid Appliciation (aka. TwilightZone Flow)
//
//      f(x,y,z,n,t)=[c(0,0,0,n)+x*(c(1,0,0,n)+x*(c(2,0,0,n)+x*(c(3,0,0,n)+c(4,0,0,n)*x)))
//                             +y*(c(0,1,0,n)+y*(c(0,2,0,n)+y*(c(0,3,0,n)+c(0,4,0,n)*y)))
//                             +z*(c(0,0,1,n)+z*(c(0,0,2,n)+z*(c(0,0,3,n)+c(0,0,4,n)*z)))]
//           *(a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n))))));
//
//    Here n is the component number and t is a parameter "time"
//    The matrices "c" and "a" defining the coefficients can be assigned
//    in the constructor or else default values can be used.
//
//===========================================================================================
class OGPolyFunction : public OGFunction
{

 private:
  enum
  {
    maximumDegreeX=10,
    maximumDegreeT=10,
    maxOptimizedDegree=6 // max degree for optimized evaluations
  };
  int numberOfComponents,degreeX,degreeY,degreeZ,degreeT,numberOfDimensions;

  RealArray cc;
  RealArray a;
  void checkArguments(const Index & N);

 public:
  

 //  
 // Create a polynomial with the given degree, number Of Space Dimensions and for
 // a maximum number of components
 // 
 OGPolyFunction(const int & degreeSpaceOfPolynomial=2, 
                const int & numberOfDimensions=3,
	        const int & numberOfComponents=10,
                const int & degreeOfTimePolynomial=1 );
  
  OGPolyFunction(const OGPolyFunction & ogp );
  OGPolyFunction & operator=(const OGPolyFunction & ogp );
  ~OGPolyFunction(){}
  
  //
  //  Supply coefficients to use
  //   note that arrays should be dimensioned
  //        c(0:d,0:d,0:d,0:nc-1) and a(0:d,0:nc-1)   where nc=numberOfComponents
  // 
  void setCoefficients( const RealArray & c, const RealArray & a );

  // get the current values of the coefficients
  void getCoefficients( RealArray & c, RealArray & a ) const;


  // ========= Here are versions with a new naming convention ===========

  // (default arguments not allowed on operators)
  virtual real operator()(const real , const real , const real);
  virtual real operator()(const real , const real , const real , const int);
  virtual real operator()(const real , const real , const real , const int, const real);

  virtual real t(const real , const real , const real , const int=0, const real=0.);
  virtual real x(const real , const real , const real , const int=0, const real=0.);
  virtual real y(const real , const real , const real , const int=0, const real=0.);
  virtual real xx(const real , const real , const real , const int=0, const real=0.);
  virtual real xy(const real , const real , const real , const int=0, const real=0.);
  virtual real yy(const real , const real , const real , const int=0, const real=0.);
  virtual real z(const real , const real , const real , const int=0, const real=0.);
  virtual real xz(const real , const real , const real , const int=0, const real=0.);
  virtual real yz(const real , const real , const real , const int=0, const real=0.);
  virtual real zz(const real , const real , const real , const int=0, const real=0.);
  virtual real xxx(const real , const real , const real , const int=0, const real=0.);
  virtual real xxxx(const real , const real , const real , const int=0, const real=0.);
  // obtain a general derivative
  virtual real gd(const int & ntd, const int & nxd, const int & nyd, const int & nzd,
                  const real , const real , const real , const int=0, const real=0. );


  virtual RealDistributedArray operator()(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N);
  virtual RealDistributedArray operator()(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
/* ---
  RealDistributedArray operator()(const MappedGrid & c, const Index & I1, const Index & I2, 
				  const Index & I3);
  RealDistributedArray operator()(const MappedGrid & c, const Index & I1, const Index & I2, 
				  const Index & I3, const int n);
  RealDistributedArray operator()(const MappedGrid & c, const Index & I1, const Index & I2, 
				  const Index & I3, const int n, const real t,
				  const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering );
----- */

  virtual RealDistributedArray t(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray x(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray y(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray z(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray xx(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray yy(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray zz(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray xy(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray xz(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray yz(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray laplacian(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray xxx(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray xxxx(const MappedGrid & c, const Index & I1, const Index & I2, 
              const Index & I3, const Index & N, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering) ;
  virtual RealDistributedArray gd(const int & ntd, const int & nxd, const int & nyd, const int & nzd,
				  const MappedGrid & c, const Index & I1, const Index & I2, 
				  const Index & I3, const Index & N, const real t=0.,
				  const GridFunctionParameters::GridFunctionType & centering
				  =GridFunctionParameters::defaultCentering);

  // obtain a general derivative using serial arrays
  virtual realSerialArray& gd( realSerialArray & result,   // put result here
			     const realSerialArray & x,  // coordinates to use if isRectangular==true
                             const int numberOfDimensions,
                             const bool isRectangular,
                             const int & ntd, const int & nxd, const int & nyd, const int & nzd,
			     const Index & I1, const Index & I2, 
			     const Index & I3, const Index & N, 
                             const real t=0., int option =0  );

  realCompositeGridFunction operator()(CompositeGrid & cg);
  realCompositeGridFunction operator()(CompositeGrid & cg, const Index & N=nullIndex);
  realCompositeGridFunction operator()(CompositeGrid & cg, const Index & N, const real t,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction t(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction x(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction y(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction z(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction xx(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction xy(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction xz(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction yy(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction yz(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction zz(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction laplacian(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction xxx(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);
  realCompositeGridFunction xxxx(CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
           const GridFunctionParameters::GridFunctionType & centering=GridFunctionParameters::defaultCentering);

  realCompositeGridFunction gd(const int & ntd, const int & nxd, const int & nyd, const int & nzd,
			       CompositeGrid & cg, const Index & N=nullIndex, const real t=0.,
			       const GridFunctionParameters::GridFunctionType & centering
                                 =GridFunctionParameters::defaultCentering);
};


#undef TIME

#endif 
