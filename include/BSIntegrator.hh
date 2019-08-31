#ifndef BSIntegrator_h
#define BSIntegrator_h 1
#include"globals.hh"
#include<vector>


class BSEquation;

class BSIntegrator
{
public:
 // constructor, destructor  
  
  
  ~BSIntegrator();
   static BSIntegrator* getInstance();
   
//methods   
   inline void SetEps(double val ) {eps=val;}
   inline void SetHmin(double val ) {hmin=val;}
   inline void SetHmax(double val ) {hmax=val;}
   inline void SetDxsav(double val ) {dxsav=val;}
   inline void SetKmax(int nval ) {kmax=nval;}
   inline void SetMaxstep(int nval ) {maxstep=nval;}
   inline void SetXmax(double val ) {xmax=val;}
   inline void SetIsTheEquationOfMotionForG4Tracking(bool val )
   				{IsTheEquationOfMotionForG4Tracking=val;}
   
   
   
   inline double GetEps( ) {return eps;}
   inline double GetHmin( ) {return hmin;}
   inline int GetMaxstep() {return maxstep;}
   inline double GetHmax( ) {return hmax;}
   inline double GetXlast() {return xlast;}
   inline int GetKmax( ) {return kmax;}
   inline int GetKount( ) {return kount;}
    
   
   inline std::vector<double> GetXp_p(){return xp_p;}
   inline G4double* GetYlast(){return ylast;}
   inline std::vector< std::vector<double> > GetYp_p(){return yp_p;}
   
   inline void SetEquationTobeIntegrated (BSEquation* anEquation) 
                             {theEquation = anEquation; }
   inline BSEquation*  GetEquationTobeIntegrated() 
                             {return  theEquation;}
			     
   bool do_integration
         (G4double ystart[], const G4double x1, const G4double h1, G4int &nok,
	                                                            G4int &nbad);
   
   void bsstep(G4double y[], G4double dydx[], G4double &xx, const G4double htry,
	 const G4double yscal[], G4double &hdid, G4double &hnext);
   
   void mmid(const G4double y[], const G4double dydx[], const G4double xs,
             const G4double htot, const G4int nstep, G4double yout[]);
	
   void pzextr(const int iest, const G4double xest, const G4double yest[], 
               G4double yz[], G4double dy[]);			 
   	 
private:
   static BSIntegrator* instance;
   BSIntegrator();
   BSIntegrator(BSEquation* anEquation);
     
     
   
   
   BSEquation* theEquation;
   G4double eps,hmin,hmax,dxsav,xlast,xmax;
   G4int kmax, kount, maxstep,last_nvar;
  // G4double *xp_p;
   G4double ylast[10];
  // G4double** yp_p;
   std::vector<double> xp_p;
   std::vector< std::vector<double> > yp_p; 
   
 
 // private for extrapolation  
   /*G4double* x_p;
   G4double** d_p;*/
   std::vector< double > x_p;
   std::vector< std::vector<double> > d_p;
   
   G4int verbose;
   
 // parameter that tell if the equation of motion is BsEquationInG4
   G4bool IsTheEquationOfMotionForG4Tracking;     
   
};
#endif
