#include "MAGCOSEquationOfMotion.hh"


////////////////////////////////////////////////////////////////////////////////
//
MAGCOSEquationOfMotion::MAGCOSEquationOfMotion( G4MagneticField* MagField )
                              : G4Mag_EqRhs( MagField ) 
{ reverse_time=false;
  direction=1.;
  Equation=&MAGCOSEquationOfMotion::LorentzMotion;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::EvaluateRhsGivenB( const G4double y[],
			                        const G4double B[3],
				                      G4double dydx[] ) const
{ (this->*Equation)( y,  B, dydx);
  return ;
}
/////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::MotionAlongBfieldLine( const G4double*,
			                            const G4double B[3],
				                          G4double dydx[] ) const
{
  G4double Bmag=direction*std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  dydx[0] = B[0]/Bmag;       
  dydx[1] = B[1]/Bmag;       
  dydx[2] = B[2]/Bmag;
   
  dydx[3]=0.;
  dydx[4]=0.;
  dydx[5]=0.; 
  return ;
}
/////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::LorentzMotion( const G4double y[],
			                    const G4double B[3],
				            G4double dydx[] ) const
{
  G4double momentum_mag_square = y[3]*y[3] +  y[4]*y[4] +  y[5]*y[5];
  G4double inv_momentum_magnitude = 
            direction * 1.0 / sqrt( momentum_mag_square );

  G4double cof =   FCof()*inv_momentum_magnitude;

  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

  dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   

  return ;
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::SetReverseTimeMode(G4bool abool)
{ reverse_time=abool;
  direction=1.;
  if (reverse_time) direction= -1.;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::SetEquationType(G4String aString)
{ if (aString=="LORENTZ_MOTION")  
          Equation=&MAGCOSEquationOfMotion::LorentzMotion;
  if (aString=="BFIELD_LINE")
           Equation=&MAGCOSEquationOfMotion::MotionAlongBfieldLine;
}




