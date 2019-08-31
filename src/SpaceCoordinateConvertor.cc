#include "SpaceCoordinateConvertor.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "magneto_fsubroutine_def.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "MAGCOSUnits.hh"
SpaceCoordinateConvertor* SpaceCoordinateConvertor::instance = 0;

SpaceCoordinateConvertor::SpaceCoordinateConvertor()
{  // Set the start date
  ReferenceDate = DateAndTime(2001,1,1,12,0,0);
  
  //pointers
  p_hh_nm=0;
  p_gg_nm=0;
  

  verbose =0;
  Initialise();
  
}

////////////////////////////////////////////////////////////////////////

SpaceCoordinateConvertor::~SpaceCoordinateConvertor()
{
   ;
}

////////////////////////////////////////////////////////////////////////////////
//
SpaceCoordinateConvertor* SpaceCoordinateConvertor::getInstance()
{if (instance == 0) instance = new SpaceCoordinateConvertor;
 return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::SetSystemInAndOut
                                            (G4String sys_in, G4String sys_out)
{ if (sys_in != "GSM" &&  sys_in != "GEO" && sys_in != "SM" &&  sys_in != "MAG"
     && sys_in != "GSE" && sys_in != "GEI")
                      G4cout<<sys_in<<" is not a good system of coordinate"<<G4endl;  
 
  else if (sys_out != "GSM" &&  sys_out != "GEO" && sys_out != "SM" &&  sys_out != "MAG"
          && sys_in != "GSE" && sys_in != "GEI")
                       G4cout<<sys_out<<" is not a good system of coordinate"<<G4endl;  
    
  else{ 
  	system_in=sys_in;
  	system_out=sys_out;    
 
  	ComputeSelectedMatrix();
  }  
}

////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinateConvertor::Transform(G4ThreeVector vec_in,
	                           G4String sys_in,
		                   G4String sys_out)
{

 if (sys_in != "GSM" &&  sys_in != "GEO" && sys_in != "SM" &&
     sys_in != "MAG" &&  sys_in != "GEI" && sys_in != "GSE"){
   		G4cout<<sys_in<<" is not a good system of coordinate"<<G4endl;  
    		return G4ThreeVector(vec_in);
 }
 
 if (sys_out != "GSM" &&  sys_out != "GEO" && sys_out != "SM" &&
     sys_out != "MAG" &&  sys_out != "GEI" && sys_out != "GSE"){
 	G4cout<<sys_out<<" is not a good system of coordinate"<<G4endl;  
    	return G4ThreeVector(vec_in);
 } 
 
 system_in=sys_in;
 system_out=sys_out;    
 
 ComputeSelectedMatrix();
 
 return Transform(vec_in);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinateConvertor::Transform(G4ThreeVector vec_in){
	if (system_in == system_out) return G4ThreeVector(vec_in);
 	else return SelectedMatrix*vec_in;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinateConvertor::TransformGEOinGSM(G4ThreeVector vec_in){
	return GEOtoGSM.inverse()*vec_in;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinateConvertor::TransformGSMinGEO(G4ThreeVector vec_in){
	return GEOtoGSM*vec_in;
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGEOPositionFromGEOID
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4ThreeVector& position)
{ G4double delta=90.*degree-latitude;
  G4double sind=sin(delta);
  G4double cosd=cos(delta);
  G4double p=longitude;
  G4double eta=sqrt(GeoidLargeSemiaxis2*sind*sind
                                  +GeoidSmallSemiaxis2*cosd*cosd);
  G4double xi=(GeoidLargeSemiaxis2-GeoidSmallSemiaxis2)*sind*cosd/eta;
  G4double r=sqrt(xi*xi+(eta+altitude)*(eta+altitude));
  G4double sbeta=xi/r;
  G4double cbeta=(eta+altitude)/r;
  G4double sth=cbeta*sind+sbeta*cosd;
  G4double cth=-sbeta*sind+cbeta*cosd;
  G4double sp=sin(p);
  G4double cp=cos(p);
  position=G4ThreeVector(sth*cp,sth*sp,cth);
  //G4cout<<position<<std::endl;
  position.setMag(r);
  //G4cout<<position/Re<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGEODirectionAndPositionFromGEOID
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4double zenith, G4double azimuth,
		     G4ThreeVector& position, G4ThreeVector& direction)
{ G4double delta=90.*degree-latitude;
  G4double sind=sin(delta);
  G4double cosd=cos(delta);
  G4double p=longitude;
  G4double eta=sqrt(GeoidLargeSemiaxis2*sind*sind
                                  +GeoidSmallSemiaxis2*cosd*cosd);
  G4double xi=(GeoidLargeSemiaxis2-GeoidSmallSemiaxis2)*sind*cosd/eta;
  G4double r=sqrt(xi*xi+(eta+altitude)*(eta+altitude));
  G4double sbeta=xi/r;
  G4double cbeta=(eta+altitude)/r;
  G4double sth=cbeta*sind+sbeta*cosd;
  G4double cth=-sbeta*sind+cbeta*cosd;
  G4double sp=sin(p);
  G4double cp=cos(p);
  position=G4ThreeVector(sth*cp,sth*sp,cth);
  position.setMag(r);
  G4double t=position.theta();
  G4ThreeVector VGeoid=G4ThreeVector(0.,0.,1.);
  VGeoid.setTheta(zenith);
  VGeoid.setPhi(180.*degree-azimuth);
  G4double VgeoR,VgeoT,VgeoP;
  VgeoR=cbeta*VGeoid.z()+sbeta*VGeoid.x();
  VgeoT=-sbeta*VGeoid.z()+cbeta*VGeoid.x();
  VgeoP=VGeoid.y();
  G4ThreeVector VGeoSphere=G4ThreeVector(VgeoT,VgeoP,VgeoR);
  direction=-VGeoSphere.rotateY(t).rotateZ(p);
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGEOIDCoordinatesFromGEOPosition
	             (const G4ThreeVector GEOposition, 
		      G4double& altitude, G4double& longitude, G4double& latitude) 	 	
	
{ G4double rho=GEOposition.rho();
  G4double rho2=rho*rho;
  G4double z=std::abs(GEOposition.z());
  G4double z2=z*z;
  G4double r1=std::sqrt(GeoidLargeSemiaxis2);
  G4double r2=std::sqrt(GeoidSmallSemiaxis2);
 
  G4double precision=.0000001;
  G4double cosa=rho/GEOposition.r();
  G4double sina=z/GEOposition.r();
  G4double alpha=std::acos(cosa); 
 
 
 //compute first b and t 
 
  G4double a=cosa*cosa/GeoidLargeSemiaxis2;
  a+= sina*sina/GeoidSmallSemiaxis2;
  G4double b = cosa*rho /  GeoidLargeSemiaxis2;
  b+= sina*z/GeoidSmallSemiaxis2;
  b*=2.;
 
  G4double c = rho2 /  GeoidLargeSemiaxis2;
  c+= z2/GeoidSmallSemiaxis2;
  c-=1.;
 
  G4double sqrt_rhot=std::sqrt( b*b  - (4.*a*c)) ;
  G4double t= ( - b + sqrt_rhot)/(2.*a);   
  G4double cosae= (rho + t*cosa) /r1;
  G4double sinae= (z + t*sina)/r2;
  G4double   cosb=r2*cosae;
  G4double   sinb=r1*sinae;
  G4double rb=std::sqrt(cosb *cosb + sinb *sinb);
  cosb/=rb;
  sinb/=rb;
  G4double  beta=(std::acos(cosb));
 
 
 
  if (t >= 0)   { // point is below the earth surface
   	cosa=cosb;
   	sina=sinb;
   	while( std::abs(alpha-beta)>precision){ 
   		alpha=std::acos(cosa); 
    		a=cosa*cosa/GeoidLargeSemiaxis2;
    		a+= sina*sina/GeoidSmallSemiaxis2;
    		b = cosa*rho /  GeoidLargeSemiaxis2;
    		b+= sina*z/GeoidSmallSemiaxis2;
    		b*=2.;
    		sqrt_rhot=std::sqrt( b*b  - (4.*a*c)) ;
    		t= ( - b + sqrt_rhot)/(2.*a);   
    		cosae= (rho + t*cosa) /r1;
    		sinae= (z + t*sina)/r2;
    		cosb=r2*cosae;
    		sinb=r1*sinae;
    		rb=std::sqrt(cosb *cosb + sinb *sinb);
    		cosb/=rb;
    		sinb/=rb;
    		beta=(std::acos(cosb));
    		cosa=cosb;
    		sina=sinb;
   	}
 }
 else{//point is above the EarthSurface 
 	G4ThreeVector direction1 =
        G4ThreeVector(rho,0.,z)-G4ThreeVector(r1,0.,0.);
 	G4ThreeVector direction2 =
        G4ThreeVector(rho,0.,z)-G4ThreeVector(0.,0.,r2);
  
  	G4double max_alpha=std::max(90.*degree-direction1.theta(),
                              90.*degree-direction2.theta());
  
  	G4double alpha1,alpha2;
  	alpha1=alpha;
    
    
    	if (beta > max_alpha){
      		cosa=cos(max_alpha);
       		sina=sin(max_alpha);
       		alpha2=max_alpha;
      	}
    	else {
     		cosa=cosb;
      		sina=sinb;
      		alpha2=beta;
     	}  
   	while( std::abs(alpha-beta)>precision){ 
   		alpha=alpha2; 
    		a=cosa*cosa/GeoidLargeSemiaxis2;
    		a+= sina*sina/GeoidSmallSemiaxis2;
    		b = cosa*rho /  GeoidLargeSemiaxis2;
    		b+= sina*z/GeoidSmallSemiaxis2;
    		b*=2.;
    		sqrt_rhot=std::sqrt( b*b  - (4.*a*c)) ;
    		t= ( - b + sqrt_rhot)/(2.*a);   
    		cosae= (rho + t*cosa) /r1;
    		sinae= (z + t*sina)/r2;
    		cosb=r2*cosae;
    		sinb=r1*sinae;
    		rb=std::sqrt(cosb *cosb + sinb *sinb);
    		cosb/=rb;
    		sinb/=rb;
    		beta=(std::acos(cosb));
    		if (alpha <beta){
     			if (beta > max_alpha){
      				cosa=cos(max_alpha);
       				sina=sin(max_alpha);
       				alpha1=alpha;
       				alpha2=alpha+max_alpha;
       			}
      			else {
      				cosa=cosb;
       				sina=sinb;
       				alpha1=alpha;
       				alpha2=beta;
       			}
     		}
    		else{
     			alpha2=(alpha1+alpha2)/2.;
       			cosa=cos(alpha2);
       			sina=sin(alpha2);
    
     		}
    	}	   
  }
  
  altitude = -t;
  latitude = beta;
  if (GEOposition.z() < 0.)  latitude = -latitude;
  longitude=GEOposition.phi();
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGEOIDCoordinatesFromGEOPositionAndDirection
	             (const G4ThreeVector GEOposition, const G4ThreeVector GEOdirection,
		      G4double& altitude, G4double& longitude, G4double& latitude,
		      G4double& zenith,G4double& azimuth) 	 	
	
{ ComputeGEOIDCoordinatesFromGEOPosition
	             (GEOposition, 
		      altitude, 
		      longitude, 
		      latitude);
  
  
  G4ThreeVector vertical(std::cos(latitude),0.,std::sin(latitude));
  vertical.setPhi(longitude);
  
  G4ThreeVector direction = GEOdirection;
  direction.rotateZ(-vertical.phi());
  direction.rotateY(-vertical.theta());
  direction=-direction;
  //zenith and azimuth give the view angle of the incomimg direction
  zenith=direction.theta();
  azimuth=direction.phi();
  
  
  		       	 	
	
}
/////////////////////////////////////////////////////////////////////////////
//
double SpaceCoordinateConvertor::
	ComputeGEODipoleInvariantLatitude(G4String sys, G4ThreeVector position)
{G4ThreeVector GEOPosition = Transform(position,sys,G4String("GEO"));
 GEOPosition = GEOPosition - geodipole_schift_in_GEO;
 G4ThreeVector MAGPosition = Transform(GEOPosition,
 				       G4String("GEO"),
				       G4String("MAG"));
				       
 G4double r = MAGPosition.mag();
 G4double sin_lat =  MAGPosition.z()/r;
 G4double cos_lat2 = 1. - sin_lat* sin_lat;
 G4double cos_Ilat2 = cos_lat2*Re/r; 
 if (cos_Ilat2 > 1 ) return -999999.*degree;  
 else return std::acos(std::sqrt(cos_Ilat2)); 				       
  
}	
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeSunPosition()
{ 
  G4int day_of_year = ReferenceDate.DayOfYear() + 1;
  //G4cout<<day_of_year<<std::endl;
  G4int year= ReferenceDate.year;
  G4int hour= ReferenceDate.hour;
  G4int minute= ReferenceDate.min;
  G4int second= ReferenceDate.sec;
  if (verbose >0) {
  	G4cout<<"Date "<<ReferenceDate.year<<'\t';
  	G4cout<<ReferenceDate.month<<'\t';
  	G4cout<<ReferenceDate.day<<'\t';
  	G4cout<<ReferenceDate.hour<<'\t';
  	G4cout<<minute<<'\t';
  	G4cout<<second<<std::endl;
  }
  if ( year < 1901 || year > 2099){
  	G4cout<<"The year is not valid "<<G4endl;
      	G4cout<<"It should be after 1900 and before 2100"<<G4endl;
      	return;
  }
  
  //New formula for gst based on Astronomical almanach
  DateAndTime dt1 = DateAndTime(2000,1,1,12,0,0);
  DateAndTime dt2=  DateAndTime(ReferenceDate.year,
  				ReferenceDate.month,
				ReferenceDate.day,0,0,0);
  				
  G4double T0 = dt2.DifferenceInDays(dt1)/36525.;
  //G4cout<<"T0 "<<T0<<std::endl;
  //G4double T=ReferenceDate.DifferenceInDays(dt1)/36525.;
  G4double fday = ReferenceDate.DifferenceInDays(dt2);
  //G4cout<<"fday "<<fday<<std::endl;
  G4double UT = fday*24.;
  //G4cout<<"UT "<<UT<<std::endl;
  
  gst = 6.697374558 + 1.0027379093* UT +
  	((8640184.812866+(0.093104-6.2e-6*T0)*T0)*T0)/3600.;
  gst *= 15.*degree;
  gst = std::fmod (gst,360.*degree);
  if (gst <0)	gst+=360.*degree;
  //G4cout<<"Almanach gst "<<gst/degree<<std::endl;	
  
  //Old formula for gst from Russel
    
  G4double fraction_of_day= (hour*3600. + minute *60. + second)/86400.; 
  G4double dj = 365. * (year -1900) + (year -1901)/4
                                    +  day_of_year - 0.5 + fraction_of_day ;                 
  G4double t = dj / 36525.;
  G4double vl = fmod(279.696678+0.9856473354*dj , 360. );
  
  // gst = fmod(279.690983+.9856473354*dj+360.*fraction_of_day+180., 360.) * degree;
 // G4cout<<"Russel gst "<<fmod(279.690983+.9856473354*dj+360.*fraction_of_day+180.,360.)<<std::endl;
  G4double g =
          fmod(358.475845+0.985600267*dj, 360.) * degree;
  
  sun_long = (vl + (1.91946-0.004789*t)*sin(g)
                          +  0.020094*sin(2.*g)) * degree;
  if (sun_long > 6.2831853) sun_long-=6.2831853;
  if (sun_long < 0) sun_long+=6.2831853;
  obliq = (23.45229-0.0130125*t) * degree;
  G4double sob = sin(obliq);
  G4double slp = sun_long - 9.924E-5;
  // the last constant is a correction for the angular aberration  due to
  //   the orbital motion of the earth  

  G4double sind = sob*sin(slp);
  G4double cosd = sqrt(1.-sind * sind);
  G4double sc = sind / cosd;
  sun_dec = atan(sc);
  sun_rasn = 3.141592654-atan2(cos(obliq)/sob*sc,-cos(slp)/cosd);
  if (verbose >0 ) {
  	G4cout.precision(10);
  	G4cout<<"gst "<< gst/degree << G4endl;
  	G4cout<<"sun_rasn "<< sun_rasn/degree << G4endl;
  	G4cout<<"sun_long "<< sun_long/degree << G4endl;
  	G4cout<<"sun_dec "<< sun_dec/degree << G4endl;
  }
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::SetReferenceDate(G4int year,G4int month ,G4int day,
	              G4int hour,G4int minute,G4int second)
{ ReferenceDate =DateAndTime(year,month,day,hour,minute,second);
  Initialise();
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::SetReferenceDate(DateAndTime ref_date)
{ ReferenceDate = ref_date;
  Initialise();
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGeodipoleParameters()
{
 G4double f1=0,f2=0,g10=0,h11=0,g11=0,g20=0,g21=0,g22=0,h21=0,h22=0,dt=0; //1st and 2nd order IGRF parameters 
 G4int year = ReferenceDate.year;
  G4int day_of_year = ReferenceDate.DayOfYear();
 if (year < 1965 || year >2005)
     G4cout<<"The year should be in the period 1965-2005"<<G4endl;    
  
 else if (year < 1970)                      // 1965-1975
     {f2=(double(year)+double(day_of_year)/365.-1965.)/5.;
      f1=1.-f2;
      g10=-30334.*f1-30220.*f2; 
      g11=-2119.*f1-2068.*f2;
      h11=5776.*f1+5737.*f2;
      g20=-1662.*f1-1781.*f2;
      g21= 2997.*f1+3000.0*f2;
      h21=-2016.*f1-2047.*f2;
      g22=1594.*f1+1611.*f2;
      h22=114.*f1+25.*f2;}
 
 else if (year < 1975)                       //1970-1975
     {f2=(double(year)+double (day_of_year)/365.-1970.)/5.;
      f1=1.-f2;
      g10=-30220.*f1- 30100.*f2; 
      g11=-2068.*f1 - 2013.*f2;
      h11=5737.*f1+5675.*f2;
      g20=-1781.*f1-1902.*f2;
      g21= 3000.0*f1+3010.*f2;
      h21=-2047.*f1-2067.*f2;
      g22=1611.*f1 + 1632.*f2;
      h22=25.*f1 - 68.*f2;}
      
 else if (year < 1980)                     //1975-1980
     {f2=(double(year)+double (day_of_year)/365.-1975.)/5.;
      f1=1.-f2;
      g10=-30100.*f1-29992.*f2; 
      g11=-2013.*f1-1956.*f2;
      h11=5675.*f1+5604.*f2;
      g20=-1902.*f1 -1997. *f2;
      g21= 3010.*f1 + 3027. *f2 ;
      h21=-2067.*f1 - 2129. *f2 ;
      g22=1632.*f1 + 1663. * f2;
      h22=- 68.*f1 - 200. * f2 ;}
       
 else if (year < 1985)                       //1980-1985
     {f2=(double(year)+double (day_of_year)/365.-1980.)/5.;
      f1=1.-f2;
      g10=-29992.*f1-29873.*f2; 
      g11=-1956.*f1-1905.*f2;
      h11=5604.*f1+5500.*f2;
      g20=-1997. *f1 -2072. * f2;
      g21= 3027. *f1 +3044. * f2;
      h21=-2129. *f1 -2197. * f2;
      g22=1663. * f1 +1687. * f2;
      h22=-200. * f1 -306 * f2;}
	 
 else if (year < 1990)                      //1985-1990
     {f2=(double(year)+float(day_of_year)/365.-1985.)/5.;
      f1=1.-f2;
      g10=-29873.*f1-29775.*f2;
      g11=-1905.*f1-1848.*f2;
      h11=5500.*f1+5406.*f2;
      g20=-2072. * f1 -2131.*f2 ;
      g21= 3044. * f1 + 3059 *f2;
      h21=-2197. * f1 -2279. * f2;
      g22=1687. * f1 + 1686. * f2;
      h22=-306 * f1 -373. * f2;}

 else if (year < 1995)                     //1990-1995
     {f2=(double(year)+double(day_of_year)/365.-1990.)/5.;
      f1=1.-f2;
      g10=-29775.*f1-29682.*f2;
      g11=-1848.*f1-1789.*f2;
      h11=5406.*f1+5318.*f2;
      g20=-2197.*f2-2131.*f1;
      g21=3074.*f2+3059.*f1;
      h21=-2356.*f2-2279.*f1;
      g22=1685.*f2+1686.*f1; 
      h22=-425.*f2-373.*f1;}
 
 else if (year < 2000)                      //1995-2000
     {f2=(double(year)+double(day_of_year)/365.-1995.)/5.;
      f1=1.-f2;
      g10=-29682.*f1-29615.*f2;
      g11=-1789.*f1-1728.*f2;
      h11=5318.*f1+5186.*f2;
      g20=-2197.*f1-2267.*f2;
      g21=3074.*f1+3072.*f2;
      h21=-2356.*f1-2478.*f2;
      g22=1685.*f1+1672.*f2; 
      h22=-425.*f1-458.*f2;}
      

     
                                    
//
//  linear extrapolation beyond 2000 
//
 else {dt=double(year)+double(day_of_year)/366.-2000.;
       g10=-29615.+14.6*dt;
       g11=-1728.+10.7*dt;
       h11=5186.-22.5*dt;
       g20=-2267.-12.4*dt;
       g21=3072.+1.1*dt;
       h21=-2478.-20.6*dt;
       g22=1672.-1.1*dt;
       h22=-458.-9.6*dt;
       }
       
       
//  geomagnetic dipole and magnitude
     
     geodipole_axis_in_GEO=-G4ThreeVector(g11,h11,g10)*tesla*1.e-9;
     geodipole_B0=geodipole_axis_in_GEO.mag();
     geodipole_axis_in_GEO/=geodipole_B0;
     G4double rd2=geodipole_B0*geodipole_B0/(nT*nT);
     
     
    /* G4cout<<"DIP B0 "<<geodipole_B0<<G4endl;
     G4cout<<"DIP THETA "<<geodipole_axis_in_GEO.theta()/degree<<G4endl;
     G4cout<<"DIP PHI "<<geodipole_axis_in_GEO.phi()/degree<<G4endl;*/
//computation of DipoleSchift in GEO according to 
//Langel, Geomagnetism 1, p 381-390 with a correction (-sign) in equation 228b p386
//Fraser-Smith, Rev. Geophys., 25, 1,1-16,1987 there the formula seems correct
//Spenvis bacgroubd url: 
//   http://www.spenvis.oma.be/spenvis/help/background/magfield/cd.html    
     G4double sqr3=std::sqrt(3.); 
     G4double L0,L1,L2,T;
     L0=2.*g10*g20+sqr3*(g11*g21+h11*h21);
     L1=-g11*g20+sqr3*(g10*g21+g11*g22+h11*h22);
     L2=-h11*g20+sqr3*(g10*h21+g11*h22-h11*g22);
     T=(L0*g10+L1*g11+L2*h11)/(4.*rd2);
     geodipole_schift_in_GEO=G4ThreeVector(re*(L1-g11*T)/(3.*rd2),
                            re*(L2-h11*T)/(3.*rd2),
			    re*(L0-g10*T)/(3.*rd2));
        
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeGeodipoleParameters1()
{
 G4double g10,h11,g11,g20,g21,g22,h21,h22; //1st and 2nd order IGRF parameters 
 g10=-29615.;
 h11=5186.;
 g11=-1728.;
 g20=-2267;
 g21=3072;
 g22=1672;
 h21=-2478;
 h22=-458;



 if (p_gg_nm )
      {g10=(*p_gg_nm)[0];
       g11=(*p_gg_nm)[1];
       h11=(*p_hh_nm)[1];
       g20=(*p_gg_nm)[2]/1.5;
       g21=(*p_gg_nm)[3]/(std::sqrt(4./3.) * 1.5);
       h21=(*p_hh_nm)[3]/(std::sqrt(4./3.) * 1.5);
       g22=(*p_gg_nm)[4]/(std::sqrt(1./3.) * 1.5);
       h22=(*p_hh_nm)[4]/(std::sqrt(1./3.) * 1.5);
      }

//  geomagnetic dipole and magnitude
     
 geodipole_axis_in_GEO=-G4ThreeVector(g11,h11,g10)*tesla*1.e-9;
 geodipole_B0=geodipole_axis_in_GEO.mag();
 geodipole_axis_in_GEO/=geodipole_B0;
 G4double rd2=geodipole_B0*geodipole_B0/(nT*nT);
     
     
 
//computation of DipoleSchift in GEO according to 
//Langel, Geomagnetism 1, p 381-390 with a correction (-sign) in equation 228b p386
//Fraser-Smith, Rev. Geophys., 25, 1,1-16,1987 there the formula seems correct
//See Spenvis background url: 
//   http://www.spenvis.oma.be/spenvis/help/background/magfield/cd.html   
   
 G4double sqr3=std::sqrt(3.); 
 G4double L0,L1,L2,T;
 L0=2.*g10*g20+sqr3*(g11*g21+h11*h21);
 L1=-g11*g20+sqr3*(g10*g21+g11*g22+h11*h22);
 L2=-h11*g20+sqr3*(g10*h21+g11*h22-h11*g22);
 T=(L0*g10+L1*g11+L2*h11)/(4.*rd2);
 geodipole_schift_in_GEO=G4ThreeVector(re*(L1-g11*T)/(3.*rd2),
                                       re*(L2-h11*T)/(3.*rd2),
			               re*(L0-g10*T)/(3.*rd2));
        
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeTranformationMatrices()
{
  // Caution!!
  //G4RotationMatrix defined by euler angles phi,theta,psi represents:
  //1st a  rotation of angle -phi around Z
  //2nd a rotation of angle -theta around X
  //3rd a rotation of angle -psi  around Z
 
 
 
  //Caution!!
  // If a rotation matrix M12 rotates  the cartesian axis of system 1 such that 
  // that after the rotation these axes coincide with the axis of system 2 then
  // for any vector V we have  
  // V1 =M12  * V2 and  V2= M12.inverse() * V1 
  //  where V1 and V2 are 
  // G4ThreeVector objects  defining the components of V in the system 1 and 2
  //,respectively.
  // this should be changed
  
 
 
  // GEI and GEO system have the same Z axis. 
  // The Greenwich meridian is in the the XZ plane of the GEO coordinate system.
  // GEI X axis is defined by the first point of Aries.  
  // gst represents the angle between Greenwich meridian and the 1st point of 
  // Aries measured eastward from the 1st point of Aries
  // Therefore the Matrix which transform   the GEI axes into  the GEO 
  // axes is a rotation  of gst  around Z 
    
    
  GEItoGEO = G4RotationMatrix(G4ThreeVector(0.,0.,1.),gst);



  // GEO and MAG
  // MAG Z axis is the gomagnetic dipole axis
  // MAG Y axis is perpendicular to the Z GEO axis
   
  G4ThreeVector YMAG_in_GEO=
               G4ThreeVector(0.,0.,1.).cross(geodipole_axis_in_GEO);
   
  YMAG_in_GEO/=YMAG_in_GEO.mag();
  G4ThreeVector XMAG_in_GEO= YMAG_in_GEO.cross(geodipole_axis_in_GEO);
  GEOtoMAG = G4RotationMatrix(XMAG_in_GEO,YMAG_in_GEO,geodipole_axis_in_GEO);
  
  // GEOtoMAG = G4RotationMatrix(0.,geodipole_axis_in_GEO.theta(),
  //                                90.*degree- geodipole_axis_in_GEO.phi());  
 
  // GEI and GSE
  // GSE X axis earth_sun direction
  // GSE Y axis in the ecliptic plane toward dusk 
  // GSE Z axis parallel to ecliptic pole  
 
  G4ThreeVector XGSE=G4ThreeVector(cos(sun_rasn)*cos(sun_dec),
                                      sin(sun_rasn)*cos(sun_dec),
				      sin(sun_dec));
				      
     
  G4cout<<"XGSE in GEI "<<XGSE<<std::endl;
  G4ThreeVector ZGSE=G4ThreeVector(0,-sin(obliq),cos(obliq));
  G4ThreeVector YGSE=ZGSE.cross(XGSE);
  GEItoGSE =G4RotationMatrix(XGSE,YGSE,ZGSE);
     
  //GEO to GSE   
  GEOtoGSE = GEItoGEO.inverse()*GEItoGSE;
    
  // GSM 
  // GSM X axis earth_sun direction
  // GSM Y axis perpendicular to the geomagnetic dipole axis
  // GSM xz plane contains geomagnetic dipole axis
  
  
  
  G4ThreeVector XGSM_in_GEI=G4ThreeVector(cos(sun_rasn)*cos(sun_dec),
                                      sin(sun_rasn)*cos(sun_dec),
				      sin(sun_dec));
  G4ThreeVector XGSM_in_GEO=GEItoGEO.inverse()*XGSM_in_GEI;
  G4ThreeVector YGSM_in_GEO= geodipole_axis_in_GEO.cross(XGSM_in_GEO);
  YGSM_in_GEO/=YGSM_in_GEO.mag();
  G4ThreeVector ZGSM_in_GEO= XGSM_in_GEO.cross(YGSM_in_GEO);
   
  GEOtoGSM=G4RotationMatrix(XGSM_in_GEO,YGSM_in_GEO,ZGSM_in_GEO);
  
  
  // Tilt angle 
  
  tilt_angle= geodipole_axis_in_GEO.angle(ZGSM_in_GEO);
  G4double angle1 = geodipole_axis_in_GEO.angle(XGSM_in_GEO);
  if ((angle1/degree) > 90.) tilt_angle = -tilt_angle; 
  G4cout<<geodipole_axis_in_GEO<<G4endl;
  G4cout<<tilt_angle/degree<<G4endl; 
   
   
  //SM
  //SM Z axis  geamagnetic dipole axis 
  //SM Y axis perpendicular to sun earth line = GSM Y axis
  //SM xz plane conatins sun earth line
   
  G4ThreeVector XSM_in_GEO=YGSM_in_GEO.cross(geodipole_axis_in_GEO);
  XSM_in_GEO/=XSM_in_GEO.mag();
  GEOtoSM=G4RotationMatrix(XSM_in_GEO,YGSM_in_GEO,geodipole_axis_in_GEO);
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::ComputeSelectedMatrix()
{
 G4RotationMatrix GEO_to_in = G4RotationMatrix();
 G4RotationMatrix out_to_GEO = G4RotationMatrix();
 SelectedMatrix = G4RotationMatrix(); 
 
 if (system_in == system_out) return; 
 
 if (system_in == "MAG") GEO_to_in =GEOtoMAG;
 else if  (system_in == "GSM") GEO_to_in =GEOtoGSM;
 else if  (system_in == "SM") GEO_to_in =GEOtoSM;
 else if  (system_in == "GSE") GEO_to_in =GEOtoGSE;
 else if  (system_in == "GEI") GEO_to_in =GEItoGEO.inverse();
 
 if (system_out == "MAG") out_to_GEO =GEOtoMAG.inverse();
 else if  (system_out == "GSM") out_to_GEO =GEOtoGSM.inverse();
 else if  (system_out == "SM") out_to_GEO =GEOtoSM.inverse();
 else if  (system_out == "GSE") out_to_GEO =GEOtoGSE.inverse();
 else if  (system_out == "GEI") out_to_GEO =GEItoGEO;
 
 
 SelectedMatrix=out_to_GEO*GEO_to_in;
 
}
/////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinateConvertor::Initialise()
{ComputeSunPosition();
 ComputeGeodipoleParameters1();
 ComputeTranformationMatrices();
 ComputeSelectedMatrix();
}







