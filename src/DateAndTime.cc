#include "DateAndTime.hh"
#include "G4ios.hh"
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime()
{ year=2000;
  month=1;
  day=1;
  hour=0;
  min=0;
  sec=0;
  msec=0;
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d)
{ year=y;
  month=mo;
  day=d;
  hour=0;
  min=0;
  sec=0;
  msec=0;
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d, int h, int mi, int s)
{year=y;
 month=mo;
 day=d;
 hour=h;
 min=mi;
 sec=s;
 msec=0;
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d, int h, int mi, int s, int ms)
{ year=y;
  month=mo;
  day=d;
  hour=h;
  min=mi;
  sec=s;
  msec=ms;
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(double julian_date)
{ ConvertToJulianDate(julian_date);
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::~DateAndTime()
{;
}
////////////////////////////////////////////////////////////////////////////////
//	
double DateAndTime::DifferenceInDays(const DateAndTime aDate) const
{ return JulianDate()-aDate.JulianDate();
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInHours(const DateAndTime aDate) const
{ return DifferenceInDays(aDate) * 24.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInMinutes(const DateAndTime aDate) const
{ return DifferenceInDays(aDate) * 24. * 60.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInSeconds(const DateAndTime aDate) const
{ return DifferenceInDays(aDate) * 24. * 3600.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInMilliseconds(const DateAndTime aDate)const
{ return DifferenceInDays(aDate) * 24. *  3600000.;
}
////////////////////////////////////////////////////////////////////////////////
//
long DateAndTime::JulianDay() const
{// The julian day 0 is 1 January 4713 BC    
 // Algorithm for computing julian day from Year Month and Day of Gregorian
 // calendar, taken from 
 //"explanatory supplement to the Astronomical Almanach", 1992, p.604 
 
  long ly=long(year);
  long lm=long(month);
  long ld=long(day);
  
  // Adjust BC years
  if ( ly < 0 ) ly++;
  
  long JD= ld - 32075L +
           1461L * (ly + 4800L + (lm - 14L) / 12L) / 4L +
           367L * (lm - 2L - 12L * ((lm - 14L) / 12L)) / 12L -
           3L * ((ly + 4900L + (lm - 14L) / 12L) / 100L) / 4L;

  return JD;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::JulianDate() const
{// The julian date 0.0 starts at 1 January 4713 BC at Greenwhich mean noon   
 
 
 double dms=double(msec);
 double dsec=double(sec);
 double dmin=double(min);
 double dhour=double(hour);
 
 double fraction_day = 
          (((dms / 1000. + dsec) / 60. + dmin) / 60. + dhour) / 24. - 0.5;
  	  
 double JulDate = double (JulianDay()) + fraction_day;
 return JulDate;
}
////////////////////////////////////////////////////////////////////////////////
//
int DateAndTime::DayOfYear() const
{// return the day of the year with 1st january =day 0
   int day_of_year = int (DifferenceInDays(DateAndTime(year,1,1))); 
   return day_of_year;  
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::ConvertToJulianDate(double aJulianDate)
{// set the date according to a given julian date

 //Julian day
 // 0.5 is because Julian Date are calciulated form noon on 1st January BC 4713
  
  long JD=long(aJulianDate + 0.5);
  double fraction_of_day =  (aJulianDate + 0.5) - double(JD);
  
  //calculating year month day with algorithm from
  //"explanatory supplement to the Astronomical Almanach", 1992, p.604 
   
  long L = JD + 68569L;
  long N = (4 * L) / 146097L;
  L = L - ((146097 * N + 3) / 4L);
  long I = (4000L * (L + 1))/1461001L;
  L = L - (1461L * I) / 4L + 31L;
  long J = (80L * L) / 2447L;
  day = int ( L - (2447L * J) / 80L );
  L  = J / 11L;
  month = int (J + 2 - (12 * L));
  year = int (100L * (N-49) + I + L);
   
  //calculating hour, min, sec , msec
  
  double dhour  = fraction_of_day * 24.;
  hour = int (dhour);
  double dmin = (dhour - double (hour)) * 60.;
  min = int (dmin);
  double dsec = (dmin - double(min)) * 60.;
  sec = int(dsec);
  msec = int ((dsec - double (sec)) * 1000.); 
    
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::ConvertToSpenvisModifiedJulianDay(double aModJulianDay)
{ ConvertToJulianDate(aModJulianDay + 2433282.5);    
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::operator=(const DateAndTime aDate)
{ year =aDate.year;
  //G4cout<<"date "<<aDate.year<<std::endl;
  month=aDate.month;
  day=aDate.day;
  hour=aDate.hour;
  min=aDate.min;
  sec=aDate.sec;
  msec=aDate.msec;
}
