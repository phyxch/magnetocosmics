#include "DurationManager.hh"
#include "G4Timer.hh"
#include "MAGCOSUnits.hh"
#include "DurationMessenger.hh"


DurationManager* DurationManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
DurationManager::DurationManager()  
{ 
  TotalDuration =0.;
  RunDuration =0.;
  EventDuration =0.;
  MaxTotalDuration =1.e300;
  MaxRunDuration =1.e300;
  MaxEventDuration =1.e300;
  theTimer=new G4Timer;
  theTimer->Start();	
  theDurationMessenger = new DurationMessenger(this);
}
////////////////////////////////////////////////////////////////////////////////
//
DurationManager::~DurationManager() 
{ 
  if (theTimer) delete  theTimer;
}
////////////////////////////////////////////////////////////////////////////////
//
DurationManager* DurationManager::getInstance()
{
  if (instance == 0) instance = new DurationManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
bool DurationManager::CheckDurationAtBeginOfEvent()
{
  ComputeDuration();
  EventDuration = 0.;
  return CheckDuration();
}
////////////////////////////////////////////////////////////////////////////////
//
bool DurationManager::CheckDurationAtBeginOfRun()
{ 
  ComputeDuration();
  RunDuration = 0.;
  EventDuration = 0.;
  return CheckDuration();
}

////////////////////////////////////////////////////////////////////////////////
//
void DurationManager::SetMaxRunDuration(G4double time)
{ 
  MaxRunDuration=time;
}
////////////////////////////////////////////////////////////////////////////////
//
void DurationManager::SetMaxEventDuration(G4double time)
{ 
  MaxEventDuration=time;
}
////////////////////////////////////////////////////////////////////////////////
//
void DurationManager::SetMaxTotalDuration(G4double time)
{ 
  MaxTotalDuration=time;
}
////////////////////////////////////////////////////////////////////////////////
//
bool DurationManager::CheckDuration()
{ 
  ComputeDuration();
  if (RunDuration >= MaxRunDuration || 
      EventDuration >= MaxEventDuration ||
      TotalDuration >= MaxTotalDuration) {
    
    if (  TotalDuration >= MaxTotalDuration) {
      G4cout<<TotalDuration<<" s total time is higher than the maximum "
	    <<"allowed total duration of "<< MaxTotalDuration<<" s"<<std::endl;
    }
    if (  RunDuration >= MaxRunDuration) {
      G4cout<<RunDuration<<" s run time is higher than the maximum "
	    <<"allowed run duration of "<< MaxRunDuration<<" s"<<std::endl;
    }
    if (  EventDuration >= MaxEventDuration) {
      G4cout<<EventDuration<<" s event time is higher than the maximum "
	    <<"allowed event duration of "<< MaxEventDuration<<" s"<<std::endl;
    }	
    
    
    return false;
  }
  return true;    
  
}
////////////////////////////////////////////////////////////////////////////////
//
void DurationManager::ComputeDuration()
{ 
  theTimer->Stop();
  RunDuration += theTimer->GetRealElapsed();
  EventDuration += theTimer->GetRealElapsed();
  TotalDuration += theTimer->GetRealElapsed();
  theTimer->Start();
}
