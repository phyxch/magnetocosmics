#ifndef DurationManager_HH
#define DurationManage_HH 1

#include"globals.hh"

class G4Timer;
class DurationMessenger;
class DurationManager
{
public:

   ~DurationManager(); 
   static DurationManager* getInstance();
   void SetMaxRunDuration(G4double time);
   void SetMaxEventDuration(G4double time);
   void SetMaxTotalDuration(G4double time);
   bool CheckDurationAtBeginOfEvent();
   bool CheckDurationAtBeginOfRun();
   bool CheckDuration();
   void ComputeDuration();
   
   
 
   
   
private:
   static DurationManager* instance;
   
   DurationManager();  
   G4Timer* theTimer;
   G4double RunDuration,TotalDuration,EventDuration;
   G4double MaxRunDuration,MaxTotalDuration,MaxEventDuration;
   DurationMessenger* theDurationMessenger;
   
   
 
};

#endif




