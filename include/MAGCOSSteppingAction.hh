#ifndef MAGCOSSteppingAction_h
#define MAGCOSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class MAGCOSEventAction;

class MAGCOSSteppingAction : public G4UserSteppingAction
{
  public:
    MAGCOSSteppingAction();
    virtual ~MAGCOSSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);
    inline void SetStopAltitude(G4double alt){stop_altitude = alt;}
    
  private:
    MAGCOSEventAction* eventAction;
    G4double stop_altitude;
};

#endif
