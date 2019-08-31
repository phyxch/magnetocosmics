#ifndef MAGCOSTrackingAction_h
#define MAGCOSTrackingAction_h 1
#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "globals.hh"

////////////////////////////////////////////////////////////////////////////////
//
class MAGCOSTrackingAction : public G4UserTrackingAction
{
  public:
    MAGCOSTrackingAction ();
    ~MAGCOSTrackingAction ();
   
    void PreUserTrackingAction (const G4Track* theTrack);
    void PostUserTrackingAction (const G4Track* theTrack);

  private:
    G4bool RegisterLastPoint;
    
  public:
    inline void SetRegisterLastPoint(G4bool abool)
                                      {RegisterLastPoint=abool;}  
};
////////////////////////////////////////////////////////////////////////////////
#endif
