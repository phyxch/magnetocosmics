#ifndef MAGCOSVisManager_h
#define MAGCOSVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MAGCOSVisManager: public G4VisManager {

public:

  MAGCOSVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
