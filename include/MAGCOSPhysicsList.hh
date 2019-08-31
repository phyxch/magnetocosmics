#ifndef MAGCOSPhysicsList_h
#define MAGCOSPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class MAGCOSPhysicsList: public G4VUserPhysicsList
{
  public:
    MAGCOSPhysicsList();
    virtual ~MAGCOSPhysicsList(); 
    
  protected:
    virtual void ConstructParticle();
    virtual void ConstructProcess(); 
    virtual void SetCuts();
   
  protected: 
    virtual void ConstructBosons();
    virtual void ConstructLeptons();
    virtual void ConstructBarions();
    
   // virtual void ConstructGeneral();
  private:
    void AddMyTransportation();
    
    //MAGCOSPhysicsMessenger* myPhysicsMessenger;
    
};
#endif



