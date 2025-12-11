#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh" 
#include "PrimaryGenerator.hh"
#include "run.hh"
#include "event.hh"
#include "stepping.hh"
#include "TrackingAction.hh"
class MyActionInitialization : public G4VUserActionInitialization{

    public:
        MyActionInitialization();
        ~MyActionInitialization();
        virtual void Build() const;
        virtual void BuildForMaster() const;
    private:
    
    };

#endif