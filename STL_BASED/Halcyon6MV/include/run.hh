#ifndef RUN_HH
#define RUN_HH 
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh" 
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"  
class MyRunAction : public G4UserRunAction
{
public:
    MyRunAction();
    ~MyRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
private: 
};

#endif