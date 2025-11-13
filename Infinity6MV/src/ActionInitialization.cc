#include "ActionInitialization.hh"
MyActionInitialization::MyActionInitialization(): G4VUserActionInitialization(){}
MyActionInitialization::~MyActionInitialization(){}



void MyActionInitialization::Build() const{

    MyPrimaryGenerator *generator =new MyPrimaryGenerator();
    SetUserAction(generator);
    MyRunAction *runAction= new MyRunAction();
    SetUserAction(runAction);
    MyEventAction *eventAction = new MyEventAction( );
    SetUserAction(eventAction);
    MySteppingAction *steppingAction = new MySteppingAction(eventAction);
    SetUserAction(steppingAction);
    TrackingAction *trackingAction = new TrackingAction();
SetUserAction(trackingAction);
}
void MyActionInitialization::BuildForMaster() const{
   
}