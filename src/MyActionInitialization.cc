#include "MyActionInitialization.hh"
#include "MyPrimaryGenerator.hh"
#include "MyEventAction.hh"
#include "MyStackingAction.hh"
#include "MySteppingAction.hh"

 MyActionInitialization:: MyActionInitialization()
{}

 MyActionInitialization::~MyActionInitialization()
{}

void MyActionInitialization::BuildForMaster() const
{
}

void  MyActionInitialization::Build() const
{
     MyPrimaryGenerator *generator = new MyPrimaryGenerator();
    SetUserAction(generator);

    MyEventAction *eventAction = new MyEventAction();
    SetUserAction(eventAction);

    // Optical photon counting and detection
    MyStackingAction *stackingAction = new MyStackingAction();
    SetUserAction(stackingAction);

    MySteppingAction *steppingAction = new MySteppingAction();
    SetUserAction(steppingAction);
}
