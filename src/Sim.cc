#include <iostream>
#include <cstdlib>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "MyDetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "MyActionInitialization.hh"

int main(int argc, char** argv)
{
    // Create output directory once (not per-thread)
    (void)system("mkdir -p output");

    G4RunManager *runManager = new G4RunManager();

    MyDetectorConstruction* detector = new MyDetectorConstruction();
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new MyPhysicsList());
    runManager->SetUserInitialization(new MyActionInitialization());
    // NOTE: Do NOT call runManager->Initialize() here.
    // The macro file must set /detector/mode BEFORE calling /run/initialize.

    G4UIExecutive *ui = nullptr;
    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (ui)
    {
        // Interactive mode: initialize visualization
        G4VisManager *visManager = new G4VisExecutive();
        visManager->Initialize();

        UImanager->ApplyCommand("/control/execute run.mac");
        ui->SessionStart();

        delete visManager;
        delete ui;
    }
    else
    {
        // Batch mode: no visualization needed
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    delete runManager;
    return 0;
}
