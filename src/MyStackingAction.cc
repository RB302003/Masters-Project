/**
 * @file MyStackingAction.cc
 * @brief Handles the tracking of newly generated particles.
 * 
 * WHAT THIS CLASS DOES:
 * This class acts as a filter for every new particle created in the simulation.
 * We use it to specifically look for new optical photons. When a photon is created,
 * this class checks if it was created by "Cerenkov" (Cherenkov radiation) or 
 * "Scintillation". It then increments the respective counters so we know exactly 
 * how many photons of each type were produced in the event.
 */

#include "MyStackingAction.hh"

#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

thread_local MyStackingAction* MyStackingAction::sInstance = nullptr;

MyStackingAction::MyStackingAction()
    : G4UserStackingAction(),
      fCherenkovCount(0),
      fScintillationCount(0)
{
    sInstance = this;
}

MyStackingAction::~MyStackingAction()
{
    if (sInstance == this) sInstance = nullptr;
}

G4ClassificationOfNewTrack
MyStackingAction::ClassifyNewTrack(const G4Track* track)
{
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        const G4VProcess* creator = track->GetCreatorProcess();
        if (creator) {
            G4String processName = creator->GetProcessName();
            if (processName == "Cerenkov") {
                fCherenkovCount++;
            } else if (processName == "Scintillation") {
                fScintillationCount++;
            }
        }
    }

    return fUrgent;
}

void MyStackingAction::NewStage()
{
}

void MyStackingAction::PrepareNewEvent()
{
    fCherenkovCount = 0;
    fScintillationCount = 0;
}
