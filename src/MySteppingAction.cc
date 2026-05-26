/**
 * @file MySteppingAction.cc
 * @brief Tracks the step-by-step movement of particles.
 * 
 * WHAT THIS CLASS DOES:
 * This class observes particles as they move through the detector geometry.
 * We use it to detect when an optical photon reaches the back face of our 
 * PbWO4 crystal (representing a photodetector sensor). When a photon hits 
 * this boundary, the class records its properties (wavelength, arrival time, 
 * and whether it was a Cherenkov or Scintillation photon), saves this data 
 * to the ROOT file, and then stops tracking the photon.
 */

#include "MySteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

thread_local MySteppingAction* MySteppingAction::sInstance = nullptr;

MySteppingAction::MySteppingAction()
    : G4UserSteppingAction(),
      fDetCherenkov(0),
      fDetScintillation(0),
      fDetTotal(0)
{
    sInstance = this;
}

MySteppingAction::~MySteppingAction()
{
    if (sInstance == this) sInstance = nullptr;
}

void MySteppingAction::Reset()
{
    fDetectedPhotons.clear();
    fDetCherenkov = 0;
    fDetScintillation = 0;
    fDetTotal = 0;
}

void MySteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return;

    G4StepPoint* preStep  = step->GetPreStepPoint();
    G4StepPoint* postStep = step->GetPostStepPoint();

    G4VPhysicalVolume* preVol = preStep->GetPhysicalVolume();
    if (!preVol) return;
    G4String preVolName = preVol->GetName();
    if (preVolName != "Detector") return;

    G4VPhysicalVolume* postVol = postStep->GetPhysicalVolume();
    if (!postVol) return;

    G4String postVolName = postVol->GetName();

    if (postVolName == preVolName) return;

    G4ThreeVector pos = postStep->GetPosition();
    G4ThreeVector mom = track->GetMomentumDirection();

    if (mom.z() > 0.0) {
        DetectedPhoton photon;

        G4double energy = track->GetKineticEnergy();
        photon.wavelength = (h_Planck * c_light / energy) / nm;
        photon.time = postStep->GetGlobalTime() / ns;
        photon.posX = pos.x() / mm;
        photon.posY = pos.y() / mm;
        photon.posZ = pos.z() / mm;

        const G4VProcess* creator = track->GetCreatorProcess();
        if (creator) {
            G4String pname = creator->GetProcessName();
            if (pname == "Cerenkov") {
                photon.creatorType = 0;
                fDetCherenkov++;
            } else if (pname == "Scintillation") {
                photon.creatorType = 1;
                fDetScintillation++;
            } else {
                photon.creatorType = 2;
            }
        } else {
            photon.creatorType = 2;
        }

        photon.detectorID = preStep->GetTouchable()->GetCopyNumber();

        fDetectedPhotons.push_back(photon);
        fDetTotal++;

        track->SetTrackStatus(fStopAndKill);
    }
}
