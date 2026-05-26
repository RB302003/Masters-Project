#ifndef MYEVENTACTION_HH
#define MYEVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include <TFile.h>
#include <TTree.h>

class MyEventAction : public G4UserEventAction {
public:
    MyEventAction();
    virtual ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*) override;
    virtual void EndOfEventAction(const G4Event*) override;

private:
    // Hits tree branches
    G4int fTrackID = -1;
    G4double fEdep = 0.0;
    G4double fScint = 0.0;
    G4double fCher  = 0.0;
    G4int fEvent = -1;
    G4int fDetectorID = -1;

    // Tracks tree branches
    G4int bTrackID = -1;
    G4double bInitKE = 0.0;

    // Optical summary tree branches (per-event)
    G4int oNcherenkov = 0;
    G4int oNscintillation = 0;
    G4int oNdetected = 0;
    G4int oNdetCherenkov = 0;
    G4int oNdetScintillation = 0;

    // Optical photons tree branches (per-detected-photon)
    G4double pWavelength = 0.0;
    G4double pTime = 0.0;
    G4double pPosX = 0.0, pPosY = 0.0, pPosZ = 0.0;
    G4int pCreatorType = 0;
    G4int pDetectorID = 0;

    // ROOT IO
    TTree *tree = nullptr;
    TTree *treeTracks = nullptr;
    TTree *treeOptical = nullptr;
    TTree *treeOptPhotons = nullptr;
    TFile *outputFile = nullptr;
    bool initialized = false;
};

#endif
