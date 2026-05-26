/**
 * @file MyEventAction.cc
 * @brief Manages the start and end of individual events.
 * 
 * WHAT THIS CLASS DOES:
 * An "Event" in Geant4 is a single particle shower (e.g., one electron hitting 
 * the detector). This class is responsible for data management per event.
 * At the start of an event, it resets all the data variables. At the end of 
 * the event, it gathers all the counted photons from the Stacking and Stepping 
 * actions and writes them out into our ROOT output file so we can analyze 
 * the detector's optical performance later.
 */

#include "MyEventAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "MyHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4ios.hh"
#include "MySensitiveDetector.hh"
#include "MyStackingAction.hh"
#include "MySteppingAction.hh"

MyEventAction::MyEventAction()
    : tree(nullptr), treeTracks(nullptr), treeOptical(nullptr),
      treeOptPhotons(nullptr), outputFile(nullptr), initialized(false)
{
    G4int tid = G4Threading::G4GetThreadId();

    outputFile = new TFile(Form("output/output_thread%d.root", tid), "RECREATE");

    tree = new TTree("Hits", "Hit Data Tree");

    tree->Branch("fEvent", &fEvent, "fEvent/I");
    tree->Branch("fEdep", &fEdep, "fEdep/D");
    tree->Branch("fScint", &fScint, "fScint/D");
    tree->Branch("fCher",  &fCher,  "fCher/D");
    tree->Branch("fDetectorID", &fDetectorID, "fDetectorID/I");
    tree->Branch("fTrackID", &fTrackID, "fTrackID/I");

    treeTracks = new TTree("Tracks", "Initial track kinetic energy per event");

    treeTracks->Branch("fEvent", &fEvent, "fEvent/I");
    treeTracks->Branch("fTrackID", &bTrackID, "fTrackID/I");
    treeTracks->Branch("fInitKE", &bInitKE, "fInitKE/D");

    treeOptical = new TTree("Optical", "Optical photon summary per event");

    treeOptical->Branch("fEvent",           &fEvent,            "fEvent/I");
    treeOptical->Branch("nCherenkov",       &oNcherenkov,       "nCherenkov/I");
    treeOptical->Branch("nScintillation",   &oNscintillation,   "nScintillation/I");
    treeOptical->Branch("nDetected",        &oNdetected,        "nDetected/I");
    treeOptical->Branch("nDetCherenkov",    &oNdetCherenkov,    "nDetCherenkov/I");
    treeOptical->Branch("nDetScintillation",&oNdetScintillation, "nDetScintillation/I");

    treeOptPhotons = new TTree("OpticalPhotons", "Detected optical photon data");

    treeOptPhotons->Branch("fEvent",       &fEvent,      "fEvent/I");
    treeOptPhotons->Branch("wavelength",   &pWavelength, "wavelength/D");
    treeOptPhotons->Branch("time",         &pTime,       "time/D");
    treeOptPhotons->Branch("posX",         &pPosX,       "posX/D");
    treeOptPhotons->Branch("posY",         &pPosY,       "posY/D");
    treeOptPhotons->Branch("posZ",         &pPosZ,       "posZ/D");
    treeOptPhotons->Branch("creatorProcess",&pCreatorType,"creatorProcess/I");
    treeOptPhotons->Branch("detectorID",   &pDetectorID, "detectorID/I");

    initialized = true;
}

MyEventAction::~MyEventAction()
{
    if (outputFile) {
        outputFile->cd();
        if (tree) tree->Write();
        if (treeTracks) treeTracks->Write();
        if (treeOptical) treeOptical->Write();
        if (treeOptPhotons) treeOptPhotons->Write();
        outputFile->Write();
        outputFile->Close();
        delete outputFile;
        outputFile = nullptr;
    }
}

void MyEventAction::BeginOfEventAction(const G4Event* event)
{
    fEvent = event->GetEventID();

    MySteppingAction* stepping = MySteppingAction::Instance();
    if (stepping) stepping->Reset();
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
    fEvent = event->GetEventID();

    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    if (sdManager) {
        G4int hcID = sdManager->GetCollectionID("MyHitCollection");
        if (hcID >= 0) {
            G4HCofThisEvent* hce = event->GetHCofThisEvent();
            if (hce) {
                MyHitCollection* hitsCollection =
                    static_cast<MyHitCollection*>(hce->GetHC(hcID));

                if (hitsCollection) {
                    for (size_t i = 0; i < hitsCollection->GetSize(); ++i) {
                        MyHit* hit = static_cast<MyHit*>((*hitsCollection)[i]);
                        if (!hit) continue;

                        fEdep = hit->GetEdep();
                        fScint = hit->GetScint();
                        fCher  = hit->GetCher();

                        fDetectorID = hit->GetID();
                        fTrackID = hit->GetTrackID();

                        tree->Fill();
                    }
                }
            }
        }
    }

    MySensitiveDetector* sd = MySensitiveDetector::Instance();
    if (sd) {
        const auto &trackMap = sd->GetTrackInitKEMap();
        for (const auto &kv : trackMap) {
            bTrackID = kv.first;
            bInitKE  = kv.second;
            fEvent   = event->GetEventID();
            treeTracks->Fill();
        }
    }

    MyStackingAction* stacking = MyStackingAction::Instance();
    MySteppingAction* stepping = MySteppingAction::Instance();

    oNcherenkov       = stacking ? stacking->GetCherenkovCount()     : 0;
    oNscintillation   = stacking ? stacking->GetScintillationCount() : 0;
    oNdetected        = stepping ? stepping->GetDetTotalCount()          : 0;
    oNdetCherenkov    = stepping ? stepping->GetDetCherenkovCount()      : 0;
    oNdetScintillation= stepping ? stepping->GetDetScintillationCount()  : 0;

    treeOptical->Fill();

    if (stepping) {
        const auto& photons = stepping->GetDetectedPhotons();
        for (const auto& ph : photons) {
            pWavelength  = ph.wavelength;
            pTime        = ph.time;
            pPosX        = ph.posX;
            pPosY        = ph.posY;
            pPosZ        = ph.posZ;
            pCreatorType = ph.creatorType;
            pDetectorID  = ph.detectorID;
            treeOptPhotons->Fill();
        }
    }
}
