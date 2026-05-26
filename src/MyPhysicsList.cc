#include "MyPhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"

MyPhysicsList::MyPhysicsList()
{
    RegisterPhysics(new G4EmStandardPhysics());
    RegisterPhysics(new G4HadronElasticPhysics());
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
    RegisterPhysics(new G4DecayPhysics());

    // Optical physics: Cherenkov, Scintillation, Absorption, Rayleigh, Boundary
    RegisterPhysics(new G4OpticalPhysics());

    // Configure optical parameters via the singleton
    G4OpticalParameters* optParams = G4OpticalParameters::Instance();
    optParams->SetCerenkovMaxPhotonsPerStep(20);
    optParams->SetCerenkovTrackSecondariesFirst(true);
    optParams->SetScintTrackSecondariesFirst(true);
}

MyPhysicsList::~MyPhysicsList()
{}
