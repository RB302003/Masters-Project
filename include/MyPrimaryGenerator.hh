#ifndef MYPRIMARYGENERATOR_HH
#define MYPRIMARYGENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator();

    virtual void GeneratePrimaries(G4Event*);

private:
    G4ParticleGun* fParticleGun;
};

#endif
