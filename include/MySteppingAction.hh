#ifndef MYSTEPPINGACTION_HH
#define MYSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

class G4Step;

/**
 * MySteppingAction:
 * Detects optical photons that exit from the rear face of PbWO4 crystals.
 * Records per-photon data: wavelength, time, position, creator process, crystal ID.
 */

struct DetectedPhoton {
    G4double wavelength;  // nm
    G4double time;        // ns
    G4double posX, posY, posZ;
    G4int    creatorType; // 0 = Cherenkov, 1 = Scintillation, 2 = other
    G4int    detectorID;  // copy number of crystal
};

class MySteppingAction : public G4UserSteppingAction {
public:
    MySteppingAction();
    ~MySteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;

    // Accessors for EventAction
    const std::vector<DetectedPhoton>& GetDetectedPhotons() const { return fDetectedPhotons; }

    G4int GetDetCherenkovCount()     const { return fDetCherenkov; }
    G4int GetDetScintillationCount() const { return fDetScintillation; }
    G4int GetDetTotalCount()         const { return fDetTotal; }

    void Reset();

    // Thread-local singleton
    static MySteppingAction* Instance() { return sInstance; }

private:
    std::vector<DetectedPhoton> fDetectedPhotons;
    G4int fDetCherenkov;
    G4int fDetScintillation;
    G4int fDetTotal;

    static thread_local MySteppingAction* sInstance;
};

#endif
