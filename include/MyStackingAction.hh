#ifndef MYSTACKINGACTION_HH
#define MYSTACKINGACTION_HH

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "globals.hh"

class G4Track;

/**
 * MyStackingAction:
 * Counts optical photons as they are generated (Cherenkov vs Scintillation).
 * Provides per-event counters that EventAction reads at end of event.
 */
class MyStackingAction : public G4UserStackingAction {
public:
    MyStackingAction();
    ~MyStackingAction() override;

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
    void NewStage() override;
    void PrepareNewEvent() override;

    // Accessors for EventAction
    G4int GetCherenkovCount()     const { return fCherenkovCount; }
    G4int GetScintillationCount() const { return fScintillationCount; }

    // Thread-local singleton for EventAction access
    static MyStackingAction* Instance() { return sInstance; }

private:
    G4int fCherenkovCount;
    G4int fScintillationCount;

    static thread_local MyStackingAction* sInstance;
};

#endif
