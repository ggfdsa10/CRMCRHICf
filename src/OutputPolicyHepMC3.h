#ifndef _OutputPolicyHepMC3_h_
#define _OutputPolicyHepMC3_h_

#include <ctime>
#include <iostream>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include "OutputPolicyNone.h"
#include "CRMChepevt.h"
#include "CRMChepmc3.h"
#include "CRMCstat.h"

#include "TRandom3.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH2Poly.h"

using namespace std;

class CRMCoptions;

class OutputPolicyHepMC3 : public OutputPolicyNone
{
    enum runtypeIndex{
        kTL = 0,
        kTS = 1, 
        kTOP = 2,
        kALL = 3
    };

    public:
        OutputPolicyHepMC3(){}
        void InitOutput(const CRMCoptions& cfg) override;
        void FillRHICfEvent(const CRMCoptions& cfg, const int nEvent, int& passEventNum) override;
        void CloseOutput(const CRMCoptions& cfg) override;

    private:
        void PrintEvent();
        void InitVertexFluctuation();
        void InitRHICfGeometry();
        bool IsNeutralParticle(int pid);
        int GetRHICfGeoHit(double posX, double posY, double posZ, double px, double py, double pz, double e);

        CRMChepevt<HepMC3::GenParticlePtr,
            HepMC3::GenVertexPtr,
            HepMC3::FourVector,
            HepMC3::GenEvent> _hepevt;
        CRMChepmc3 _hepmc3;
        HepMC3::GenEvent _event;

        TFile* fFile;
        TTree* fEventTree;
        TTree* fRunTree;
        TClonesArray* fParticleArray;
        TParticle* fParticle;
        Int_t fRHICfRunType;
        Int_t fModelIdx;
        Int_t fProcessID;

        // ====== vertex fluctuation parameters =======
        TRandom3* fRandom;
        double fVertexMean[3]; // mm [x, y, z]
        double fVertexSigma[3]; // mm [x, y, z]

        // ======== RHICf Geometry =======
        TH2Poly* fRHICfPoly; // only west
        const double fRHICfDetZ = 17800.; // [mm]

    protected:

};


#endif

