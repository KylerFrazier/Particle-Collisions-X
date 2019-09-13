/*
----------------------------------------------------------------------------------
Kyler Frazier | July 2019
    
Description: 
    This macro takes in a Delphes output file, and checks a variety of conditions
    to see if there is a possible canditate for the "X" particle. A histogram of
    the invarient mass of the "X" particle is outputted. 

Details - For each event in the Delphes output file, the macro does the following:
  * Check that the even has exactly 1 photon and at least 2 jets (at least 2 jets
    are needed to combine to 1 Z Bozon).
  * Get the 4-vectors of all jets.
  * Combine these 4-vectors in every possible combination in groupings of 2; e.g.
    with 3 jets, let them be named A, B, and C respectively, it will combine them
    in groupings of 2 such as ( {A},{B,C} ), ( {B},{A,C} ), and ( {C},{A,B} ).
  * Each pair of groupings is checked to see if the both of the 2 groupings can
    combine to form a Z-boson. If the invarient mass of the grouping is within
    10 GeV/c^2 of the accepted invarient mass of the Z-boson, then it continues.
  * If both groups in one pair could be Z-bosons, then they are combined and the
    invarient mass of this combination is put into the histogram. 

Running:
    root 
    .X Runner6.cpp("data.root");
----------------------------------------------------------------------------------
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include <algorithm>
#include <iostream>
#include <string>

// Function that generates all combinations C(N,r) and puts it into a vector
void comb(int N, int r, vector<vector<int>> &combList)
{
    string s(r, 1);
    s.resize(N, 0);
    do {
        vector<int> v;
        for (int i = 0; i < N; ++i)
        {
            if (s[i]) 
            {
                v.push_back(i);
            }
        }
        combList.push_back(v);
    } while (prev_permutation(s.begin(), s.end()));
}

// Main
void Runner6(const char *inputFile)
{
    gSystem->Load("libDelphes");

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    TH1F *hist = new TH1F("hist", "Particle \"X\" Invarient Mass Histogram", 150, 0, 3000);
    hist -> GetXaxis() -> SetTitle("Invarient Mass (GeV/c^{2})");
    hist -> GetYaxis() -> SetTitle("Instances");

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");

    Photon *photon;
    Track *track;
    Tower *tower;

    Jet *jet;
    TObject *object;

    Long64_t entry;

    Int_t i, j, k;

    // Can change the margin to accept Z-boson masses below
    // Can change the sample size of the data to loop over for debugging
    double ZMass = 91.1876;
    double margin = 10.0;
    double sample_percent = 100;
    Long64_t allEntries = sample_percent*(treeReader->GetEntries())/100;
    
    int proper_events = 0;
    int jetN = 0;
    
    cout << "  * Chain contains " << allEntries << " events" << endl;
    cout << "  * Margin for Z Mass: " << margin << endl;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {        
        // Progress bar
        if( entry % (allEntries/100) == 0)
        {
            cout << "Progress: " << 100.0*entry/allEntries << "%" << endl;
        }
        
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        // Skip event conditions
        if (branchPhoton-> GetEntriesFast() != 1 ||
            branchJet->    GetEntriesFast() <  2  )
        {
            continue;
        }
        
        // jetN is the number of jets per event, jetsLV is an array of Lorentz Vectors
        jetN = branchJet->GetEntriesFast();
        TLorentzVector jetsLV[jetN];

        // Loop over all jets in event and assigning Lorentz Vectors
        for(i = 0; i < branchJet->GetEntriesFast(); ++i)
        {
            jet = (Jet*) branchJet->At(i);
            jetsLV[i] = jet->P4();
        }

        // Creates all possible combinations and groupings and puts them in combList
        vector<vector<int>> combList;
        for(i = 1;i <= jetN/2;++i)
        {
            comb(jetN,i,combList);
        }

        // Future implementation:
        //     Reorder combList such that 2 | 2, 1 | 2, and 1 | 3 appear first in that order.
        //     Can use if statement to search for each comb, checking for existance, then pop and re-add at index

        // Loops over all combinations of jets and checks if they add to Z-bosons
        for (i = 0; i < combList.size(); ++i)
        {
            TLorentzVector Z1,Z2,X;
            for (j = 0; j < combList[i].size(); ++j) 
            {
                //cout << combList[i][j];
                Z1 += jetsLV[combList[i][j]];
            }
            //cout << " | ";
            for (j = 0; j < jetN; ++j) 
            {
                if(find(combList[i].begin(), combList[i].end(), j) == combList[i].end()) 
                {
                    //cout << j;
                    Z2 += jetsLV[j];
                }
            }
            if((ZMass-margin) < Z1.M() && Z1.M() < (ZMass+margin) && 
               (ZMass-margin) < Z2.M() && Z2.M() < (ZMass+margin) )
            {
                //cout << "\n\nZ1 Mass: " << Z1.M() << "\nZ2 Mass: " << Z2.M() << "\n\n" << endl;
                X += Z1 + Z2;
                hist -> Fill(X.M());
                proper_events += 1;
                //break;
            } 
        }
    }

//------------------------------------------------------------------------------

    cout << "Progress: 100%" << endl;
    cout << "Number of events that match criteria: " << proper_events << endl;
    cout << "Percentage of such events:            " << 100.0*proper_events/allEntries << "%" << endl;
    hist -> Draw();

    cout << "  * Exiting..." << endl;

    delete treeReader;
    delete chain;
}