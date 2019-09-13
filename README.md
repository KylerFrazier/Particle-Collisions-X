Kyler Frazier | July 2019
    
Description: 
    This macro takes in a Delphes output file and checks a variety of conditions
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