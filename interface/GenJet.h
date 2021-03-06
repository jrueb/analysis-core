#ifndef Analysis_Core_GenJet_h
#define Analysis_Core_GenJet_h 1

// -*- C++ -*-
//
// Package:    Analysis/Core
// Class:      GenJet
// 
/**\class GenJet GenJet.cc Analysis/Core/src/GenJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <memory>
// 
// user include files
#include "Analysis/Core/interface/Candidate.h"
//
// class declaration
//

namespace analysis {
   namespace tools {

      class GenJet : public Candidate {
         public:
            GenJet();
            GenJet(const float & pt, const float & eta, const float & phi, const float & e, const float & q);
           ~GenJet();
//            using Candidate::set; // in case needed to overload the function set
      
         private:
            // ----------member data ---------------------------
            
            // 
      };
   }
}

#endif  // Analysis_Core_GenJet_h
