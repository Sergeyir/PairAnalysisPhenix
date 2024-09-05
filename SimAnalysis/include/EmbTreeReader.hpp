// $HEADER$
//------------------------------------------------------------------------------------------------
//                     EmbTreeReader and EventReader classes declaration
//------------------------------------------------------------------------------------------------
// EmbTreeReader - embedding tree reader
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic classes for simple post simulation embedding TTree reading and for similar interaction as PHCentralTrack for variables retrieval for embedding efficiencies evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef EMB_TREE_READER_HPP
#define EMB_TREE_READER_HPP

#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

class EventReader
{
   public :

   EventReader(TTreeReader &reader);
   
   //Event variables
   Short_t type();
   Int_t geant_id();
   Double_t cent();
   Double_t bbcz();

   ~EventReader();

   private:
   
   //Event variables branches
   std::unique_ptr<TTreeReaderValue<short>> b_type;
   std::unique_ptr<TTreeReaderValue<int>> b_geant_id;
   std::unique_ptr<TTreeReaderValue<float>> b_cent;
   std::unique_ptr<TTreeReaderValue<float>> b_bbcz;
   
};

class EmbTreeReader
{
   public :

   EmbTreeReader(TTreeReader &reader, const std::string input_type);
   
   //DC-PC1 variables
   Double_t zed();
   Double_t the0();
   Double_t mom();
   Double_t phi();
   Double_t alpha();
   Short_t qual();

   //PC1 variables
   Double_t ppc1x();
   Double_t ppc1y();
   Double_t ppc1z();

   //PC2 variables
   Double_t ppc2x();
   Double_t ppc2y();
   Double_t ppc2z();
   Double_t pc2dz();
   Double_t pc2dphi();
   Double_t pc2sdz();
   Double_t pc2sdphi();
   
   //PC3 variables
   Double_t ppc3x();
   Double_t ppc3y();
   Double_t ppc3z();
   Double_t pc3dz();
   Double_t pc3dphi();
   Double_t pc3sdz();
   Double_t pc3sdphi();
   
   //EMCal variables
   Double_t emcsdz();
   Double_t emcsdphi();
   Double_t pemcy();
   Double_t pemcz();
   Double_t ecore();

   //TOFe variables
   Double_t tofsdz();
   Double_t tofsdphi();
   Double_t tofe();
   Double_t ptofy();
   Double_t ptofz();
   //Double_t slat(); //there is no slat variable in embedding output :(
   Double_t pltof();
   Double_t ttof();
   
   //TOFw variables
   Double_t tofwdz();
   Double_t tofwdphi();
   Int_t striptofw();
   Double_t pltofw();

   ~EmbTreeReader();

   private:
   
   //DC-PC1 variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_the0;
   std::unique_ptr<TTreeReaderValue<float>> b_zed;
   std::unique_ptr<TTreeReaderValue<float>> b_mom;
   std::unique_ptr<TTreeReaderValue<float>> b_phi;
   std::unique_ptr<TTreeReaderValue<float>> b_alpha;
   std::unique_ptr<TTreeReaderValue<short>> b_qual;

   //PC1 variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_ppc1x;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc1y;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc1z;
   
   //PC2 variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_ppc2x;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc2y;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc2z;
   std::unique_ptr<TTreeReaderValue<float>> b_pc2dz;
   std::unique_ptr<TTreeReaderValue<float>> b_pc2dphi;
   std::unique_ptr<TTreeReaderValue<float>> b_pc2sdz;
   std::unique_ptr<TTreeReaderValue<float>> b_pc2sdphi;
   
   //PC3 variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_ppc3x;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc3y;
   std::unique_ptr<TTreeReaderValue<float>> b_ppc3z;
   std::unique_ptr<TTreeReaderValue<float>> b_pc3dz;
   std::unique_ptr<TTreeReaderValue<float>> b_pc3dphi;
   std::unique_ptr<TTreeReaderValue<float>> b_pc3sdz;
   std::unique_ptr<TTreeReaderValue<float>> b_pc3sdphi;
   
   //EMCal variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_emcsdz;
   std::unique_ptr<TTreeReaderValue<float>> b_emcsdphi;
   std::unique_ptr<TTreeReaderValue<float>> b_pemcy;
   std::unique_ptr<TTreeReaderValue<float>> b_pemcz;
   std::unique_ptr<TTreeReaderValue<float>> b_ecore;

   //TOFe variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_tofsdz;
   std::unique_ptr<TTreeReaderValue<float>> b_tofsdphi;
   std::unique_ptr<TTreeReaderValue<float>> b_tofe;
   std::unique_ptr<TTreeReaderValue<float>> b_ptofy;
   std::unique_ptr<TTreeReaderValue<float>> b_ptofz;
   //there is no slat variable in embedding output :(
   //std::unique_ptr<TTreeReaderValue<float>> b_slat; 
   std::unique_ptr<TTreeReaderValue<float>> b_pltof;
   std::unique_ptr<TTreeReaderValue<float>> b_ttof;

   //TOFw variables branches
   std::unique_ptr<TTreeReaderValue<float>> b_tofwdz;
   std::unique_ptr<TTreeReaderValue<float>> b_tofwdphi;
   std::unique_ptr<TTreeReaderValue<int>> b_striptofw;
   std::unique_ptr<TTreeReaderValue<float>> b_pltofw;
};

#endif /*EMB_TREE_READER_HPP*/
