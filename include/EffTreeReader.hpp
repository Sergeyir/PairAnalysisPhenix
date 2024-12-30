// $HEADER$
//------------------------------------------------------------------------------------------------
//                            EffTreeReader class declaration
//------------------------------------------------------------------------------------------------
// EffTreeReader - efficiency tree reader
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic class for simple post simulation TTree reading and for similar interaction as PHCentralTrack for variables retrieval for efficiencies and systematic uncertainties evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef EFF_TREE_READER_HPP
#define EFF_TREE_READER_HPP

#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

class EffTreeReader
{
   public :

   EffTreeReader(TTreeReader &reader);
   
   //Event variables
   Int_t nch();
   Int_t bbcz();
   Double_t mom_orig(Int_t i);
   Int_t parent_id(Int_t i);
   Int_t primary_id(Int_t i);
   Int_t particle_id(Int_t i);

   //DC-PC1 variables
   Int_t qual(Int_t i);
   Int_t charge(Int_t i);
   Double_t the0(Int_t i);
   Double_t phi0(Int_t i);
   Double_t phi(Int_t i);
   Double_t alpha(Int_t i);
   Double_t mom(Int_t i);
   Double_t zed(Int_t i);

   //TOFe variables
   Double_t tofdz(Int_t i);
   Double_t tofdphi(Int_t i);
   Double_t ptofy(Int_t i);
   Double_t ptofz(Int_t i);
   Double_t slat(Int_t i);
   Double_t pltof(Int_t i);
   Double_t ttof(Int_t i);
   Double_t etof(Int_t i);

   //TOFw
   Double_t tofwdz(Int_t i);
   Double_t tofwdphi(Int_t i);
   Double_t ptofwy(Int_t i);
   Double_t ptofwz(Int_t i);
   Double_t ttofw(Int_t i);
   Double_t pltofw(Int_t i);
   Double_t striptofw(Int_t i);

   //EMCal
   Double_t emcdz(Int_t i);
   Double_t emcdphi(Int_t i);
   Double_t ecore(Int_t i);
   Int_t sect(Int_t i);
   Double_t pemcy(Int_t i);
   Double_t pemcz(Int_t i);
   Double_t plemc(Int_t i);
   Double_t temc(Int_t i);

   //PC1
   Double_t ppc1x(Int_t i);
   Double_t ppc1y(Int_t i);
   Double_t ppc1z(Int_t i);

   //PC2
   Double_t pc2dz(Int_t i);
   Double_t pc2dphi(Int_t i);

   Double_t ppc2x(Int_t i);
   Double_t ppc2y(Int_t i);
   Double_t ppc2z(Int_t i);

   //PC3
   Double_t pc3dz(Int_t i);
   Double_t pc3dphi(Int_t i);
   Double_t ppc3x(Int_t i);
   Double_t ppc3y(Int_t i);
   Double_t ppc3z(Int_t i);

   //RICH
   Int_t n0(Int_t i);
   Double_t disp(Int_t i);
   Double_t chi2(Int_t i);
   Double_t npe0(Int_t i);

   ~EffTreeReader();
   
   private:

   //Event variables branches
   std::unique_ptr<TTreeReaderValue<int>> b_nch;
   std::unique_ptr<TTreeReaderValue<float>> b_bbcz;
   std::unique_ptr<TTreeReaderArray<float>> b_mom_orig;
   std::unique_ptr<TTreeReaderArray<int>> b_parent_id;
   std::unique_ptr<TTreeReaderArray<int>> b_primary_id;
   std::unique_ptr<TTreeReaderArray<int>> b_particle_id;
   
   //DC-PC1 variables branches
   std::unique_ptr<TTreeReaderArray<int>> b_qual;
   std::unique_ptr<TTreeReaderArray<int>> b_charge;
   std::unique_ptr<TTreeReaderArray<float>> b_the0;
   std::unique_ptr<TTreeReaderArray<float>> b_phi0;
   std::unique_ptr<TTreeReaderArray<float>> b_phi;
   std::unique_ptr<TTreeReaderArray<float>> b_alpha;
   std::unique_ptr<TTreeReaderArray<float>> b_mom;
   std::unique_ptr<TTreeReaderArray<float>> b_zed;

   //TOFe variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_tofdz;
   std::unique_ptr<TTreeReaderArray<float>> b_tofdphi;
   std::unique_ptr<TTreeReaderArray<float>> b_ptofy;
   std::unique_ptr<TTreeReaderArray<float>> b_ptofz;
   std::unique_ptr<TTreeReaderArray<int>> b_slat;
   std::unique_ptr<TTreeReaderArray<float>> b_pltof;
   std::unique_ptr<TTreeReaderArray<float>> b_ttof;
   std::unique_ptr<TTreeReaderArray<float>> b_etof;

   //TOFw variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_tofwdz;
   std::unique_ptr<TTreeReaderArray<float>> b_tofwdphi;
   std::unique_ptr<TTreeReaderArray<float>> b_ptofwy;
   std::unique_ptr<TTreeReaderArray<float>> b_ptofwz;
   std::unique_ptr<TTreeReaderArray<float>> b_ttofw;
   std::unique_ptr<TTreeReaderArray<float>> b_pltofw;
   std::unique_ptr<TTreeReaderArray<int>> b_striptofw;
   
   //EMCal variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_emcdz;
   std::unique_ptr<TTreeReaderArray<float>> b_emcdphi;
   std::unique_ptr<TTreeReaderArray<float>> b_pemcy;
   std::unique_ptr<TTreeReaderArray<float>> b_pemcz;
   std::unique_ptr<TTreeReaderArray<float>> b_plemc;
   std::unique_ptr<TTreeReaderArray<float>> b_temc;
   std::unique_ptr<TTreeReaderArray<float>> b_ecore;
   std::unique_ptr<TTreeReaderArray<int>> b_sect;

   //PC1 variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1x;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1y;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1z;
   
   //PC2 variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_pc2dz;
   std::unique_ptr<TTreeReaderArray<float>> b_pc2dphi;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2x;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2y;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2z;

   //PC3 variables branches
   std::unique_ptr<TTreeReaderArray<float>> b_pc3dz;
   std::unique_ptr<TTreeReaderArray<float>> b_pc3dphi;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3x;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3y;
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3z;
   
   //RICH variables branches
   std::unique_ptr<TTreeReaderArray<int>> b_n0;
   std::unique_ptr<TTreeReaderArray<float>> b_disp;
   std::unique_ptr<TTreeReaderArray<float>> b_chi2;
   std::unique_ptr<TTreeReaderArray<float>> b_npe0;      
};

#endif /*EFF_TREE_READER_HPP*/
