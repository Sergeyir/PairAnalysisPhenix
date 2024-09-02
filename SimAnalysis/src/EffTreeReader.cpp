// $SOURCE$
//------------------------------------------------------------------------------------------------
//                            EffTreeReader class realisation
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

#ifndef EFF_TREE_READER_CPP
#define EFF_TREE_READER_CPP

#include "../include/EffTreeReader.hpp"

EffTreeReader::EffTreeReader(TTreeReader &reader)
{
   //Event variables branches
   b_nch = std::make_unique<TTreeReaderValue<int>>(
      TTreeReaderValue<int>(reader, "nch"));
   b_bbcz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, "bbcz"));
   b_mom_orig = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "mom_orig"));
   b_parent_id = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "parent_id"));
   b_primary_id = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "primary_id"));
   b_particle_id = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "particle_id"));
   
   //DC-PC1 variables branches
   b_qual = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "qual"));
   b_charge = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "charge"));
   b_the0 = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "the0"));
   b_phi0 = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "phi0"));
   b_phi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "phi"));
   b_alpha = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "alpha"));
   b_mom = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "mom"));
   b_zed = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "zed"));

   //TOFe variables branches
   b_tofdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofdz"));
   b_tofdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofdphi"));
   b_tofsdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofsdz"));
   b_tofsdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofsdphi"));
   b_ptofy = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ptofy"));
   b_ptofz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ptofz"));

   //TOFw variables branches
   b_tofwdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofwdz"));
   b_tofwdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofwdphi"));
   b_tofwsdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofwsdz"));
   b_tofwsdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "tofwsdphi"));
   b_ttofw = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ttofw"));
   b_pltofw = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pltofw"));
   b_striptofw = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "striptofw"));		
   b_slat = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "slat"));
   b_pltof = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pltof"));
   b_ttof = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ttof"));
   b_etof = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "etof"));
   b_ptofz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ptofz"));
   
   //EMCal variables branches
   b_emcdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "emcdz"));
   b_emcdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "emcdphi"));
   b_emcsdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "emcsdz"));
   b_emcsdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "emcsdphi"));
   b_pemcy = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pemcy"));
   b_pemcz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pemcz"));
   b_plemc = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "plemc"));
   b_temc = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "temc"));
   
   //PC1 variables branches
   b_ppc1x = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc1x"));
   b_ppc1y = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc1y"));
   b_ppc1z = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc1z"));
   
   //PC2 variables branches
   b_pc2dz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc2dz"));
   b_pc2dphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc2dphi"));
   b_pc2sdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc2sdz"));
   b_pc2sdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc2sdphi"));
   b_ppc2x = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc2x"));
   b_ppc2y = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc2y"));
   b_ppc2z = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc2z"));

   //PC3 variables branches
   b_pc3dz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc3dz"));
   b_pc3dphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc3dphi"));
   b_pc3sdz = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc3sdz"));
   b_pc3sdphi = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "pc3sdphi"));
   b_ppc3x = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc3x"));
   b_ppc3y = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc3y"));
   b_ppc3z = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ppc3z"));
   
   //RICH variables branches
   b_n0 = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "n0"));
   b_disp = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "disp"));
   b_chi2 = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "chi2"));
   b_npe0 = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "npe0"));
   b_ecore = std::make_unique<TTreeReaderArray<float>>(
      TTreeReaderArray<float>(reader, "ecore"));
   b_sect = std::make_unique<TTreeReaderArray<int>>(
      TTreeReaderArray<int>(reader, "sect"));
}

//Event variables
Int_t EffTreeReader::nch() {return *(b_nch->Get());}
Int_t EffTreeReader::bbcz() {return *(b_bbcz->Get());}
Double_t EffTreeReader::mom_orig(Int_t i) {return (*b_mom_orig)[i];}
Int_t EffTreeReader::parent_id(Int_t i) {return (*b_parent_id)[i];}
Int_t EffTreeReader::primary_id(Int_t i) {return (*b_primary_id)[i];}
Int_t EffTreeReader::particle_id(Int_t i) {return (*b_particle_id)[i];}

//PC1-DC variables
Int_t EffTreeReader::qual(Int_t i) {return (*b_qual)[i];}
Int_t EffTreeReader::charge(Int_t i) {return (*b_charge)[i];}
Double_t EffTreeReader::the0(Int_t i) {return (*b_the0)[i];}
Double_t EffTreeReader::phi0(Int_t i) {return (*b_phi0)[i];}
Double_t EffTreeReader::phi(Int_t i) {return (*b_phi)[i];}
Double_t EffTreeReader::alpha(Int_t i) {return (*b_alpha)[i];}
Double_t EffTreeReader::mom(Int_t i) {return (*b_mom)[i];}
Double_t EffTreeReader::zed(Int_t i) {return (*b_zed)[i];}

//TOFe variables
Double_t EffTreeReader::tofdz(Int_t i) {return (*b_tofdz)[i];}
Double_t EffTreeReader::tofdphi(Int_t i) {return (*b_tofdphi)[i];}
Double_t EffTreeReader::tofsdz(Int_t i) {return (*b_tofsdz)[i];}
Double_t EffTreeReader::tofsdphi(Int_t i) {return (*b_tofsdphi)[i];}
Double_t EffTreeReader::ptofy(Int_t i) {return (*b_ptofy)[i];}
Double_t EffTreeReader::ptofz(Int_t i) {return (*b_ptofz)[i];}

//TOFw
Double_t EffTreeReader::tofwdz(Int_t i) {return (*b_tofwdz)[i];}
Double_t EffTreeReader::tofwdphi(Int_t i) {return (*b_tofwdphi)[i];}
Double_t EffTreeReader::tofwsdz(Int_t i) {return (*b_tofwsdz)[i];}
Double_t EffTreeReader::tofwsdphi(Int_t i) {return (*b_tofwsdphi)[i];}
Double_t EffTreeReader::ttofw(Int_t i) {return (*b_ttofw)[i];}
Double_t EffTreeReader::pltofw(Int_t i) {return (*b_pltofw)[i];}
Double_t EffTreeReader::striptofw(Int_t i) {return (*b_striptofw)[i];}
Double_t EffTreeReader::slat(Int_t i) {return (*b_slat)[i];}
Double_t EffTreeReader::pltof(Int_t i) {return (*b_pltof)[i];}
Double_t EffTreeReader::ttof(Int_t i) {return (*b_ttof)[i];}
Double_t EffTreeReader::etof(Int_t i) {return (*b_etof)[i];}

//EMCal
Double_t EffTreeReader::emcdz(Int_t i) {return (*b_emcdz)[i];}
Double_t EffTreeReader::emcdphi(Int_t i) {return (*b_emcdphi)[i];}
Double_t EffTreeReader::emcsdz(Int_t i) {return (*b_emcsdz)[i];}
Double_t EffTreeReader::emcsdphi(Int_t i) {return (*b_emcsdphi)[i];}
Double_t EffTreeReader::ecore(Int_t i) {return (*b_ecore)[i];}
Int_t EffTreeReader::sect(Int_t i) {return (*b_sect)[i];}
Double_t EffTreeReader::pemcy(Int_t i) {return (*b_pemcy)[i];}
Double_t EffTreeReader::pemcz(Int_t i) {return (*b_pemcz)[i];}
Double_t EffTreeReader::plemc(Int_t i) {return (*b_plemc)[i];}
Double_t EffTreeReader::temc(Int_t i) {return (*b_temc)[i];}

//PC1

Double_t EffTreeReader::ppc1x(Int_t i) {return (*b_ppc1x)[i];}
Double_t EffTreeReader::ppc1y(Int_t i) {return (*b_ppc1y)[i];}
Double_t EffTreeReader::ppc1z(Int_t i) {return (*b_ppc1z)[i];}

//PC2
Double_t EffTreeReader::pc2dz(Int_t i) {return (*b_pc2dz)[i];}
Double_t EffTreeReader::pc2dphi(Int_t i) {return (*b_pc2dphi)[i];}
Double_t EffTreeReader::pc2sdz(Int_t i) {return (*b_pc2sdz)[i];}
Double_t EffTreeReader::pc2sdphi(Int_t i) {return (*b_pc2sdphi)[i];}

Double_t EffTreeReader::ppc2x(Int_t i) {return (*b_ppc2x)[i];}
Double_t EffTreeReader::ppc2y(Int_t i) {return (*b_ppc2y)[i];}
Double_t EffTreeReader::ppc2z(Int_t i) {return (*b_ppc2z)[i];}

//PC3
Double_t EffTreeReader::pc3dz(Int_t i) {return (*b_pc3dz)[i];}
Double_t EffTreeReader::pc3dphi(Int_t i) {return (*b_pc3dphi)[i];}
Double_t EffTreeReader::pc3sdz(Int_t i) {return (*b_pc3sdz)[i];}
Double_t EffTreeReader::pc3sdphi(Int_t i) {return (*b_pc3sdphi)[i];}
Double_t EffTreeReader::ppc3x(Int_t i) {return (*b_ppc3x)[i];}
Double_t EffTreeReader::ppc3y(Int_t i) {return (*b_ppc3y)[i];}
Double_t EffTreeReader::ppc3z(Int_t i) {return (*b_ppc3z)[i];}

//RICH
Int_t EffTreeReader::n0(Int_t i) {return (*b_n0)[i];}
Double_t EffTreeReader::disp(Int_t i) {return (*b_disp)[i];}
Double_t EffTreeReader::chi2(Int_t i) {return (*b_chi2)[i];}
Double_t EffTreeReader::npe0(Int_t i) {return (*b_npe0)[i];}

EffTreeReader::~EffTreeReader() {}

#endif
