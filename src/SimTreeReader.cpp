/** 
 *  @file   SimTreeReader.cpp
 *  @brief  Contains realisation of class SimTreeReader that can be used to read data from Trees obtained from PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SIM_TREE_READER_CPP
#define SIM_TREE_READER_CPP

#include "../include/SimTreeReader.hpp"

SimTreeReader::SimTreeReader(TTreeReader &reader)
{
   // branch for original momentum
   b_mom_orig = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "mom_orig"));
   // branches for event variables
   b_nch = std::make_unique<TTreeReaderValue<int>>
      (TTreeReaderValue<int>(reader, "nch"));
   b_bbcz = std::make_unique<TTreeReaderValue<float>>
      (TTreeReaderValue<float>(reader, "bbcz"));
   // branches for charged track variables
   b_dcarm = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "dcarm"));
   b_phi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "phi"));
   b_alpha = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "alpha"));
   b_zed = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "zed"));
   b_mom = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "mom"));
   b_the0 = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "the0"));
   b_phi0 = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "phi0"));
   b_nx1hits = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "nx1hits"));
   b_nx2hits = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "nx2hits"));
   b_qual = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "qual"));
   b_charge = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "charge"));
   b_parent_id = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "parent_id"));
   b_primary_id = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "primary_id"));
   b_particle_id = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "particle_id"));
   b_ttof = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ttof"));
   b_ttofw = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ttofw"));
   b_temc = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "temc"));
   b_pltof = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pltof"));
   b_pltofw = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pltofw"));
   b_plemc = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "plemc"));
   b_ptofx = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofx"));
   b_ptofy = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofy"));
   b_ptofz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofz"));
   b_ptofwx = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofwx"));
   b_ptofwy = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofwy"));
   b_ptofwz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptofwz"));
   b_pemcx = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pemcx"));
   b_pemcy = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pemcy"));
   b_pemcz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pemcz"));
   b_ppc1x = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc1x"));
   b_ppc1y = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc1y"));
   b_ppc1z = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc1z"));
   b_ppc2x = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc2x"));
   b_ppc2y = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc2y"));
   b_ppc2z = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc2z"));
   b_ppc3x = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc3x"));
   b_ppc3y = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc3y"));
   b_ppc3z = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ppc3z"));
   b_ptecx = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptecx"));
   b_ptecy = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptecy"));
   b_ptecz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ptecz"));
   b_tofdz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "tofdz"));
   b_tofdphi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "tofdphi"));
   b_tofwdz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "tofwdz"));
   b_tofwdphi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "tofwdphi"));
   b_emcdz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emcdz"));
   b_emcdphi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emcdphi"));
   b_pc2dz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pc2dz"));
   b_pc2dphi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pc2dphi"));
   b_pc3dz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pc3dz"));
   b_pc3dphi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "pc3dphi"));
   b_striptofw = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "striptofw"));
   b_slat = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "slat"));
   b_etof = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "etof"));
   b_ecore = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ecore"));
   b_emce = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emce"));
   b_ecent = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "ecent"));
   b_e9 = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "e9"));
   b_emcchi2 = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emcchi2"));
   b_twrhit = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "twrhit"));
   b_emcdispy = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emcdispy"));
   b_emcdispz = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "emcdispz"));
   b_prob = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "prob"));
   b_sect = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "sect"));
   b_ysect = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "ysect"));
   b_zsect = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "zsect"));
   b_n0 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "n0"));
   b_npe0 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "npe0"));
   b_n1 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "n1"));
   b_npe1 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "npe1"));
   b_n2 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "n2"));
   b_npe2 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "npe2"));
   b_n3 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "n3"));
   b_npe3 = std::make_unique<TTreeReaderArray<short>>
      (TTreeReaderArray<short>(reader, "npe3"));
   b_center_phi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "center_phi"));
   b_center_z = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "center_z"));
   b_cross_phi = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "cross_phi"));
   b_cross_z = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "cross_z"));
   b_disp = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "disp"));
   b_chi2 = std::make_unique<TTreeReaderArray<float>>
      (TTreeReaderArray<float>(reader, "chi2"));
}

// original momentum
float SimTreeReader::mom_orig(const int i) const {return (*b_mom_orig)[i];}
// event variables
int SimTreeReader::nch() const {return *(b_nch->Get());}
float SimTreeReader::bbcz() const {return *(b_bbcz->Get());}
// charged track variables
short SimTreeReader::dcarm(const int i) const {return (*b_dcarm)[i];}
float SimTreeReader::phi(const int i) const {return (*b_phi)[i];}
float SimTreeReader::alpha(const int i) const {return (*b_alpha)[i];}
float SimTreeReader::zed(const int i) const {return (*b_zed)[i];}
float SimTreeReader::mom(const int i) const {return (*b_mom)[i];}
float SimTreeReader::the0(const int i) const {return (*b_the0)[i];}
float SimTreeReader::phi0(const int i) const {return (*b_phi0)[i];}
short SimTreeReader::nx1hits(const int i) const {return (*b_nx1hits)[i];}
short SimTreeReader::nx2hits(const int i) const {return (*b_nx2hits)[i];}
short SimTreeReader::qual(const int i) const {return (*b_qual)[i];}
short SimTreeReader::charge(const int i) const {return (*b_charge)[i];}
short SimTreeReader::parent_id(const int i) const {return (*b_parent_id)[i];}
short SimTreeReader::primary_id(const int i) const {return (*b_primary_id)[i];}
short SimTreeReader::particle_id(const int i) const {return (*b_particle_id)[i];}
float SimTreeReader::ttof(const int i) const {return (*b_ttof)[i];}
float SimTreeReader::ttofw(const int i) const {return (*b_ttofw)[i];}
float SimTreeReader::temc(const int i) const {return (*b_temc)[i];}
float SimTreeReader::pltof(const int i) const {return (*b_pltof)[i];}
float SimTreeReader::pltofw(const int i) const {return (*b_pltofw)[i];}
float SimTreeReader::plemc(const int i) const {return (*b_plemc)[i];}
float SimTreeReader::ptofx(const int i) const {return (*b_ptofx)[i];}
float SimTreeReader::ptofy(const int i) const {return (*b_ptofy)[i];}
float SimTreeReader::ptofz(const int i) const {return (*b_ptofz)[i];}
float SimTreeReader::ptofwx(const int i) const {return (*b_ptofwx)[i];}
float SimTreeReader::ptofwy(const int i) const {return (*b_ptofwy)[i];}
float SimTreeReader::ptofwz(const int i) const {return (*b_ptofwz)[i];}
float SimTreeReader::pemcx(const int i) const {return (*b_pemcx)[i];}
float SimTreeReader::pemcy(const int i) const {return (*b_pemcy)[i];}
float SimTreeReader::pemcz(const int i) const {return (*b_pemcz)[i];}
float SimTreeReader::ppc1x(const int i) const {return (*b_ppc1x)[i];}
float SimTreeReader::ppc1y(const int i) const {return (*b_ppc1y)[i];}
float SimTreeReader::ppc1z(const int i) const {return (*b_ppc1z)[i];}
float SimTreeReader::ppc2x(const int i) const {return (*b_ppc2x)[i];}
float SimTreeReader::ppc2y(const int i) const {return (*b_ppc2y)[i];}
float SimTreeReader::ppc2z(const int i) const {return (*b_ppc2z)[i];}
float SimTreeReader::ppc3x(const int i) const {return (*b_ppc3x)[i];}
float SimTreeReader::ppc3y(const int i) const {return (*b_ppc3y)[i];}
float SimTreeReader::ppc3z(const int i) const {return (*b_ppc3z)[i];}
float SimTreeReader::ptecx(const int i) const {return (*b_ptecx)[i];}
float SimTreeReader::ptecy(const int i) const {return (*b_ptecy)[i];}
float SimTreeReader::ptecz(const int i) const {return (*b_ptecz)[i];}
float SimTreeReader::tofdz(const int i) const {return (*b_tofdz)[i];}
float SimTreeReader::tofdphi(const int i) const {return (*b_tofdphi)[i];}
float SimTreeReader::tofwdz(const int i) const {return (*b_tofwdz)[i];}
float SimTreeReader::tofwdphi(const int i) const {return (*b_tofwdphi)[i];}
float SimTreeReader::emcdz(const int i) const {return (*b_emcdz)[i];}
float SimTreeReader::emcdphi(const int i) const {return (*b_emcdphi)[i];}
float SimTreeReader::pc2dz(const int i) const {return (*b_pc2dz)[i];}
float SimTreeReader::pc2dphi(const int i) const {return (*b_pc2dphi)[i];}
float SimTreeReader::pc3dz(const int i) const {return (*b_pc3dz)[i];}
float SimTreeReader::pc3dphi(const int i) const {return (*b_pc3dphi)[i];}
short SimTreeReader::striptofw(const int i) const {return (*b_striptofw)[i];}
short SimTreeReader::slat(const int i) const {return (*b_slat)[i];}
float SimTreeReader::etof(const int i) const {return (*b_etof)[i];}
float SimTreeReader::ecore(const int i) const {return (*b_ecore)[i];}
float SimTreeReader::emce(const int i) const {return (*b_emce)[i];}
float SimTreeReader::ecent(const int i) const {return (*b_ecent)[i];}
float SimTreeReader::e9(const int i) const {return (*b_e9)[i];}
float SimTreeReader::emcchi2(const int i) const {return (*b_emcchi2)[i];}
short SimTreeReader::twrhit(const int i) const {return (*b_twrhit)[i];}
float SimTreeReader::emcdispy(const int i) const {return (*b_emcdispy)[i];}
float SimTreeReader::emcdispz(const int i) const {return (*b_emcdispz)[i];}
float SimTreeReader::prob(const int i) const {return (*b_prob)[i];}
short SimTreeReader::sect(const int i) const {return (*b_sect)[i];}
short SimTreeReader::ysect(const int i) const {return (*b_ysect)[i];}
short SimTreeReader::zsect(const int i) const {return (*b_zsect)[i];}
short SimTreeReader::n0(const int i) const {return (*b_n0)[i];}
short SimTreeReader::npe0(const int i) const {return (*b_npe0)[i];}
short SimTreeReader::n1(const int i) const {return (*b_n1)[i];}
short SimTreeReader::npe1(const int i) const {return (*b_npe1)[i];}
short SimTreeReader::n2(const int i) const {return (*b_n2)[i];}
short SimTreeReader::npe2(const int i) const {return (*b_npe2)[i];}
short SimTreeReader::n3(const int i) const {return (*b_n3)[i];}
short SimTreeReader::npe3(const int i) const {return (*b_npe3)[i];}
float SimTreeReader::center_phi(const int i) const {return (*b_center_phi)[i];}
float SimTreeReader::center_z(const int i) const {return (*b_center_z)[i];}
float SimTreeReader::cross_phi(const int i) const {return (*b_cross_phi)[i];}
float SimTreeReader::cross_z(const int i) const {return (*b_cross_z)[i];}
float SimTreeReader::disp(const int i) const {return (*b_disp)[i];}
float SimTreeReader::chi2(const int i) const {return (*b_chi2)[i];}

SimTreeReader::~SimTreeReader() {}

#endif /* SIM_TREE_READER_CPP */
