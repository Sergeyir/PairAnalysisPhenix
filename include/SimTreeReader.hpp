/** 
 *  @file   SimTreeReader.hpp
 *  @brief  Contains declaration of class SimTreeReader that can be used to read data from Trees obtained from PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/

#ifndef SIM_TREE_READER_HPP
#define SIM_TREE_READER_HPP

#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

/*! @class SimTreeReader
 * @brief Class SimTreeReadercan be used to read data from Trees that were obtained from PHENIX simulation 
 *
 * See AnalyzeSingleTrack.cpp on its usage
 */
class SimTreeReader
{
   public :
   /*! @brief Constructor with parameters
    * 
    * @param[in] reader TTreeReader variable that will be used to point branches of a TTree to
    */
   SimTreeReader(TTreeReader &reader);
   /// original momentum (i=0,1,2 for px,py,pz respectively)
   float mom_orig(int i);
   /// number of charged particles
   int nch(); 
   /// bbcz (same as in PHCentralTrack)
   float bbcz();
   /// dcarm (same as in PHCentralTrack)
   short dcarm(int i);
   /// phi (same as in PHCentralTrack)
   float phi(int i);
   /// alpha (same as in PHCentralTrack)
   float alpha(int i);
   /// zed (same as in PHCentralTrack)
   float zed(int i);
   /// mom (same as in PHCentralTrack)
   float mom(int i);
   /// the0 (same as in PHCentralTrack)
   float the0(int i);
   /// phi0 (same as in PHCentralTrack)
   float phi0(int i);
   /// nx1hits (same as in PHCentralTrack)
   short nx1hits(int i);
   /// nx2hits (same as in PHCentralTrack)
   short nx2hits(int i);
   /// qual (same as in PHCentralTrack)
   short qual(int i);
   /// charge (same as in PHCentralTrack)
   short charge(int i);
   /// parent_id - id of a parent particle (GEANT id scheme)
   short parent_id(int i);
   /// primary_id - id of an original particle (GEANT id scheme)
   short primary_id(int i);
   /// particle_id - id of a track (GEANT id scheme)
   short particle_id(int i);
   /// ttof (same as in PHCentralTrack)
   float ttof(int i);
   /// ttofw (same as in PHCentralTrack)
   float ttofw(int i);
   /// temc (same as in PHCentralTrack)
   float temc(int i);
   /// pltof (same as in PHCentralTrack)
   float pltof(int i);
   /// pltofw (same as in PHCentralTrack)
   float pltofw(int i);
   /// plemc (same as in PHCentralTrack)
   float plemc(int i);
   /// ptofx (same as in PHCentralTrack)
   float ptofx(int i);
   /// ptofy (same as in PHCentralTrack)
   float ptofy(int i);
   /// ptofz (same as in PHCentralTrack)
   float ptofz(int i);
   /// ptofwx (same as in PHCentralTrack)
   float ptofwx(int i);
   /// ptofwy (same as in PHCentralTrack)
   float ptofwy(int i);
   /// ptofwz (same as in PHCentralTrack)
   float ptofwz(int i);
   /// pemcx (same as in PHCentralTrack)
   float pemcx(int i);
   /// pemcy (same as in PHCentralTrack)
   float pemcy(int i);
   /// pemcz (same as in PHCentralTrack)
   float pemcz(int i);
   /// ppc1x (same as in PHCentralTrack)
   float ppc1x(int i);
   /// ppc1y (same as in PHCentralTrack)
   float ppc1y(int i);
   /// ppc1z (same as in PHCentralTrack)
   float ppc1z(int i);
   /// ppc2x (same as in PHCentralTrack)
   float ppc2x(int i);
   /// ppc2y (same as in PHCentralTrack)
   float ppc2y(int i);
   /// ppc2z (same as in PHCentralTrack)
   float ppc2z(int i);
   /// ppc3x (same as in PHCentralTrack)
   float ppc3x(int i);
   /// ppc3y (same as in PHCentralTrack)
   float ppc3y(int i);
   /// ppc3z (same as in PHCentralTrack)
   float ppc3z(int i);
   /// ptecx (same as in PHCentralTrack)
   float ptecx(int i);
   /// ptecy (same as in PHCentralTrack)
   float ptecy(int i);
   /// ptecz (same as in PHCentralTrack)
   float ptecz(int i);
   /// tofdz (same as in PHCentralTrack)
   float tofdz(int i);
   /// tofdphi (same as in PHCentralTrack)
   float tofdphi(int i);
   /// tofwdz (same as in PHCentralTrack)
   float tofwdz(int i);
   /// tofwdphi (same as in PHCentralTrack)
   float tofwdphi(int i);
   /// emcdz (same as in PHCentralTrack)
   float emcdz(int i);
   /// emcdphi (same as in PHCentralTrack)
   float emcdphi(int i);
   /// pc2dz (same as in PHCentralTrack)
   float pc2dz(int i);
   /// pc2dphi (same as in PHCentralTrack)
   float pc2dphi(int i);
   /// pc3dz (same as in PHCentralTrack)
   float pc3dz(int i);
   /// pc3dphi (same as in PHCentralTrack)
   float pc3dphi(int i);
   /// striptofw (same as in PHCentralTrack)
   short striptofw(int i);
   /// slat (same as in PHCentralTrack)
   short slat(int i);
   /// etof (same as in PHCentralTrack)
   float etof(int i);
   /// ecore (same as in PHCentralTrack)
   float ecore(int i);
   /// emce (same as in PHCentralTrack)
   float emce(int i);
   /// ecent (same as in PHCentralTrack)
   float ecent(int i);
   /// e9 (same as in PHCentralTrack)
   float e9(int i);
   /// emcchi2 (same as in PHCentralTrack)
   float emcchi2(int i);
   /// twrhit (same as in PHCentralTrack)
   short twrhit(int i);
   /// emcdispy (same as in PHCentralTrack)
   float emcdispy(int i);
   /// emcdispz (same as in PHCentralTrack)
   float emcdispz(int i);
   /// prob (same as in PHCentralTrack)
   float prob(int i);
   /// sect (same as in PHCentralTrack)
   short sect(int i);
   /// ysect (same as in PHCentralTrack)
   short ysect(int i);
   /// zsect (same as in PHCentralTrack)
   short zsect(int i);
   /// n0 (same as in PHCentralTrack)
   short n0(int i);
   /// npe0 (same as in PHCentralTrack)
   short npe0(int i);
   /// n1 (same as in PHCentralTrack)
   short n1(int i);
   /// npe1 (same as in PHCentralTrack)
   short npe1(int i);
   /// n2 (same as in PHCentralTrack)
   short n2(int i);
   /// npe2 (same as in PHCentralTrack)
   short npe2(int i);
   /// n3 (same as in PHCentralTrack)
   short n3(int i);
   /// npe3 (same as in PHCentralTrack)
   short npe3(int i);
   /// center_phi (same as in PHCentralTrack)
   float center_phi(int i);
   /// center_z (same as in PHCentralTrack)
   float center_z(int i);
   /// cross_phi (same as in PHCentralTrack)
   float cross_phi(int i);
   /// cross_z (same as in PHCentralTrack)
   float cross_z(int i);
   /// disp (same as in PHCentralTrack)
   float disp(int i);
   /// chi2 (same as in PHCentralTrack)
   float chi2(int i);
   /// Default desctructor
   ~SimTreeReader();

   private:
   /// branch for original momentum variables
   std::unique_ptr<TTreeReaderArray<float>> b_mom_orig;
   /// branch for number of charged particles variable
   std::unique_ptr<TTreeReaderValue<int>> b_nch;
   /// branch for bbcz variable
   std::unique_ptr<TTreeReaderValue<float>> b_bbcz;
   /// branch for dcarm variables
   std::unique_ptr<TTreeReaderArray<short>> b_dcarm;
   /// branch for phi variables
   std::unique_ptr<TTreeReaderArray<float>> b_phi;
   /// branch for alpha variables
   std::unique_ptr<TTreeReaderArray<float>> b_alpha;
   /// branch for zed variables
   std::unique_ptr<TTreeReaderArray<float>> b_zed;
   /// branch for mom variables
   std::unique_ptr<TTreeReaderArray<float>> b_mom;
   /// branch for the0 variables
   std::unique_ptr<TTreeReaderArray<float>> b_the0;
   /// branch for the0 variables
   std::unique_ptr<TTreeReaderArray<float>> b_phi0;
   /// branch for phi0 variables
   std::unique_ptr<TTreeReaderArray<short>> b_nx1hits;
   /// branch for nx1hits variables
   std::unique_ptr<TTreeReaderArray<short>> b_nx2hits;
   /// branch for qual variables
   std::unique_ptr<TTreeReaderArray<short>> b_qual;
   /// branch for charge variables
   std::unique_ptr<TTreeReaderArray<short>> b_charge;
   /// branch for parent_id variables
   std::unique_ptr<TTreeReaderArray<short>> b_parent_id;
   /// branch for primary_id variables
   std::unique_ptr<TTreeReaderArray<short>> b_primary_id;
   /// branch for particle_id variables
   std::unique_ptr<TTreeReaderArray<short>> b_particle_id;
   /// branch for ttof variables
   std::unique_ptr<TTreeReaderArray<float>> b_ttof;
   /// branch for ttofw variables
   std::unique_ptr<TTreeReaderArray<float>> b_ttofw;
   /// branch for temc variables
   std::unique_ptr<TTreeReaderArray<float>> b_temc;
   /// branch for pltof variables
   std::unique_ptr<TTreeReaderArray<float>> b_pltof;
   /// branch for pltofw variables
   std::unique_ptr<TTreeReaderArray<float>> b_pltofw;
   /// branch for plemc variables
   std::unique_ptr<TTreeReaderArray<float>> b_plemc;
   /// branch for ptofx variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofx;
   /// branch for ptofy variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofy;
   /// branch for ptofz variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofz;
   /// branch for ptofwx variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofwx;
   /// branch for ptofwy variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofwy;
   /// branch for ptofwz variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptofwz;
   /// branch for pemcx variables
   std::unique_ptr<TTreeReaderArray<float>> b_pemcx;
   /// branch for pemcy variables
   std::unique_ptr<TTreeReaderArray<float>> b_pemcy;
   /// branch for pemcz variables
   std::unique_ptr<TTreeReaderArray<float>> b_pemcz;
   /// branch for ppc1x variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1x;
   /// branch for ppc1y variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1y;
   /// branch for ppc1z variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc1z;
   /// branch for ppc2x variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2x;
   /// branch for ppc2y variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2y;
   /// branch for ppc2z variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc2z;
   /// branch for ppc3x variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3x;
   /// branch for ppc3y variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3y;
   /// branch for ppc3z variables
   std::unique_ptr<TTreeReaderArray<float>> b_ppc3z;
   /// branch for tecx variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptecx;
   /// branch for tecy variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptecy;
   /// branch for tecz variables
   std::unique_ptr<TTreeReaderArray<float>> b_ptecz;
   /// branch for tofdz variables
   std::unique_ptr<TTreeReaderArray<float>> b_tofdz;
   /// branch for tofdphi variables
   std::unique_ptr<TTreeReaderArray<float>> b_tofdphi;
   /// branch for tofwdz variables
   std::unique_ptr<TTreeReaderArray<float>> b_tofwdz;
   /// branch for tofwdphi variables
   std::unique_ptr<TTreeReaderArray<float>> b_tofwdphi;
   /// branch for emcdz variables
   std::unique_ptr<TTreeReaderArray<float>> b_emcdz;
   /// branch for emcdphi variables
   std::unique_ptr<TTreeReaderArray<float>> b_emcdphi;
   /// branch for pc2dz variables
   std::unique_ptr<TTreeReaderArray<float>> b_pc2dz;
   /// branch for pc3dphi variables
   std::unique_ptr<TTreeReaderArray<float>> b_pc2dphi;
   /// branch for pc3dz variables
   std::unique_ptr<TTreeReaderArray<float>> b_pc3dz;
   /// branch for pc3dphi variables
   std::unique_ptr<TTreeReaderArray<float>> b_pc3dphi;
   /// branch for striptofw variables
   std::unique_ptr<TTreeReaderArray<short>> b_striptofw;
   /// branch for slat variables
   std::unique_ptr<TTreeReaderArray<short>> b_slat;
   /// branch for etof variables
   std::unique_ptr<TTreeReaderArray<float>> b_etof;
   /// branch for ecore variables
   std::unique_ptr<TTreeReaderArray<float>> b_ecore;
   /// branch for emce variables
   std::unique_ptr<TTreeReaderArray<float>> b_emce;
   /// branch for ecent variables
   std::unique_ptr<TTreeReaderArray<float>> b_ecent;
   /// branch for e9 variables
   std::unique_ptr<TTreeReaderArray<float>> b_e9;
   /// branch for emcchi2 variables
   std::unique_ptr<TTreeReaderArray<float>> b_emcchi2;
   /// branch for twrhit variables
   std::unique_ptr<TTreeReaderArray<short>> b_twrhit;
   /// branch for emcdispy variables
   std::unique_ptr<TTreeReaderArray<float>> b_emcdispy;
   /// branch for emcdispz variables
   std::unique_ptr<TTreeReaderArray<float>> b_emcdispz;
   /// branch for prob variables
   std::unique_ptr<TTreeReaderArray<float>> b_prob;
   /// branch for sect variables
   std::unique_ptr<TTreeReaderArray<short>> b_sect;
   /// branch for ysect variables
   std::unique_ptr<TTreeReaderArray<short>> b_ysect;
   /// branch for zsect variables
   std::unique_ptr<TTreeReaderArray<short>> b_zsect;
   /// branch for n0 variables
   std::unique_ptr<TTreeReaderArray<short>> b_n0;
   /// branch for npe0 variables
   std::unique_ptr<TTreeReaderArray<short>> b_npe0;
   /// branch for n1 variables
   std::unique_ptr<TTreeReaderArray<short>> b_n1;
   /// branch for npe1 variables
   std::unique_ptr<TTreeReaderArray<short>> b_npe1;
   /// branch for n2 variables
   std::unique_ptr<TTreeReaderArray<short>> b_n2;
   /// branch for npe2 variables
   std::unique_ptr<TTreeReaderArray<short>> b_npe2;
   /// branch for n3 variables
   std::unique_ptr<TTreeReaderArray<short>> b_n3;
   /// branch for npe3 variables
   std::unique_ptr<TTreeReaderArray<short>> b_npe3;
   /// branch for center_phi variables
   std::unique_ptr<TTreeReaderArray<float>> b_center_phi;
   /// branch for center_z variables
   std::unique_ptr<TTreeReaderArray<float>> b_center_z;
   /// branch for cross_phi variables
   std::unique_ptr<TTreeReaderArray<float>> b_cross_phi;
   /// branch for cross_z variables
   std::unique_ptr<TTreeReaderArray<float>> b_cross_z;
   /// branch for disp variables
   std::unique_ptr<TTreeReaderArray<float>> b_disp;
   /// branch for chi2 variables
   std::unique_ptr<TTreeReaderArray<float>> b_chi2;
};

#endif /* SIM_TREE_READER_HPP */
