// $HEADER$
//------------------------------------------------------------------------------------------------
//                     EmbTreeReader and EventReader classes realisation
//------------------------------------------------------------------------------------------------
// EmbTreeReader - embedding tree reader
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic classes for simple post simulation TTree reading and for similar interaction as PHCentralTrack for variables retrieval for embedding efficiencies evaluation
 **/
//------------------------------------------------------------------------------------------------

#ifndef EMB_TREE_READER_CPP
#define EMB_TREE_READER_CPP

#include "../include/EmbTreeReader.hpp"

EventReader::EventReader(TTreeReader &reader)
{
   //event variables branches
   b_type = std::make_unique<TTreeReaderValue<short>>(
      TTreeReaderValue<short>(reader, "type"));
   b_geant_id = std::make_unique<TTreeReaderValue<int>>(
      TTreeReaderValue<int>(reader, "partidG"));
   b_cent = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, "bbccent"));
   b_bbcz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, "bbcvtx"));
}

Short_t EventReader::type() {return *(b_type->Get());}
Int_t EventReader::geant_id() {return *(b_geant_id->Get());}
Double_t EventReader::cent() {return *(b_cent->Get());}
Double_t EventReader::bbcz() {return *(b_bbcz->Get());}

EventReader::~EventReader() {};

EmbTreeReader::EmbTreeReader(TTreeReader &reader, const std::string input_type)
{
   //DC-PC1 variables branches
   b_the0 = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("the0" + input_type).c_str()));
   b_zed = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("zed" + input_type).c_str()));
   b_mom = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("mom" + input_type).c_str()));
   b_phi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("phi" + input_type).c_str()));
   b_alpha = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("alpha" + input_type).c_str()));
   b_qual = std::make_unique<TTreeReaderValue<short>>(
      TTreeReaderValue<short>(reader, ("dctrkQual" + input_type).c_str()));
   
   //PC1 variables branches
   b_ppc1x = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc1x" + input_type).c_str()));
   b_ppc1y = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc1y" + input_type).c_str()));
   b_ppc1z = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc1z" + input_type).c_str()));

   //PC2 variables branches
   b_ppc2x = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc2x" + input_type).c_str()));
   b_ppc2y = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc2y" + input_type).c_str()));
   b_ppc2z = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc2z" + input_type).c_str()));
   b_pc2dz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc2dz" + input_type).c_str()));
   b_pc2dphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc2dphi" + input_type).c_str()));
   b_pc2sdz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc2sdz" + input_type).c_str()));
   b_pc2sdphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc2sdphi" + input_type).c_str()));

   //PC3 variables branches
   b_ppc3x = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc3x" + input_type).c_str()));
   b_ppc3y = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc3y" + input_type).c_str()));
   b_ppc3z = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ppc3z" + input_type).c_str()));
   b_pc3dz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc3dz" + input_type).c_str()));
   b_pc3dphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc3dphi" + input_type).c_str()));
   b_pc3sdz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc3sdz" + input_type).c_str()));
   b_pc3sdphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pc3sdphi" + input_type).c_str()));

   //EMCal variables branches
   b_pemcy = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pemcy" + input_type).c_str()));
   b_pemcz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pemcz" + input_type).c_str()));
   b_emcsdz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("emcsdz" + input_type).c_str()));
   b_emcsdphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("emcsdphi" + input_type).c_str()));
   b_ecore = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("emcecore" + input_type).c_str()));

   //TOFe variables branches
   b_tofsdz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("tofsdz" + input_type).c_str()));
   b_tofsdphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("tofsdphi" + input_type).c_str()));
   b_tofe = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("tofe" + input_type).c_str()));
   b_ptofy = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ptofy" + input_type).c_str()));
   b_ptofz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("ptofz" + input_type).c_str()));
   //there is no slat variable in embedding output :(
   /*
   b_slat = std::make_unique<TTreeReaderValue<float>>(
         TTreeReaderValue<float>(reader, ("slat" + input_type).c_str()));
   */
   b_pltof = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pltof" + input_type).c_str()));
   b_ttof = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("toft" + input_type).c_str()));

   //TOFw variables branches
   b_tofwdz = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("tofwdz" + input_type).c_str()));
   b_tofwdphi = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("tofwdphi" + input_type).c_str()));
   b_striptofw = std::make_unique<TTreeReaderValue<int>>(
      TTreeReaderValue<int>(reader, ("tofwstrip" + input_type).c_str()));
   b_pltofw = std::make_unique<TTreeReaderValue<float>>(
      TTreeReaderValue<float>(reader, ("pltofw" + input_type).c_str()));
}
	
//DC-PC1 variables
Double_t EmbTreeReader::zed() {return *(b_zed->Get());}
Double_t EmbTreeReader::the0() {return *(b_the0->Get());}
Double_t EmbTreeReader::mom() {return *(b_mom->Get());}
Double_t EmbTreeReader::phi() {return *(b_phi->Get());}
Double_t EmbTreeReader::alpha() {return *(b_alpha->Get());}
Short_t EmbTreeReader::qual() {return *(b_qual->Get());}

//PC1 variables
Double_t EmbTreeReader::ppc1x() {return *(b_ppc1x->Get());}
Double_t EmbTreeReader::ppc1y() {return *(b_ppc1y->Get());}
Double_t EmbTreeReader::ppc1z() {return *(b_ppc1z->Get());}

//PC2 variables
Double_t EmbTreeReader::ppc2x() {return *(b_ppc2x->Get());}
Double_t EmbTreeReader::ppc2y() {return *(b_ppc2y->Get());}
Double_t EmbTreeReader::ppc2z() {return *(b_ppc2z->Get());}
Double_t EmbTreeReader::pc2dz() {return *(b_pc2dz->Get());}
Double_t EmbTreeReader::pc2dphi() {return *(b_pc2dphi->Get());}
Double_t EmbTreeReader::pc2sdz() {return *(b_pc2sdz->Get());}
Double_t EmbTreeReader::pc2sdphi() {return *(b_pc2sdphi->Get());}

//PC3 variables
Double_t EmbTreeReader::ppc3x() {return *(b_ppc3x->Get());}
Double_t EmbTreeReader::ppc3y() {return *(b_ppc3y->Get());}
Double_t EmbTreeReader::ppc3z() {return *(b_ppc3z->Get());}
Double_t EmbTreeReader::pc3dz() {return *(b_pc3dz->Get());}
Double_t EmbTreeReader::pc3dphi() {return *(b_pc3dphi->Get());}
Double_t EmbTreeReader::pc3sdz() {return *(b_pc3sdz->Get());}
Double_t EmbTreeReader::pc3sdphi() {return *(b_pc3sdphi->Get());}

//EMCal variables
Double_t EmbTreeReader::emcsdz() {return *(b_emcsdz->Get());}
Double_t EmbTreeReader::emcsdphi() {return *(b_emcsdphi->Get());}
Double_t EmbTreeReader::pemcy() {return *(b_pemcy->Get());}
Double_t EmbTreeReader::pemcz() {return *(b_pemcz->Get());}
Double_t EmbTreeReader::ecore() {return *(b_ecore->Get());}

//TOFe variables
Double_t EmbTreeReader::tofsdz() {return *(b_tofsdz->Get());}
Double_t EmbTreeReader::tofsdphi() {return *(b_tofsdphi->Get());}
Double_t EmbTreeReader::tofe() {return *(b_tofe->Get());}
Double_t EmbTreeReader::ptofy() {return *(b_ptofy->Get());}
Double_t EmbTreeReader::ptofz() {return *(b_ptofz->Get());}
//there is no slat variable in embedding output :(
//Double_t EmbTreeReader::slat() {return *(b_slat->Get());}
Double_t EmbTreeReader::pltof() {return *(b_pltof->Get());}
Double_t EmbTreeReader::ttof() {return *(b_ttof->Get());}

//TOFw variables
Double_t EmbTreeReader::tofwdz() {return *(b_tofwdz->Get());}
Double_t EmbTreeReader::tofwdphi() {return *(b_tofwdphi->Get());}
Int_t EmbTreeReader::striptofw() {return *(b_striptofw->Get());}
Double_t EmbTreeReader::pltofw() {return *(b_pltofw->Get());}

EmbTreeReader::~EmbTreeReader() {}

#endif
