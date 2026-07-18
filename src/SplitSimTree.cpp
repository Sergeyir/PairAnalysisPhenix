/** 
 *  @file   AnalyzeResonance.hpp
 *  @brief  Contains declarations of functions and variables that are used for analysis of a resonance from a trees acquired from the PHENIX simulation
 *
 *  This file is a part of a project PairAnalysisPhenix (https://github.com/Sergeyir/PairAnalysisPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SPLIT_SIM_TREE_CPP
#define SPLIT_SIM_TREE_CPP

#include "SplitSimTree.hpp"

// this namespace is only used so that documentation does not become a mess
// so there is no need to enforce the contents inside of it 
// being accessed only via the scope resolution operator in this file
using namespace SplitSimTree;

int main(int argc, char **argv)
{
   if (argc != 4) 
   {
      std::string errMsg = "Expected 3 parameters while " + std::to_string(argc - 1) + " ";
      errMsg += "parameter(s) were provided \n Usage: bin/SplitSimTree ";
      errMsg += "treeFileName.root pTMin pTMax";
      CppTools::PrintError(errMsg);
   }
 
   CppTools::CheckInputFile(argv[1]);
 
   std::filesystem::create_directories("tmp");

   pTMin = std::stod(argv[2]);
   pTMax = std::stod(argv[3]);

   TFile *inputFile = TFile::Open(argv[1]);

   numberOfEvents = ((static_cast<TTree *>(inputFile->Get("Tree"))->GetEntries()));

   ProgressBar pBar{"BLOCK"};

   TTreeReader reader("Tree", inputFile);
   SimTreeReader simCNT(reader);

   const std::string outputFileName = "tmp/split_tree_" + std::string(argv[2]) + "-" + 
                                      std::string(argv[3]) + ".root";

   TFile outputFile(outputFileName.c_str(), "RECREATE");
   outputFile.SetCompressionLevel(6);

   while (reader.Next())
   {
      pBar.Print(static_cast<double>(numberOfCalls)/
                 static_cast<double>(numberOfEvents));

      numberOfCalls++;
      const double origPT = sqrt(pow(simCNT.mom_orig(0), 2) + pow(simCNT.mom_orig(1), 2));

      if (origPT < pTMin || origPT > pTMax) continue;

      int nch = simCNT.nch();

      Container::nch.SetValue(nch);
      Container::bbcz.SetValue(simCNT.bbcz());

      for (int i = 0; i < 3; i++) {Container::mom_orig.SetValue(simCNT.mom_orig(i), i);}

      distrOrigPT.Fill(origPT);
      
      for (int i = 0; i < nch; i++)
      {
         distrOrigPTVsPT.Fill(origPT, simCNT.mom(0)*sin(simCNT.the0(0)));

         Container::dcarm.SetValue(simCNT.dcarm(i), i);
         Container::phi.SetValue(simCNT.phi(i), i);
         Container::alpha.SetValue(simCNT.alpha(i), i);
         Container::zed.SetValue(simCNT.zed(i), i);
         Container::mom.SetValue(simCNT.mom(i), i);

         Container::phi0.SetValue(simCNT.phi0(i), i);
         Container::the0.SetValue(simCNT.the0(i), i);

         Container::nx1hits.SetValue(simCNT.nx1hits(i), i);
         Container::nx2hits.SetValue(simCNT.nx2hits(i), i);
         
         Container::qual.SetValue(simCNT.qual(i), i);
         Container::charge .SetValue(simCNT.charge(i), i);

         Container::parent_id.SetValue(simCNT.parent_id(i), i);
         Container::particle_id.SetValue(simCNT.particle_id(i), i);
         Container::primary_id.SetValue(simCNT.primary_id(i), i);

         Container::ttof.SetValue(simCNT.ttof(i), i);
         Container::ttofw.SetValue(simCNT.ttofw(i), i);
         Container::temc.SetValue(simCNT.temc(i), i);

         Container::pltof.SetValue(simCNT.pltof(i), i);
         Container::pltofw.SetValue(simCNT.pltofw(i), i);
         Container::plemc.SetValue(simCNT.plemc(i), i);

         Container::ptofx.SetValue(simCNT.ptofx(i), i);
         Container::ptofy.SetValue(simCNT.ptofy(i), i);
         Container::ptofz.SetValue(simCNT.ptofz(i), i);

         Container::ptofwx.SetValue(simCNT.ptofwx(i), i);
         Container::ptofwy.SetValue(simCNT.ptofwy(i), i);
         Container::ptofwz.SetValue(simCNT.ptofwz(i), i);

         Container::pemcx.SetValue(simCNT.pemcx(i), i);
         Container::pemcy.SetValue(simCNT.pemcy(i), i);
         Container::pemcz.SetValue(simCNT.pemcz(i), i);

         Container::ppc1x.SetValue(simCNT.ppc1x(i), i);
         Container::ppc1y.SetValue(simCNT.ppc1y(i), i);
         Container::ppc1z.SetValue(simCNT.ppc1z(i), i);

         Container::ppc2x.SetValue(simCNT.ppc2x(i), i);
         Container::ppc2y.SetValue(simCNT.ppc2y(i), i);
         Container::ppc2z.SetValue(simCNT.ppc2z(i), i);

         Container::ppc3x.SetValue(simCNT.ppc3x(i), i);
         Container::ppc3y.SetValue(simCNT.ppc3y(i), i);
         Container::ppc3z.SetValue(simCNT.ppc3z(i), i);

         Container::ptecx.SetValue(simCNT.ptecx(i), i);
         Container::ptecy.SetValue(simCNT.ptecy(i), i);
         Container::ptecz.SetValue(simCNT.ptecz(i), i);

         Container::tofdz.SetValue(simCNT.tofdz(i), i);
         Container::tofdphi.SetValue(simCNT.tofdphi(i), i);

         Container::tofwdz.SetValue(simCNT.tofwdz(i), i);
         Container::tofwdphi.SetValue(simCNT.tofwdphi(i), i);

         Container::emcdz.SetValue(simCNT.emcdz(i), i);
         Container::emcdphi.SetValue(simCNT.emcdphi(i), i);
         
         Container::pc2dphi.SetValue(simCNT.pc2dphi(i), i);
         Container::pc2dz.SetValue(simCNT.pc2dz(i), i);
         
         Container::pc3dphi.SetValue(simCNT.pc3dphi(i), i);
         Container::pc3dz.SetValue(simCNT.pc3dz(i), i);

         Container::striptofw.SetValue(simCNT.striptofw(i), i);
         Container::slat.SetValue(simCNT.slat(i), i);

         Container::etof.SetValue(simCNT.etof(i), i);
         Container::ecore.SetValue(simCNT.ecore(i), i);
         Container::emce.SetValue(simCNT.emce(i), i);
         Container::ecent.SetValue(simCNT.ecent(i), i);
         Container::e9.SetValue(simCNT.e9(i), i);
         Container::emcchi2.SetValue(simCNT.emcchi2(i), i);
         Container::twrhit.SetValue(simCNT.twrhit(i), i);
         Container::emcdispy.SetValue(simCNT.emcdispy(i), i);
         Container::emcdispz.SetValue(simCNT.emcdispz(i), i);
         Container::prob.SetValue(simCNT.prob(i), i);

         Container::sect.SetValue(simCNT.sect(i), i);
         Container::ysect.SetValue(simCNT.ysect(i), i);
         Container::zsect.SetValue(simCNT.zsect(i), i);
         
         Container::n0.SetValue(simCNT.n0(i), i);
         Container::npe0.SetValue(simCNT.npe0(i), i);
         Container::n1.SetValue(simCNT.n1(i), i);
         Container::npe1.SetValue(simCNT.npe1(i), i);
         Container::n2.SetValue(simCNT.n2(i), i);
         Container::npe2.SetValue(simCNT.npe2(i), i);
         Container::n3.SetValue(simCNT.n3(i), i);
         Container::npe3.SetValue(simCNT.npe3(i), i);
         Container::center_phi.SetValue(simCNT.center_phi(i), i);
         Container::center_z.SetValue(simCNT.center_z(i), i);
         Container::cross_phi.SetValue(simCNT.cross_phi(i), i);
         Container::cross_z.SetValue(simCNT.cross_z(i), i);
         Container::disp.SetValue(simCNT.disp(i), i);
         Container::chi2.SetValue(simCNT.chi2(i), i);
      }

      tree.Fill();

   }

   HistHolder::Write();
   tree.Write();
   distrOrigPT.Write();
   distrOrigPTVsPT.Write();

   outputFile.Close();
 
   pBar.Finish();
   CppTools::PrintInfo("SplitSimTree has finished running; "\
                       "split tree was written as " + outputFileName);

   return 0;
}

void SplitSimTree::HistHolder::Write()
{
   for (unsigned int i = 0; i < vec.size(); i++)
   {
      vec[i]->Write();
   }
}

void SplitSimTree::HistHolder::Clear()
{
   for (unsigned int i = 0; i < vec.size(); i++) {delete vec[i];}
   vec.clear();
}

unsigned int SplitSimTree::HistHolder::CurrentPos() 
{
   return vec.size();
}

SplitSimTree::SimVarInt::SimVarInt(const std::string& name, const double xMin, 
                                   const double xMax, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "/I").c_str());
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos() - 1;
}

void SplitSimTree::SimVarInt::SetValue(const int newVal)
{
   val = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

SplitSimTree::SimVarShort::SimVarShort(const std::string& name, const double xMin, 
                                       const double xMax, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "/S").c_str());
   
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos() - 1;
}

void SplitSimTree::SimVarShort::SetValue(const short newVal)
{
   val = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

SplitSimTree::SimVarF::SimVarF(const std::string& name, const double xMin, 
                               const double xMax, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "/F").c_str());
   
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos()-1;
}

void SplitSimTree::SimVarF::SetValue(const float newVal)
{
   val = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

SplitSimTree::SimVarIntArr::SimVarIntArr(const std::string& name, 
                                         const double xMin, const double xMax, 
                                         const std::string& size, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "[" + size + "]/I").c_str());
   
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos()-1;
}
void SplitSimTree::SimVarIntArr::SetValue(const int newVal, const unsigned int index)
{
   val[index] = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

SplitSimTree::SimVarShortArr::SimVarShortArr(const std::string& name, 
                                             const double xMin, const double xMax, 
                                             const std::string& size, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "[" + size + "]/S").c_str());
   
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos()-1;
}

void SplitSimTree::SimVarShortArr::SetValue(const short newVal, const unsigned int index)
{
   val[index] = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

SplitSimTree::SimVarFArr::SimVarFArr(const std::string& name, 
                                     const double xMin, const double xMax, 
                                     const std::string& size, const int nbins)
{
   tree.Branch(name.c_str(), &val, (name + "[" + size + "]/F").c_str());
   
   HistHolder::vec.push_back(new TH1D(name.c_str(), name.c_str(), nbins, xMin, xMax));
   histHolderIndex = HistHolder::CurrentPos()-1;
}

void SplitSimTree::SimVarFArr::SetValue(const float newVal, const unsigned int index)
{
   val[index] = newVal;
   HistHolder::vec[histHolderIndex]->Fill(static_cast<double>(newVal));
}

#endif /* SPLIT_SIM_TREE_CPP */
