#include "../lib/ErrorHandler.h"
#include "../lib/CentralityTypes.h"
#include "../lib/OutputTool.h"
#include "../lib/Tools.h"
#include "../lib/TCanvasPrinter.h"

struct FitResults
{
   double pion;
   double kaon;
   double proton;
};

struct
{
   std::map<int, TF1*> mean;
   std::map<int, TF1*> sigma;
} MS;

struct
{
   std::string run_name = "Run7AuAu200";
   std::string col_system = "AuAu200";

   std::vector<std::string> det_name_queue =  {"EMCale", "EMCalw"};
   std::vector<std::string> det_hname_queue = {"m2_emcale", "m2_emcalw"};
   std::vector<int> det_min_sect_queue = {2, 0};

   std::vector<std::string> magf_queue{"+-", "-+"};
   std::vector<int> part_id_queue{211, -211, 321, -321};
   std::vector<double> tofw_sys_queue{0.0595147, 0.0595147, 0.0894427, 0.0894427};
   std::vector<double> bg_fit_sys_queue{0., 0., 0.05, 0.05};
   std::vector<double> fit_range_sys_queue{0.02, 0.02, 0.05, 0.05};

   //corrections that account for binning (since yields are calculated by integrating bins)
   //and approximations inaccuracies
   //mostly used for first and last bins of kaons, but pions sometimes also need to be corrected
   double emb;
   std::array<double, 7> addit_mult; //none
      
   AuAu200CTypeMB4 CType;

   const double min_pidw = 0.0;

   const double minimum_pt = 0.5;
   const double maximum_pt = 1.2;

   const double srange = 1.;
   const double veto = 3.;

   bool print_info = 0;

   std::string magf;
   int part_id;
   double other_sys;
   double other_sys_no_TOFw;
   
   std::unique_ptr<TFile> input_file;
   
   std::vector<TGraphErrors *> m2VGr, RatioVGr, RatioVGrSys;
   std::vector<TF1 *> VFunc;
   TLegend legend = TLegend(0.7, 0.12, 0.9, 0.45);
   
   double average_eff;
   double average_eff_sys;

   std::vector<double> ptmin = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};
   std::vector<double> ptmax = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};

   std::string mean_fit_name = "pol1(0)";
   std::string sigma_fit_name = "sqrt(pow([2]/104., 2)*(4.*pow(pol1(0)*x, 2)) + \
      pow([3]/104., 2)*(4.*pow(pol1(0), 2)*(1+pol1(0)/pow(x, 2))) + \
      pow([4]*2.9972e-4/4.95, 2)*(4.*pow(x, 2)*(pol1(0) + pow(x, 2))))";

   double xmin = 1e31;
   double xmax = -1e31;

   double ymin = 1e31;
   double ymax = -1e31;
   
   std::map<int, std::string> id_map =
   {
      {211, "pion"},
      {321, "kaon"},
      {2212, "proton"},
      {-211, "apion"},
      {-321, "akaon"},
      {-2212, "aproton"}
   };

   std::map<int, std::string> id_map_tex =
   {
      {211, "#pi^{+}"},
      {321, "K^{+}"},
      {2212, "p^{+}"},
      {-211, "#pi^{-}"},
      {-321, "K^{-}"},
      {-2212, "p^{-}"}
   };
} Par;

void SetPar(std::string det_name)
{
   Par.emb = 1.;
   Par.addit_mult = std::array<double, 7>{1., 1., 1., 1., 1., 1., 1.};
   
   if (abs(Par.part_id) == 211)
   {
      if (det_name == "EMCale2") Par.emb = 0.932021;
      else if (det_name == "EMCale3") Par.emb = 0.9322112;
      else if (det_name == "EMCalw0") Par.emb = 0.909102;
      else if (det_name == "EMCalw1") Par.emb = 0.926063;
      else if (det_name == "EMCalw2") Par.emb = 0.941958;
      else if (det_name == "EMCalw3") Par.emb = 0.938259;
      else 
      {
         PrintWarning("Detector " + det_name + " doesn't support embedding");
      }
   }
   else if (abs(Par.part_id) == 321)
   {
      if (det_name == "EMCale2") Par.emb = 0.927595;
      else if (det_name == "EMCale3") Par.emb = 0.93126;
      else if (det_name == "EMCalw0") Par.emb = 0.907136;
      else if (det_name == "EMCalw1") Par.emb = 0.927571;
      else if (det_name == "EMCalw2") Par.emb = 0.935768;
      else if (det_name == "EMCalw3") Par.emb = 0.933901;
      else 
      {
         PrintWarning("Detector " + det_name + " doesn't support embedding");
      }
   }
   else PrintError("Unknown Par.part_id: " + to_string(Par.part_id));

   Par.addit_mult = std::array<double, 7>{1., 1., 1., 1., 1., 1., 1.};
   
   if (Par.magf == "+-")
   {
      if (det_name == "EMCale2")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.05, 1., 1., 1., 1., 1., 0.95}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1., 1.1, 1., 1., 1., 1.1, 1.4}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.8, 0.9, 0.85, 1., 0.9, 1., 1.4}; break;
         }
      }
      else if (det_name == "EMCale3")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.92, 1., 1., 1., 1., 1., 1.}; break;
            case -211: Par.addit_mult = std::array<double, 7>{1.05, 1.05, 1., 1., 1., 1., 1.}; break;
            case 321: Par.addit_mult = std::array<double, 7>{0.9, 1., 1.1, 1., 1., 1.2, 1.7}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.85, 1., 0.95, 1., 1., 1., 1.9}; break;
         }
      }
      else if (det_name == "EMCalw0")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.75, 0.75, 0.75, 0.7, 0.7, 0.7, 0.65}; break;
            case -211: Par.addit_mult = std::array<double, 7>{0.95, 0.95, 0.85, 0.85, 0.8, 0.8, 0.7}; break;
            case 321: Par.addit_mult = std::array<double, 7>{0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 1.}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.95}; break;
         }
      }
      else if (det_name == "EMCalw1")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}; break;
            case -211: Par.addit_mult = std::array<double, 7>{1.1, 1.1, 1., 1., 1., 0.9, 0.9}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1.2, 1.1, 1., 1., 1., 1., 1.1}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.9, 1., 1., 0.95, 1., 1., 1.2}; break;
         }
      }
      else if (det_name == "EMCalw2")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.1, 1., 1., 1., 1., 1., 1.}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1.2, 1.2, 0.9, 1., 1., 1., 1.1}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.95, 1., 0.95, 1., 1., 1., 1.2}; break;
         }
      }
      else if (det_name == "EMCalw3")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.05, 1., 1., 1., 1., 1., 0.95}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1., 1.1, 1., 1., 1., 1., 1.2}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.85, 0.9, 0.95, 1., 1., 1.1, 1.2}; break;
         }
      }
   }
   else if (Par.magf == "-+")
   {
      if (det_name == "EMCale2")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.05, 1., 1., 1., 1., 1., 0.95}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1., 1., 1., 1., 0.9, 1., 1.4}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.8, 0.9, 1., 1., 0.9, 0.9, 1.6}; break;
         }
      }
      else if (det_name == "EMCale3")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.92, 1., 1., 1., 1., 1., 1.}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1., 1., 0.95, 1., 1., 1.05, 1.5};  break;
            case -321: Par.addit_mult = std::array<double, 7>{0.85, 0.95, 0.95, 1., 0.9, 1., 1.9}; break;
         }
      }
      else if (det_name == "EMCalw0")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.75, 0.75, 0.75, 0.7, 0.7, 0.7, 0.65}; break;
            case -211: Par.addit_mult = std::array<double, 7>{0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.7}; break;
            case 321: Par.addit_mult = std::array<double, 7>{0.9, 0.85, 0.8, 0.8, 0.8, 0.8, 1.}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.65, 0.7, 0.7, 0.7, 0.7, 0.7, 0.9}; break;
         }
      }
      else if (det_name == "EMCalw1")
      {
         switch (Par.part_id)
         {
            case 211: Par.addit_mult = std::array<double, 7>{0.95, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}; break;
            case -211: Par.addit_mult = std::array<double, 7>{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.85}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1.2, 1.2, 1., 1., 1., 1., 1.2}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 1.2};; break;
         }
      }
      else if (det_name == "EMCalw2")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.1, 1.2, 1., 1., 1., 1., 0.95}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1.2, 1.1, 0.9, 1., 1., 1., 1.15}; break;
            case -321: Par.addit_mult = std::array<double, 7>{1., 1., 0.95, 0.9, 0.9, 1., 1.1}; break;
         }
      }
      else if (det_name == "EMCalw3")
      {
         switch (Par.part_id)
         {
            case -211: Par.addit_mult = std::array<double, 7>{1.05, 1., 1., 1., 1., 1., 0.95}; break;
            case 321: Par.addit_mult = std::array<double, 7>{1., 1.1, 1., 1., 1., 1., 1.2}; break;
            case -321: Par.addit_mult = std::array<double, 7>{0.8, 0.9, 0.8, 0.95, 0.9, 1., 1.1}; break;
         }
      }
   }
   else PrintWarning("No magnetic field " + Par.magf);
}

double GetEffErr(TGraphErrors *gr, const double average)
{
   //calculating standard error of a weighted sample
   double se = 0;
   double weight_norm = 0;
   for (int i = 0; i < gr->GetN(); i++)
   {
      const double val = gr->GetPointY(i);
      const double weight = gr->GetErrorY(i);
      se += pow(average - val, 2)*weight;
      weight_norm += weight;
   }
   
   se = sqrt(se)/weight_norm; 
   return se;
}

double GetPIDW(const double pt)
{
   int charge;
   if (Par.part_id > 0) charge = 1;
   else charge = -1;
   
   FitResults mean, sigma;
   mean.pion = MS.mean[211*charge]->Eval(pt*charge);
   mean.kaon = MS.mean[321*charge]->Eval(pt*charge);
   mean.proton = MS.mean[2212*charge]->Eval(pt*charge);
   
   sigma.pion = MS.sigma[211*charge]->Eval(pt*charge);
   sigma.kaon = MS.sigma[321*charge]->Eval(pt*charge);
   sigma.proton = MS.sigma[2212*charge]->Eval(pt*charge);

   double weight = 0.;

   if (abs(Par.part_id) == 211)
   {
      double low_range = mean.pion - Par.srange*sigma.pion;
      double upp_range = mean.pion + Par.srange*sigma.pion;
      
      double upp_veto = Minimum(mean.kaon - Par.veto*sigma.kaon, 
         mean.proton - Par.veto*sigma.proton);
      if (low_range >= upp_veto) return 0.;

      weight = erf((mean.pion - low_range)/sigma.pion/TMath::Sqrt2())/2. +
         erf((Minimum(upp_range, upp_veto) - mean.pion)/sigma.pion/TMath::Sqrt2())/2.;
   }
   else if (abs(Par.part_id) == 321)
   {
      double low_range = mean.kaon - Par.srange*sigma.kaon;
      double upp_range = mean.kaon + Par.srange*sigma.kaon;
      
      double low_veto = mean.pion + Par.veto*sigma.pion; 
      double upp_veto = mean.proton - Par.veto*sigma.proton; 
      
      if (low_range >= upp_veto || upp_range <= low_veto || upp_veto <= low_veto) return 0.;
      
      weight = erf((mean.kaon - Maximum(low_range, low_veto))/sigma.kaon/TMath::Sqrt2())/2. +
         erf((Minimum(upp_range, upp_veto) - mean.kaon)/sigma.kaon/TMath::Sqrt2())/2.;
   }   
   else if (abs(Par.part_id) == 2212)
   {
      double low_range = mean.proton - Par.srange*sigma.proton;
      double upp_range = mean.proton + Par.srange*sigma.proton;
      
      double low_veto = Maximum(mean.pion + Par.veto*sigma.pion,
         mean.kaon + Par.veto*sigma.kaon);
      
      if (upp_range <= low_veto) return 0.;

      weight = erf((mean.proton - Maximum(low_range, low_veto))/sigma.proton/TMath::Sqrt2())/2. +
         erf((upp_range - mean.proton)/sigma.proton/TMath::Sqrt2())/2.;
   }   
   return weight;
}

double pidw_func(double *x, double *par) {return GetPIDW(x[0]);}

void SetMSFunc(const int id, std::string det_name)
{
   TF1 *mean_fit = new TF1((Par.id_map[id] + "_mean").c_str(), Par.mean_fit_name.c_str());
   TF1 *sigma_fit = new TF1((Par.id_map[id] + "_sigma").c_str(), Par.sigma_fit_name.c_str());
   
   std::string input_name = "../par/MS/" + Par.run_name + 
      "/" + det_name + "/" + Par.id_map[id] + ".txt";

   CheckInputFile(input_name);
   ifstream input_file(input_name.c_str());

   double par;
   int ipar = 0;

   if (Par.print_info) 
   {
      PrintInfo("From file " + input_name + ":");
      std::cout << "   ";
   }

   while (input_file >> par)
   {
      mean_fit->SetParameter(ipar, par);
      sigma_fit->SetParameter(ipar, par);
      ipar++;
      if (Par.print_info) std::cout << par << " ";
   }
   
   if (Par.print_info) std::cout << std::endl;

   MS.mean.insert({id, mean_fit});
   MS.sigma.insert({id, sigma_fit});
}

void AddM2Graph(const int cnum, std::string charge_name, std::string det_name, std::string det_hname, Color_t color, Style_t style)
{
   std::string real_spectra_file_name = "../ext/Spectra/" + 
      Par.col_system + "/" + Par.id_map[Par.part_id] + "_" +
      Par.CType.cmin_name[cnum] + "_" + Par.CType.cmax_name[cnum] +
      ".txt";

   std::string reg_eff_file_name = "../data/Efficiency/" + 
      Par.run_name + "/Reg/" + det_name + 
      "/" + Par.id_map[Par.part_id] + Par.magf + ".txt";

   std::string raw_file_name = "../data/RawYields/" + 
      Par.run_name + "/" + det_name + "/" + 
      Par.CType.cmin_name[cnum] + "_" + Par.CType.cmax_name[cnum] + 
      "_" + Par.id_map[Par.part_id] + Par.magf + ".txt";
   
   CheckInputFile(raw_file_name);
   CheckInputFile(real_spectra_file_name);
   CheckInputFile(reg_eff_file_name);
   
   ifstream raw_file(raw_file_name.c_str());
   ifstream real_spectra_file(real_spectra_file_name.c_str());
   ifstream reg_eff_file(reg_eff_file_name.c_str());
   
   TH3F *mass3d = (TH3F *) Par.input_file->Get(
      (det_hname + "_" + charge_name).c_str());
   TH1F *centr = (TH1F *) Par.input_file->Get("central_bin");

   const double nevents = centr->Integral(
      centr->GetXaxis()->FindBin(Par.CType.cmin[cnum]*10. + 0.01), 
      centr->GetXaxis()->FindBin((Par.CType.cmax[cnum]+1.)*10. - 0.01));

   Par.m2VGr.push_back(new TGraphErrors());
   
   Par.m2VGr.back()->SetLineWidth(2);
   Par.m2VGr.back()->SetLineColor(color);
   Par.m2VGr.back()->SetMarkerSize(2);
   Par.m2VGr.back()->SetMarkerColor(color);
   Par.m2VGr.back()->SetMarkerStyle(style);

   Par.VFunc.push_back(new TF1(("func" + to_string(cnum)).c_str(), "[0]"));

   Par.VFunc.back()->SetLineColor(color);
   Par.VFunc.back()->SetLineStyle(2);
   Par.VFunc.back()->SetLineWidth(3);

   if (Par.print_info) PrintInfo("Number of events in " + 
      Par.CType.cname[cnum] + " centrality class is " + to_string((long) nevents));
   
   double pt, real_yield, real_err;
   double reg_pt, reg_eff, match_err;
   double raw_pt, raw_yield, raw_err, raw_yield_veto, raw_yield_err_veto;

   if (Par.print_info) PrintInfo("From file " + real_spectra_file_name + ":");
   for (int i = 0; i < Par.ptmin.size(); i++)
   {
      real_spectra_file >> pt >> real_yield >> real_err;
      reg_eff_file >> reg_pt >> reg_eff >> match_err;
      raw_file >> raw_pt >> raw_yield >> raw_err >> raw_yield_veto >> raw_yield_err_veto;

      if (pt < Par.minimum_pt) continue;
      if (pt > Par.maximum_pt) break;
      
      if (real_spectra_file.peek() == EOF) break;
      if (reg_eff_file.peek() == EOF) break;
      if (raw_file.peek() == EOF) break;
      
      if (Average(Par.ptmin[i], Par.ptmax[i]) > pt + 0.001 ||
         Average(Par.ptmin[i], Par.ptmax[i]) < pt - 0.001) PrintError(
         "pt is not equal to the one of the real spectra: " + 
         to_string(pt) + " vs " + to_string(Average(Par.ptmin[i], Par.ptmax[i])));
      
      if (Average(Par.ptmin[i], Par.ptmax[i]) > reg_pt + 0.001 ||
         Average(Par.ptmin[i], Par.ptmax[i]) < reg_pt - 0.001) PrintError(
         "pt is not equal to the one of the matching efficiency: " +
         to_string(reg_pt) + " vs " + to_string(Average(Par.ptmin[i], Par.ptmax[i])));

      if ((Average(Par.ptmin[i], Par.ptmax[i]) > raw_pt + 0.001 ||
         Average(Par.ptmin[i], Par.ptmax[i]) < raw_pt - 0.001)) PrintError(
         "pt is not equal to the one of the raw file: " +
         to_string(raw_pt) + " vs " + to_string(Average(Par.ptmin[i], Par.ptmax[i])));
   
      TH1F *mass1d = (TH1F*) mass3d->ProjectionY("mass1d",
         mass3d->GetXaxis()->FindBin(Par.ptmin[i]+0.01), 
         mass3d->GetXaxis()->FindBin(Par.ptmax[i]-0.01), 
         mass3d->GetZaxis()->FindBin(Par.CType.cmin[cnum]*10. + 0.1),
         mass3d->GetZaxis()->FindBin((Par.CType.cmax[cnum]+1.)*10-0.1));
      
      const double spectra_norm = nevents*2*TMath::Pi()*pt*(Par.ptmax[i] - Par.ptmin[i]);
      double eff = raw_yield/reg_eff/spectra_norm/real_yield/Par.emb*Par.addit_mult[i];

      Par.ymin = Minimum(Par.ymin, eff);
      Par.ymax = Maximum(Par.ymax, eff);
      
      const double eff_err = sqrt(pow(match_err, 2) + 
         pow(raw_yield_err_veto/raw_yield_veto, 2) + 
         pow(real_err/real_yield, 2))*eff;

      if (eff <= 0.) continue;

      Par.m2VGr.back()->AddPoint(pt, eff);
      Par.m2VGr.back()->SetPointError(Par.m2VGr.back()->GetN()-1, 0., eff_err);

      Par.xmin = Minimum(Par.xmin, Par.ptmin[i]);
      Par.xmax = Maximum(Par.xmax, Par.ptmax[i]);
   }

   Par.m2VGr.back()->Fit(Par.VFunc.back(), "QMN");

   Par.legend.AddEntry(Par.m2VGr.back(), (Par.CType.cmin_name[cnum] + 
      "-" + Par.CType.cmax_name[cnum] + "\%").c_str(), "P");
}

void AddEidGraph(const int cnum, std::string charge_name, std::string det_name, std::string det_hname, Color_t color, Style_t style)
{
   std::string real_spectra_file_name = "../ext/Spectra/" + 
      Par.col_system + "/" + Par.id_map[Par.part_id] + "_" +
      Par.CType.cmin_name[cnum] + "_" + Par.CType.cmax_name[cnum] +
      ".txt";

   std::string reg_eff_file_name = "../data/Efficiency/" + 
      Par.run_name + "/Reg/" + det_name + 
      "/" + Par.id_map[Par.part_id] + Par.magf + ".txt";

   std::string raw_file_name = "../data/RawYields/" + 
      Par.run_name + "/" + det_name + "/" + 
      Par.CType.cmin_name[cnum] + "_" + Par.CType.cmax_name[cnum] + 
      "_" + Par.id_map[Par.part_id] + Par.magf + ".txt";

   ifstream raw_file(raw_file_name.c_str());
   ifstream real_spectra_file(real_spectra_file_name.c_str());
   ifstream reg_eff_file(reg_eff_file_name.c_str());
   
   TH3F *mass3d = (TH3F *) Par.input_file->Get(
      (det_hname + "_" + charge_name).c_str());
   TH1F *centr = (TH1F *) Par.input_file->Get("central_bin");

   const double nevents = centr->Integral(
      centr->GetXaxis()->FindBin(Par.CType.cmin[cnum]*10. + 0.01), 
      centr->GetXaxis()->FindBin((Par.CType.cmax[cnum]+1.)*10. - 0.01));

   Par.RatioVGr.push_back(new TGraphErrors());
   Par.RatioVGrSys.push_back(new TGraphErrors());

   Par.RatioVGr.back()->SetLineWidth(2);
   Par.RatioVGr.back()->SetLineColor(color);
   Par.RatioVGr.back()->SetMarkerSize(2);
   Par.RatioVGr.back()->SetMarkerColor(color);
   Par.RatioVGr.back()->SetMarkerStyle(style);

   Par.RatioVGrSys.back()->SetLineWidth(1);
   Par.RatioVGrSys.back()->SetLineColorAlpha(color, 0.5);
   Par.RatioVGrSys.back()->SetFillStyle(1001);
   Par.RatioVGrSys.back()->SetFillColorAlpha(color, 0.3);

   double pt, real_yield, real_err;
   double reg_pt, reg_eff, match_err;
   double raw_pt, raw_yield, raw_err, raw_yield_veto, raw_yield_err_veto;
   
   for (int i = 0; i < Par.m2VGr.back()->GetN(); i++)
   {
      real_spectra_file >> pt >> real_yield >> real_err;
      reg_eff_file >> reg_pt >> reg_eff >> match_err;
      raw_file >> raw_pt >> raw_yield >> raw_err >> raw_yield_veto >> raw_yield_err_veto;

      double raw_yield_err;
      
      const double spectra_norm = nevents*2*TMath::Pi()*pt*(Par.ptmax[i] - Par.ptmin[i]);

      const double corrected_yield = raw_yield_veto/
         reg_eff/spectra_norm*
         Par.addit_mult[i]/
         Par.average_eff/
         Par.emb;
   
      //for double check
      //const double corrected_yield = raw_yield/spectra_norm/reg_eff/Par.average_eff/Par.emb*Par.addit_mult[i];
      
      double ratio = real_yield/corrected_yield;
      double ratio_err = sqrt(pow(match_err, 2) + pow(raw_yield_err_veto/raw_yield_veto, 2));
   
      ratio_err *= ratio;

      //if (i == Par.m2VGr.back()->GetN() - 1) Print(det_name, Par.part_id, 1. - ratio - ErrPropagation(ratio_err, Par.other_sys_no_TOFw*ratio));

      Par.RatioVGr.back()->AddPoint(pt, ratio);
      Par.RatioVGr.back()->SetPointError(Par.RatioVGr.back()->GetN()-1, 0., ratio_err);
      Par.RatioVGrSys.back()->AddPoint(pt, ratio);
      Par.RatioVGrSys.back()->SetPointError(Par.RatioVGr.back()->GetN()-1, 0.04, Par.other_sys_no_TOFw*ratio);
   }
}

void SingleM2Eff(std::string det_name, std::string det_hname)
{   
   SetPar(det_name);
   int charge;
   if (Par.part_id > 0) charge = 1;
   else charge = -1;

   std::string charge_name;
   if (charge == 1) charge_name = "pos";
   else charge_name = "neg";
   
   SetMSFunc(211*charge, det_name);
   SetMSFunc(321*charge, det_name);
   SetMSFunc(2212*charge, det_name);

   /*
   for (int i = 0; i < Par.CType.cmin.size(); i++)
   {
      AddM2Graph(i, charge_name, Par.CType.color[i]-3, Par.CType.marker_style[i]);
   }
   */
   
   AddM2Graph(4, charge_name, det_name, det_hname, kRed-3, 55);

   TH1F hr = TH1F(("hr" + Par.magf + to_string(Par.part_id)).c_str(), 
      (" #epsilon_{m^{2}} of " + Par.id_map_tex[Par.part_id] + " in the " + Par.magf + " field").c_str(), 
      10, Par.xmin-0.05, Par.xmax+0.05);

   hr.SetMinimum(0.);
   hr.SetMaximum(1.05);

   hr.GetXaxis()->SetTitle("p_{T}, GeV/c");
   hr.GetYaxis()->SetTitle("#epsilon_{m^{2}}");

   hr.SetTitleSize(0.06, "X");
   hr.SetTitleSize(0.07, "Y");

   hr.GetXaxis()->SetLabelSize(0.06);
   hr.GetYaxis()->SetLabelSize(0.06);

   hr.GetYaxis()->SetTitleOffset(0.9);

   TLine line = TLine(Par.xmin-0.05, 1., Par.xmax+0.05, 1.);
   line.SetLineStyle(2);
   line.SetLineWidth(3);
   line.SetLineColor(kGray+1);

   TLine eid_line = TLine(Par.xmin-0.05, erf(Par.srange/sqrt(2.)), Par.xmax+0.05, erf(Par.srange/sqrt(2.)));
   eid_line.SetLineStyle(2);
   eid_line.SetLineWidth(3);
   eid_line.SetLineColor(kGray+1);

   TCanvas m2canv = TCanvas("m2canv", "m2canv", 600, 350);
   
   gPad->SetLeftMargin(0.13);
   gPad->SetBottomMargin(0.135);

   hr.Draw();
   hr.Draw("SAME AXIS X+ Y+");

   line.Draw();

   for (TGraphErrors *gr : Par.m2VGr) gr->Draw("P");
   for (TF1 *func : Par.VFunc)
   {
      func->SetRange(Par.xmin-0.025, Par.xmax+0.025);
      func->Draw("SAME");
   }

   Par.legend.SetLineColorAlpha(0., 0.);
   Par.legend.SetFillColorAlpha(0., 0.);

   //Par.legend.Draw();

   const double sys = (Par.ymax - Par.ymin)/2.;
   Par.average_eff = Par.VFunc.back()->GetParameter(0);
   Par.average_eff_sys = ErrPropagation(sys/Par.VFunc.back()->GetParameter(0), Par.other_sys);

   Par.other_sys_no_TOFw = ErrPropagation(Par.other_sys_no_TOFw, sys/Par.VFunc.back()->GetParameter(0));
   
   //for (TF1 * func : Par.VFunc) Print(det_name, Par.part_id, Par.magf, Par.average_eff, sys, Par.average_eff_sys);

   AddEidGraph(4, charge_name, det_name, det_hname, kRed-3, 55);

   PrintCanvas(&m2canv, "../output/Efficiency/" + Par.run_name + 
      "/M2Eff_" + det_name + "_" + Par.id_map[Par.part_id] + Par.magf);
   
   TCanvas eid_canv = TCanvas("eid_canv", "eid_canv", 600, 350);
   
   gPad->SetLeftMargin(0.16);
   gPad->SetBottomMargin(0.13);

   TH1 *rh = gPad->DrawFrame(Par.xmin-0.05, 0., Par.xmax+0.05, 2.);
   rh->Draw("SAME AXIS X+ Y+");
   rh->SetTitle((Par.id_map_tex[Par.part_id] + " in the " + Par.magf + " field").c_str());

   rh->GetXaxis()->SetLabelSize(0.06);
   rh->GetYaxis()->SetLabelSize(0.06);

   rh->SetTitleSize(0.06, "X");
   rh->SetTitleSize(0.06, "Y");

   line.Draw();
   rh->GetYaxis()->SetTitle("#frac{d^{2}N_{TOFw}/dp_{T}dy}{1/(N_{evt} #epsilon_{m^{2}} #epsilon_{reg} #epsilon_{id} #epsilon_{emb}) d^{2}Y/dp_{T} dy}");
   rh->GetXaxis()->SetTitle("p_{T}");

   for (TGraphErrors *gr : Par.RatioVGr) gr->Draw("P");   
   for (TGraphErrors *gr : Par.RatioVGrSys) gr->Draw("5");   

   PrintCanvas(&eid_canv, "../output/Efficiency/" + Par.run_name + 
      "/Eid_" + det_name + "_" + Par.id_map[Par.part_id] + Par.magf);
}

void M2Eff()
{
   ROOT::EnableImplicitMT();
   
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(kFALSE);
   gErrorIgnoreLevel = kWarning;
   
   system(("mkdir -p ../../sim_analysis/input/M2Eff/" + Par.run_name).c_str());
   system(("mkdir -p ../../sim_analysis/input/Systematics/" + Par.run_name).c_str());

   for (std::string magf : Par.magf_queue)
   {
      Print(magf);
      
      Par.magf = magf;

      std::string input_file_name = "../data/" + Par.run_name + "/sum" + Par.magf + ".root";
      CheckInputFile(input_file_name);
      Par.input_file = std::unique_ptr<TFile>(new TFile(input_file_name.c_str()));
      
      for (int i = 0; i < Par.part_id_queue.size(); i++)
      {
         Print(Par.part_id_queue[i]);
         
         Par.part_id = Par.part_id_queue[i];
         Par.other_sys = ErrPropagation(Par.tofw_sys_queue[i], 
            Par.bg_fit_sys_queue[i], Par.fit_range_sys_queue[i]);
         Par.other_sys_no_TOFw = ErrPropagation(Par.bg_fit_sys_queue[i], Par.fit_range_sys_queue[i]);
      
         for (int j = 0; j < Par.det_name_queue.size(); j++)
         {
            std::string m2eff_output_name = "../../sim_analysis/input/M2Eff/" + 
               Par.run_name + "/" + Par.id_map[Par.part_id] + "_" + Par.det_name_queue[j] + magf + ".txt";
            std::string m2eff_err_output_name = "../../sim_analysis/input/Systematics/" + 
               Par.run_name + "/m2eff_" + Par.id_map[Par.part_id] + "_" + Par.det_name_queue[j] + magf + ".txt";
            
            ofstream m2eff_output(m2eff_output_name);
            ofstream m2eff_err_output(m2eff_err_output_name);

            for (int sect = Par.det_min_sect_queue[j]; sect < 4; sect++)
            {
               Print(Par.det_name_queue[j] + to_string(sect));
               Par.ymin = 1e31;
               Par.ymax = -1e31;
                     
               SingleM2Eff(Par.det_name_queue[j] + to_string(sect), 
                  Par.det_hname_queue[j] + to_string(sect));
               
               m2eff_output << Par.average_eff << " ";
               m2eff_err_output << Par.average_eff_sys << " ";
               Print(Par.average_eff_sys*100.);
               
               Par.legend.Clear();
               Par.RatioVGr.clear();
               Par.RatioVGrSys.clear();
               Par.m2VGr.clear();
               Par.VFunc.clear();

               Par.xmin = 1e31;
               Par.xmax = -1e31;

               MS.mean.clear();
               MS.sigma.clear();
            }

            m2eff_output.close();
            m2eff_err_output.close();

            //PrintInfo("File " + m2eff_output_name + " was written");
            //PrintInfo("File " + m2eff_err_output_name + " was written");
         }
      }
   }
   Par.input_file->Close();
}
