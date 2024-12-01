class DataSelector
{
   public:

   DataSelector(const std::string &datasetName)
   {
      if (datasetName == "Run7AuAu200")
      {
         using namespace Run7AuAu200;
         fIsDeadDC = &IsDeadDC;
         fIsDeadPC1 = &IsDeadPC1;
         fIsDeadPC2 = &IsDeadPC2;
         fIsDeadPC3 = &IsDeadPC3;
         fIsDeadEMCal = &IsDeadEMCal;
         fIsBadSlat = &IsBadSlat;
         fIsBadStripTOFw = &IsBadStripTOFw;
         fIsDeadTOFe = &IsDeadTOFe;
         fIsDeadTOFw = &IsDeadTOFw;
      }
   }

   bool IsDeadDC(const double phi, const double zed, const double board, const double alpha)
   {
      return fIsDeadDC(phi, zed, board, alpha);
   }
   
   private:
   
   bool (*fIsDeadDC)(const double, const double, const double, const double);
   bool (*fIsDeadPC1)(const double, const double, const double);
   bool (*fIsDeadPC2)(const double, const double);
   bool (*fIsDeadEMCal)(const double, const double, const int, const double, const double);
   bool (*fIsBadSlat)(const int);
   bool (*fIsBadStripTOFw)(const int);
   bool (*fIsDeadTOFe)(const double, const double, const double);
   bool (*fIsDeadTOFw)(const double, const double, const double);
};
