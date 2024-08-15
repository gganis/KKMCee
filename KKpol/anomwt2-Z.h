
  // Interface to code provided by Z. Was for approximate evaluation of Tau polarization in e+e- -> tau+tau- production
  // Included in TauPair.cxx, used by TauPair::GetPol
  // Author G Ganis, CERN, Aug 2024
  extern "C"{
      // /BornV/ common block filler
      void fillbornv_(int *keyelw, double *swsq, double *alfinv, double *mz, double *gammz, double GSW[]);
      // Main routine
      void dipolqqrijradcor_(int *iqed, double *Energy, double *theta, double *ReA0, double *ImA0, double *ReB, double *ImB,
                                                                       double *ReX, double *ImX, double *ReY, double *ImY, double R[][4], int *channel);
  }
