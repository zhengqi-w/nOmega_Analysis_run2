#include <iostream>
namespace config
{ 
  const std::vector<Float_t> Cen_bins = {0, 20 ,60};
  const int nCenbins = Cen_bins.size()-1;
  const std::vector<Float_t> Pt_bins = {2, 3.5 , 5, 8};
  const int nPtbins = Pt_bins.size()-1;
  const Float_t cascade_cpa_min = 0.999;
  const Float_t v0_cpa_min = 0.995;
  const Float_t cascade_dca_max = 1;
  const Float_t v0_dca_max = 1;
  const Float_t abs_kaon_nsigmatpc_max = 3;
  const Float_t abs_deuteron_nsigmatpc_max = 3;
  const Float_t abs_pion_nsigmatpc_max = 3;
  const Float_t tvtopv_min = 1, tvtopv_max = 10;
  const Float_t svtopv_min = 2, svtopv_max = 10;
  const Float_t svtotv_min = 0.0, svtotv_max = 10;
  const Float_t kaon_tpcncls_min = 70;
  const Float_t deuteron_tpcncls_min = 80;
  const Float_t pion_tpcncls_min = 70;
  const Float_t kaon_dcatopv_min = 0.05;
  const Float_t deuteron_dcatopv_min = 0.2;
  const Float_t pion_dcatopv_min = 0.2;

  
  //bins for topological cuts
  //std::vector <double> 
}