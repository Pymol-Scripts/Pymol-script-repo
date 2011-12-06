//
// $Id: pMC_mult.h,v 1.1 2007/08/03 20:59:24 sargis Exp $
//
// (C) Jens Erik Nielsen, University College Dublin 2005
//
#ifndef PKAMC_H
#define PKAMC_H

#include <string>
#include <vector>
#include <math.h>
#include <time.h>

using namespace std;

class MC {
 public:
  //
  // Functions
  //
  MC(vector<double> intpKas, 
     vector<double> lin_matrix, 
     vector<double> ab,
     vector<int> num_states,
     vector<int>charged_state): _intpKas_lin(intpKas),_lin_matrix(lin_matrix),_acid_base(ab),_num_states(num_states),_charged_state_lin(charged_state) {
    reformat_arrays();
    // Set default value for MCsteps
    _MCsteps=20000;
  };
  //
  vector<float> calc_pKas(float pH_start,float pH_end, float pH_step);
  //
  void set_MCsteps(int MCsteps) {
    _MCsteps=MCsteps;
    return;
  }
  //
  // Private functions
  //
 private:
  // Reformats the matrix
  void reformat_arrays();
  //


  double calc_pKa(vector<float> charges, vector<double> pHs,double acid_base);
  vector<float> calc_charge(float pH);
  double get_energy(float pH, vector<int> state);
  //double get_energy_fast(float pH,vector<int> state,int change_group,double old_energy);
  //
  //
  // Variables
  //
  vector<double> _intpKas_lin, _lin_matrix, _acid_base;
  vector<vector<double> > _intpKas;
  vector<vector<vector<vector<double> > > > _matrix;
  vector<int> _charged_state_lin, _num_states;
  vector<vector<int> > _charged_state;
  int _groups, _MCsteps;
  double lnten;
  //vector<vector <double>> charges;
};
#endif
