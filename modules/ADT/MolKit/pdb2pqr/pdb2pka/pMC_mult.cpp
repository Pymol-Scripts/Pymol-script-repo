//
// $Id: pMC_mult.cpp,v 1.1 2007/08/03 20:59:24 sargis Exp $
//
// (C) Jens Erik Nielsen, University College Dublin 2005
//
#include "pMC_mult.h"
//
//
//

void MC::reformat_arrays() {
  //
  // Reformat the matrix
  //
  printf ("Reformatting arrays\n");
  _groups=static_cast<int>(_num_states.size());
  printf ("Number of groups: %d\n",_groups);
  for (int group=0;group<_groups;group++) {
    printf ("Num_states for group %d is %d\n",group,_num_states[group]); 
  }

  //
  int count=0;
  for (int row=0;row<_groups;row++) {
    printf ("Constructing matrix for group %d\n",row);
    vector<vector<vector<double > > > row_vals;
    for (int column=0;column<_groups;column++) {
      printf ("Second level matrix for group %d\n",column);
      vector<vector<double > > column_vals;
      for (int group1_states=0;group1_states<_num_states[row];group1_states++) {
	vector<double> group1_s;
	for (int group2_states=0;group2_states<_num_states[column];group2_states++) {
	  group1_s.push_back(_lin_matrix[count]);
	 
	  printf ("While reformatting arrays: g1: %d, g2: %d, st1: %d st2: %d value: %5.3f\n",
		  row,column,group1_states,group2_states,_lin_matrix[count]);
	  count=count+1;
	}
	column_vals.push_back(group1_s);
      }
      row_vals.push_back(column_vals);
    }
    _matrix.push_back(row_vals);
  }
  //
  // Reformat the intrinsic pKa array and the charged_state_array
  //
  count=0;
  for (int group=0;group<_groups;group++) {
    vector<double> these_intpKas;
    vector<int> these_charged_state;
    for (int t_state=0;t_state<_num_states[group];t_state++) {
      these_intpKas.push_back(_intpKas_lin[count]);
      these_charged_state.push_back(_charged_state_lin[count]);
      count++;
    }
    _intpKas.push_back(these_intpKas);
    _charged_state.push_back(these_charged_state);
  }
  //
  // Set natural log
  //
  lnten=log(10);
  return;
}

//
// ---------------------
//

vector<float> MC::calc_pKas(float pH_start,float pH_end, float pH_step) {
  //
  // Calculate pKa values for the system
  //
  // First get charges at all pH values
  //
  float max_pH=0.0;
  vector< vector<float> > charges;
  for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
    charges.push_back(calc_charge(pH));
    max_pH=pH;
  }
  //
  // Print the intrinsic pKa values we got
  //
  for (int group=0;group<_groups;group++) {
    for (int state=0;state<_num_states[group];state++) {
      printf ("Intrinsic pKa group %d state %d is %5.2f\n",group,state,_intpKas[group][state]);
    }
  }
  // 
  // Now determine pKa values
  //
  int datapoints=11; // How many data points do we use for pKa determination?
  datapoints=(datapoints-1)/2;
  vector<float> pKas;
  for (int group=0;group<_groups;group++) {
    int count=0;
    float pKa=-999.9;
    float last_crg=charges[count][group];
    //
    for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
      if ((pH-max_pH)>0.0) {
	continue;
      }
      float this_crg=charges[count][group];
      if (_acid_base[group]==1.0) {
	if (this_crg<=0.5 && last_crg>0.5) {
	  //pKa=(last_crg-0.5)/(last_crg-this_crg)*pH_step+(pH-pH_step);
	  //
	  // Get ph,charge sets and calc pKa from those
	  //
	  vector<double> pHs_pKadet;
	  vector<float> charges_pKadet;
	  int count2=count-static_cast<int>(datapoints);
	  if (count2<0) {count2=0;}
	  for (double pH2=max(pH_start,pH-datapoints*pH_step);pH2<min(pH_end,pH+datapoints*pH_step);pH2=pH2+pH_step) {
	    pHs_pKadet.push_back(pH2);
	    charges_pKadet.push_back(charges[count2][group]);
	    count2=count2+1;
	  }
	  pKa=calc_pKa(charges_pKadet,pHs_pKadet,_acid_base[group]);
	}
      } else {
	if (this_crg<=-0.5 && last_crg>-0.5) {
	  //pKa=(last_crg-(-0.5))/(last_crg-this_crg)*pH_step+(pH-pH_step);
	  //
	  // Get ph,charge sets and calc pKa from those
	  //
	  vector<double> pHs_pKadet;
	  vector<float> charges_pKadet;
	  int count2=count-static_cast<int>(datapoints);
	  if (count2<0) {count2=0;}
	  for (double pH2=max(pH_start,pH-datapoints*pH_step);pH2<min(pH_end,pH+datapoints*pH_step);pH2=pH2+pH_step) {
	    pHs_pKadet.push_back(pH2);
	    charges_pKadet.push_back(charges[count2][group]);
	    count2=count2+1;
	  }
	  pKa=calc_pKa(charges_pKadet,pHs_pKadet,_acid_base[group]);
	}
      }
      last_crg=this_crg;
      count=count+1;
    }
    pKas.push_back(pKa);
  }   
  //
  // We also return all the charges
  //
  // The first numbers we return is the start pH, the pH step, and the number of values
  //
  int num_pHs=0;
  for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
    num_pHs++;
  }
  pKas.push_back(pH_start);
  pKas.push_back(pH_step);
  pKas.push_back(static_cast<float>(num_pHs));
  // 
  // Now add the charges
  //
  float this_crg;
  int count=0;
  for (int group=0;group<_groups;group++) {
    count=0;
    for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
      pKas.push_back(pH);
      this_crg=charges[count][group];
      pKas.push_back(this_crg);
      count=count+1;
    }
    pKas.push_back(999.0);
    pKas.push_back(-999.0);
  }
  return pKas;
}

//
// ---------------------
//

double MC::calc_pKa(vector<float> charges,vector<double> pHs,double acid_base) {
  //
  // Calculate the pKa value from a selection of charges and pH values
  //
  // Assume perfect Henderson-Hasselbalch behaviour
  //
  // acid_base = -1.0 for acids
  // acid_base = 1.0 for bases
  //
  double ratio=0.0;
  vector<double> pKas;
  double pKa=0.0;
  int points=static_cast<int>(charges.size());
  for (int count=0;count<points;count++) {
    if (acid_base!=1.0) {
      ratio=fabs(charges[count])/(1.0-fabs(charges[count]));
    } else {
      ratio=(1.0-fabs(charges[count]))/fabs(charges[count]);
    }
    pKas.push_back(pHs[count]-log10(ratio));
  }
  //
  // Find the average of the pKa values
  //
  double sum=0.0;
  for (int count=0;count<static_cast<int>(pKas.size());count++) {
    sum=sum+pKas[count];
  }
  //printf ("Done with calc_pKa\n");
  pKa=sum/static_cast<double>(pKas.size());
  return pKa;
}
    

//
// ---------------------
//

vector<float> MC::calc_charge(float pH) {
  //
  // Calculate the fractional charges at this pH
  //
  // Initialise random number generator
  //
  srand(time(NULL));
  //
  // Get a random starting state
  //
  vector<int> current_state;
  vector<int> try_state;
  vector<vector<int> > sum_state;
  for (int group=0;group<_groups;group++) {
    current_state.push_back(static_cast<int>(rand()%2));
    if (current_state[group]==2) {
      current_state[group]=1;
    }
    // 
    // Dummy initialisation of try_state and sum_state
    //
    try_state.push_back(0);
    //sum_state.push_back(0);
    vector<int> dummy;
    sum_state.push_back(dummy);
  }
  //
  // Get the energy of the starting state 
  //
  double current_energy=get_energy(pH,current_state);
  //
  // Count the number of charged and uncharged states
  //
  vector <int> charged_states;
  vector <int> uncharged_states;
  for (int group=0;group<_groups;group++) {
    int cha=0;
    int uncha=0;
    for (int state=0;state<_num_states[group];state++) {
      if (_charged_state[group][state]) {
	cha++;
      } else {
	uncha++;
      }
    }
    charged_states.push_back(cha);
    uncharged_states.push_back(uncha);
  }
  // Need to incorporate direct switching between neutral states.. no tie now
  //
  //
  // Start the MC loop
  //
  int eqsteps=_MCsteps/10;
  int keep=0;
  double tilf=0.0;
  double try_energy_new=0.0;
  double diff=0.0;
  for (int step=0;step<_MCsteps;step++) {
    //
    // Copy the current state to trystate
    //
    for (int count=0;count<_groups;count++) {
      try_state[count]=current_state[count];
    }
    //
    // Change a random group
    //
    int rand_group=static_cast<int>(rand()%_groups);
    //
    // Change to a random state
    //
    // Current charge state
    int cur_charge_state=_charged_state[rand_group][current_state[rand_group]];
    //
    //printf ("Group %d, current state: %d with charge %d\n",rand_group,current_state[rand_group],cur_charge_state);
    int rand_group_state=static_cast<int>(rand()%_num_states[rand_group]);
    while (_charged_state[rand_group][rand_group_state]==cur_charge_state) {
      //
      rand_group_state=static_cast<int>(rand()%_num_states[rand_group]);
      //printf ("New charged state: %d with charge %d. Numstates is %d \n",rand_group_state,
      //      _charged_state[rand_group][rand_group_state],_num_states[rand_group]);
    }
    //printf ("Changing group %d to state %d\n\n",rand_group,rand_group_state);
    // Change state
    try_state[rand_group]=abs(rand_group_state);
    //
    // Get the energy of the new state
    //
    //try_energy_new=get_energy_fast(pH,try_state,rand_group,current_energy);
    try_energy_new=get_energy(pH,try_state);
    //
    // Keep or reject?
    //
    diff=try_energy_new-current_energy;
    if (diff<0.0) {
      //
      // Keep
      //
      for (int count=0;count<_groups;count++) {
	current_state[count]=try_state[count];
      }
      current_energy=try_energy_new;
    } else {
      if (diff<20.0) {
	tilf=static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)+static_cast<double>(1));
	if (tilf<exp(-diff)) {
	  //
	  // Keep
	  //
	  for (int count=0;count<_groups;count++) {
	    current_state[count]=try_state[count];
	  }
	  current_energy=try_energy_new;
	}
      }
    }
    //
    // Record the state if we have equilibrated
    //
    if (step>eqsteps) {
      for (int count=0;count<_groups;count++) {
	sum_state[count].push_back(current_state[count]);
	//sum_state[count]=sum_state[count]+current_state[count];
      }
    }
  }
  //
  // Calculate fractional charge
  //
  int sample_steps=_MCsteps-eqsteps;
  vector<float> charges_thispH;
  for (int count=0;count<_groups;count++){
    float charge=0.0;
    charge=static_cast<float>(_acid_base[count]);
    //
    // Loop over all the steps and add the correct charge
    //
    float charge_sum=0.0;
    for (int step=0;step<sum_state[count].size();step++){
      //
      // Get the charge from the charged state
      //
      int state=sum_state[count][step];
      charge_sum=charge_sum+charge*static_cast<float>(_charged_state[count][state]);
    }
    charge=charge_sum/(static_cast<float>(sample_steps));
    charges_thispH.push_back(charge);
    printf ("pH: %5.2f, Group: %d, Charge: %5.3f\n",pH,count,charge);
  }
  return charges_thispH;
}

//
// --------------------
// 

double MC::get_energy(float pH,vector<int> state) {
  //
  // Calculate the energy of the present state
  //
  //printf ("\nCalculating energy\n");
  double pH_value=static_cast<double>(pH);
  double energy=0.0;
  for (int group1=0;group1<_groups;group1++) {
    //
    // Add the energy from the intrinsic pKa
    //
    //printf ("State in get_energy for group %d is %d\n",group1,state[group1]);
    int charge_grp1=_charged_state[group1][state[group1]];
    if (charge_grp1!=0) {
    //printf ("Intrinsic pKa for group: %d in state %d is %5.2f\n",
    // 	      group1,state[group1],
    //       _intpKas[group1][state[group1]]);
      energy=energy+_acid_base[group1]*lnten*(pH_value-_intpKas[group1][state[group1]]);
    } else {
      energy=energy+_intpKas[group1][state[group1]]*lnten;
    }
      //
      // Entropy correction
      //
      //if (_charged_state[group1][state[group1]]==0 && _acid_base[group1]==-1) {
      //energy=energy+log(static_cast<double>(_charged_state[group1].size())-1)*1.80;
	//printf("Entropy correction: states-1 %d  correction: %5.3f \n",
	//      (_charged_state[group1].size())-1,
	//      log(static_cast<double>(_charged_state[group1].size())-1)*1.80 );
      //}
      //
      // Add the charged-charged energies
      //
      //}
    if (charge_grp1!=0) {
      if (charge_grp1!=0) { //Is this a charged state?
	for (int group2=0;group2<_groups;group2++) {
	  int charge_grp2=_charged_state[group2][state[group2]];
	  if (charge_grp2!=0 and group2!=group1) {
	    // printf ("For group: %d Trying to get matrix: %d %d %d %d: %5.3f\n",
// 		    group1,
// 		    group1,
// 		    group2,
// 		    state[group1],
// 		    state[group2],
// 		    _matrix[group1][group2][state[group1]][state[group2]]/2.0);
	    energy=energy+_matrix[group1][group2][state[group1]][state[group2]]/2.0;
	  }
	}
      }
    }
  }
  //printf ("Exiting Energy calcs\n");
  return energy;
}
	
//
// --------------------
// 
  
// double MC::get_energy_fast(float pH,vector<int> state,int change_group,double old_energy) {
//   //
//   // Calculate the energy of the present state
//   //
//   int startpointer;
//   double energy=old_energy;
//   double energy_diff=0.0;
//   //
//   // Add the energy from the intrinsic pKa
//   //
//   energy_diff=_acid_base[change_group]*lnten*(pH-_intpKas[change_group]);
//   //
//   // Add the charged-charged energies
//   //
//   for (int group2=0;group2<_groups;group2++) {
//     if (state[group2]==1 && group2!=change_group) {
//       energy_diff=energy_diff+_matrix[change_group][group2];
//     }
//   }
//   //
//   // Should we add or subtract the energy
//   //
//   if (state[change_group]==1) {
//     //
//     // Group became charged - add energy
//     //
//     energy=energy+energy_diff;
//   } else {
//     // 
//     // Group became uncharged
//     //
//     energy=energy-energy_diff;
//   }
//   return energy;
// }
