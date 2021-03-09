/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(nnp2v,PairNNP2v)

#else

#ifndef LMP_PAIR_NNP2V_H
#define LMP_PAIR_NNP2V_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNNP2v : public Pair {
 public:
  PairNNP2v(class LAMMPS *);
  virtual ~PairNNP2v();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);

 private:
  bigint num_atoms;
  bigint num_atoms_a,num_atoms_b;

 protected:
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;              // added for TersoffZBL
    double ZBLcut,ZBLexpscale;
    double c5,ca1,ca4;           // added for TersoffMOD
    double powern_del;
    double c0;                   // added for TersoffMODC
  };

  Param *params;                // parameter set for an I-J-K interaction
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  double rc;
  int nweight_a,nweight_b;
  double *weight_a,*weight_b;
  int *g_a,*g_b;
  int nlayer_a,nlayer_b;
  int *layer_a,*layer_b;
  //(natom_x,num_g2_x_x)
  double *g2_a_a, *g2_a_b;
  double *g5_a_aa,*g5_a_ab,
         *g5_a_bb;
  double *g2_b_a, *g2_b_b;
  double *g5_b_aa,*g5_b_ab,
         *g5_b_bb;
  //(natom_x,natom,3,num_g2_x_x)
  double ***g2_deriv_a_a, ***g2_deriv_a_b;
  double ***g5_deriv_a_aa,***g5_deriv_a_ab,
         ***g5_deriv_a_bb;
  double ***g2_deriv_b_a, ***g2_deriv_b_b;
  double ***g5_deriv_b_aa,***g5_deriv_b_ab,
         ***g5_deriv_b_bb;
  double *g2_a_a_max, *g2_a_b_max;
  double *g5_a_aa_max,*g5_a_ab_max,
         *g5_a_bb_max;
  double *g2_b_a_max, *g2_b_b_max;
  double *g5_b_aa_max,*g5_b_ab_max,
         *g5_b_bb_max;
  double *g2_a_a_min, *g2_a_b_min;
  double *g5_a_aa_min,*g5_a_ab_min,
         *g5_a_bb_min;
  double *g2_b_a_min, *g2_b_b_min;
  double *g5_b_aa_min,*g5_b_ab_min,
         *g5_b_bb_min;
  double *eta2_a_a,*eta2_a_b;
  double *eta2_b_a,*eta2_b_b;
  double *rs2_a_a,*rs2_a_b;
  double *rs2_b_a,*rs2_b_b;
  double *eta5_a_aa,*eta5_a_ab,
         *eta5_a_bb;
  double *eta5_b_aa,*eta5_b_ab,
         *eta5_b_bb;
  double *theta5_a_aa,*theta5_a_ab,
         *theta5_a_bb;
  double *theta5_b_aa,*theta5_b_ab,
         *theta5_b_bb;
  double *zeta5_a_aa,*zeta5_a_ab,
         *zeta5_a_bb;
  double *zeta5_b_aa,*zeta5_b_ab,
         *zeta5_b_bb;
  double *lambda5_a_aa,*lambda5_a_ab,
         *lambda5_a_bb;
  double *lambda5_b_aa,*lambda5_b_ab,
         *lambda5_b_bb;
  double *input_a,*input_b;
  double *dGdx_a,*dGdy_a,*dGdz_a;
  double *dGdx_b,*dGdy_b,*dGdz_b;
  int    nhidden_a,nhidden_b;
  double *hidden_a,*hidden_b;
  //double *nnp_eng,*f_x,*f_y,*f_z,*f_s_x,*f_s_y,*f_s_z;//**riki;
  //double *x_tmp_x,*x_tmp_y,*x_tmp_z,*x_s_x,*x_s_y,*x_s_z;

  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  void arraynnp();
  void sf_g2(int, int, double,
             double, double, double, 
             double, double, double,
             double *&, double ***&, int,
             double *&, double *& );
  void sf_g5(int, int, int, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double *&, double ***&, int,
             double *&, double *&, double *& );
  void hdnnp_e_a(double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, double & );
  void hdnnp_f_a(double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 int, double *&, double &, double &, double & );
  void hdnnp_e_b(double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, double & );
  void hdnnp_f_b(double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 int, double *&, double &, double &, double & );

  virtual double fc(double);
  virtual double fc_deriv(double,double,double);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Tersoff requires atom IDs

This is a requirement to use the Tersoff potential.

E: Pair style Tersoff requires newton pair on

See the newton command.  This is a restriction to use the Tersoff
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/
