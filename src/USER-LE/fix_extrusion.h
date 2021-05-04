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

#ifdef FIX_CLASS

FixStyle(extrusion,FixExtrusion)

#else

#ifndef LMP_FIX_EXTRUSION_H
#define LMP_FIX_EXTRUSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExtrusion : public Fix {
 public:
  FixExtrusion(class LAMMPS *, int, char **);
  ~FixExtrusion();
  int setmask();
  void init();
  void setup(int);
  void post_integrate();
  void post_integrate_respa(int,int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me,nprocs;
  int btype,seed;
  double cutoff,cutsq,fraction;
  bigint lastcheck;

  int breakcount,breakcounttotal;
  int createcount,createcounttotal;//to track balance between broken and created bonds
  int nmax;
  int chain_length;
  int *bondcount;
  tagint *to_remove, *to_add, *recreate;
  tagint *final_to_remove, *final_to_add;
  double *distsq_c,*distsq_b,*probability;

  int nbreak,ncreate,maxbreak, maxcreate;//, nshift, maxshift;//number of created and broken bonds should be equal
  //tagint **shifted;
  tagint **broken, **created;
  tagint *copy;

  class RanMars *random;
  int nlevels_respa;

  int countflag, commflag;
  int nangles,ndihedrals,nimpropers;

  //unique for extrusion
  int neutral_type, ctcf_left, ctcf_right; //already have btype
  double through_prob;

  void check_ghosts();
  void update_topology();
  void rebuild_special_one(int);
  int dedup(int, int, tagint *);

  // DEBUG

  void print_bb();
  void print_copy(const char *, tagint, int, int, int, int *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid bond type in fix bond/break command

Self-explanatory.

E: Cannot use fix bond/break with non-molecular systems

Only systems with bonds that can be changed can be used.  Atom_style
template does not qualify.

E: Cannot yet use fix bond/break with this improper style

This is a current restriction in LAMMPS.

E: Fix bond/break needs ghost atoms from further away

This is because the fix needs to walk bonds to a certain distance to
acquire needed info, The comm_modify cutoff command can be used to
extend the communication range.

*/
