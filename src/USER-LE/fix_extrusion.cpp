/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_extrusion.h"
#include <mpi.h>
#include <cstring>
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include <math.h>


using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20
#define DELTA 16

/* ---------------------------------------------------------------------- */

FixExtrusion::FixExtrusion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  to_remove(NULL), to_add(NULL), final_to_remove(NULL), final_to_add(NULL), 
  distsq_c(NULL), distsq_b(NULL), probability(NULL), broken(NULL), created(NULL), copy(NULL),
  random(NULL), bondcount(NULL), recreate(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix extrusion command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix extrusion command, n_steps <= 0");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  //btype = utils::inumeric(FLERR,arg[4]);
  //cutoff = utils::numeric(FLERR,arg[5],false,lmp);
  neutral_type = utils::inumeric(FLERR,arg[4],false,lmp);
  ctcf_left = utils::inumeric(FLERR,arg[5],false,lmp);
  ctcf_right = utils::inumeric(FLERR,arg[6],false,lmp);
  

  if (neutral_type < 1 || neutral_type > atom->ntypes ||
      ctcf_left < 1 || ctcf_left > atom->ntypes ||
      ctcf_right < 1 || ctcf_right > atom->ntypes)
    error->all(FLERR,"Invalid atom type (CTCF) in fix extrusion command");

  through_prob = utils::numeric(FLERR,arg[7],false,lmp);
  if (through_prob <=0 || through_prob>1)
  error->all(FLERR,"Invalid probability to pass through CTCF in fix extrusion command");

  btype = utils::inumeric(FLERR,arg[8],false,lmp);
    if (btype < 1 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid atom type in fix extrusion command");
  
  
  chain_length = utils::inumeric(FLERR,arg[9],false,lmp);
    if (chain_length < 1 || chain_length > atom->natoms)
    error->all(FLERR,"Invalid chain_length in fix extrusion command");
  
  // error check

  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use fix extrusion with non-molecular systems");

  // initialize Marsaglia RNG with processor-unique seed
  int seed = 12345;
  random = new RanMars(lmp,seed + me);

  // perform initial allocation of atom-based arrays
  // register with Atom class
  // bondcount values will be initialized in setup()

  bondcount = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  countflag = 0;

  // set comm sizes needed by this fix
  // forward is big due to comm of broken bonds and 1-2 neighbors

  comm_forward = MAX(2,2+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  to_add = to_remove = final_to_add = final_to_remove = NULL;
  recreate = NULL;
  maxcreate = 0;
  maxbreak = 0;

  // copy = special list for one atom
  // size = ms^2 + ms is sufficient
  // b/c in rebuild_special_one() neighs of all 1-2s are added,
  //   then a dedup(), then neighs of all 1-3s are added, then final dedup()
  // this means intermediate size cannot exceed ms^2 + ms

  int maxspecial = atom->maxspecial;
  if (me==0) printf ("Attention! maxspecial = %d\n", maxspecial);
  copy = new tagint[maxspecial*maxspecial + maxspecial];

  // zero out stats
  createcount = 0;
  createcounttotal = 0;
  breakcount = 0;
  breakcounttotal = 0;
  created = NULL;
  //broken = NULL;
}

/* ---------------------------------------------------------------------- */

FixExtrusion::~FixExtrusion()
{
  delete random;

  // delete locally stored arrays
  memory->destroy(bondcount);
  memory->destroy(to_remove);
  memory->destroy(to_add);
  memory->destroy(recreate);
  memory->destroy(final_to_remove);
  memory->destroy(final_to_add);
  memory->destroy(distsq_c);
  memory->destroy(distsq_b);
  memory->destroy(broken);
  memory->destroy(created);
  delete [] copy;
}

/* ---------------------------------------------------------------------- */

int FixExtrusion::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  //mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExtrusion::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // warn if more than one fix bond/create or also a fix bond/break
  // because this fix stores per-atom state in bondcount
  //   if other fixes create/break bonds, this fix will not know about it

  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"extrusion") == 0) count++;
  }
  if (count > 1 && me == 0)
    error->warning(FLERR,"Fix extrusion is used multiple times "
                   " - may not work as expected");

  lastcheck = -1;

  // DEBUG
  //print_bb();
}


void FixExtrusion::setup(int /*vflag*/)
{
  int i,j,m;

  // compute initial bondcount if this is first run
  // can't do this earlier, in constructor or init, b/c need ghost info

  if (countflag) return;
  countflag = 1;

  // count bonds stored with each bond I own
  // if newton bond is not set, just increment count on atom I
  // if newton bond is set, also increment count on atom J even if ghost
  // bondcount is long enough to tally ghost atom counts

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

  for (i = 0; i < nall; i++) bondcount[i] = 0;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] == btype) {
        bondcount[i]++;
        if (newton_bond) {
          m = atom->map(bond_atom[i][j]);
          if (m < 0)
            error->one(FLERR,"Fix extrusion needs ghost atoms "
                       "from further away");
          bondcount[m]++;
        }
      }
    }

  // if newton_bond is set, need to sum bondcount

  commflag = 1;
  if (newton_bond) comm->reverse_comm_fix(this,1);
}


 /*---------------------------------------------------------------------- */

void FixExtrusion::post_integrate()
{
  int i,j,k,m,n,i1,i2,n1,n2,n3,type;
  int rb, lb;
  double delx,dely,delz,rsq;
  tagint *slist;
  int local_left, local_right;
  bool l_flag;

  if (update->ntimestep % nevery - 1) return; // 1 to avoid traj_recording problems

  // count bonds stored with each bond I own
  // if newton bond is not set, just increment count on atom I
  // if newton bond is set, also increment count on atom J even if ghost
  // bondcount is long enough to tally ghost atom counts

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

  for (i = 0; i < nall; i++) bondcount[i] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] == btype) {
        bondcount[i]++;
        if (bondcount[i] > 1) error->one(FLERR,"Fix extrusion, more than one bond type 2");
      }
    }
  }

  commflag = 1;
  //comm->reverse_comm_fix(this,1);

  // check that all procs have needed ghost atoms within ghost cutoff
  // only if neighbor list has changed since last check

  //if (lastcheck < neighbor->lastcall) check_ghosts();

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm

  comm->forward_comm();

  // forward comm of bondcount, so ghosts have it
  //printf ("rank %d chkpnt 0\n", me);
  commflag = 1;
  comm->forward_comm_fix(this,1);
  //printf ("rank %d chkpnt 1\n", me);
  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(to_remove);
    memory->destroy(to_add);
    memory->destroy(recreate);
    memory->destroy(final_to_remove);
    memory->destroy(final_to_add);
    memory->destroy(distsq_c);
    memory->destroy(distsq_b);
    nmax = atom->nmax;
    memory->create(to_remove,nmax,"extrusion:to_remove");
    memory->create(to_add,nmax,"extrusion:to_add");
    memory->create(recreate,nmax,"extrusion:recreate");
    memory->create(final_to_remove,nmax,"extrusion:final_to_remove");
    memory->create(final_to_add,nmax,"extrusion:final_to_add");
    memory->create(distsq_c,nmax,"extrusion:distsq_c");
    probability = distsq_c;
    memory->create(distsq_b,nmax,"extrusion:distsq_b");
  }

  //int nlocal = atom->nlocal;
  //int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    to_remove[i] = 0;
    to_add[i] = 0;
    recreate[i] = 0;
    final_to_remove[i] = 0;
    final_to_add[i] = 0;
    distsq_c[i] = BIG;
    distsq_b[i] = 0.0;
  }

  // loop over bond list
  // setup possible partner list of bonds to break

  double **x = atom->x;//
  tagint *tag = atom->tag;//
  int *mask = atom->mask;//
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int mi;

  //int **bond_type = atom->bond_type;
  //tagint **bond_atom = atom->bond_atom;//
  //int *num_bond = atom->num_bond;//
  int **nspecial = atom->nspecial;//
  tagint **special = atom->special;//

  int swap_buffer;
  bool left_end_moving, right_end_moving;
  bool alreadymet = false;

  //int newton_bond = force->newton_bond;
  //if (!newton_bond) return;

//new lines for extrusion
  for (int numbnd = 0; numbnd < nbondlist; numbnd++) {
    i1 = bondlist[numbnd][0];
    i2 = bondlist[numbnd][1]; //local atom numbers
    type = bondlist[numbnd][2];
    if (!(mask[i1] & groupbit)) continue;
    if (!(mask[i2] & groupbit)) continue;
    if (type != btype) continue;
    if (tag[i1] < tag[i2]) {
      i1 = atom->map(tag[i1]);//local atom numbers double check
      i2 = atom->map(tag[i2]);
    }
    else if (tag[i1] > tag[i2]) {
      mi = i1;
      i1 = atom->map(tag[i2]);//local atom numbers double check
      i2 = atom->map(tag[mi]);
      }
    else {
        error->one(FLERR,"Fix extrusion, bond i-i exists");
     }
 
    left_end_moving = false;
    right_end_moving = false;
    

    //   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //   !!!                                                                             !!!
    //   !!! I should add energetic criteria to shift the bond depending on the distance !!!
    //   !!!                                                                             !!!
    //   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //compare global atom numbers
      //try to move left atom
      //###############################################################################
      //###                                                                        ####
      //### Get the partner id, exchange it, in the new loop over the nall,        ####
      //### then add/remove bonds and store it in the finalpartner.                ####
      //### Only then calculate if number of creation equals to number of breaks.  ####
      //###                                                                        ####
      //###############################################################################
    //printf("i1 and i2: %d %d\n", tag[i1], tag[i2]); 
    if (
      num_bond[i1] == 0 ||
      num_bond[i2] == 0 ||
      num_bond[atom->map(tag[i2]+1)] == 0 ||
      num_bond[atom->map(tag[i1]-1)] == 0 ||
      bondcount[i1] == 0 ||
      bondcount[i2] == 0) continue; // check if bead on another proc
      
    if (
    num_bond[atom->map(tag[i1] - 1)] - bondcount[atom->map(tag[i1] - 1)] == 2 &&
    bondcount[atom->map(tag[i1] - 1)] == 0 &&
    (atom->type[atom->map(tag[i1] - 1)] != ctcf_left || through_prob > random->uniform())) {
      local_left = atom->map(tag[i1] - 1);
      if (
      num_bond[atom->map(tag[i2] + 1)] - bondcount[atom->map(tag[i2] + 1)] == 2 &&
      bondcount[atom->map(tag[i2] + 1)] == 0 &&
      atom->type[atom->map(tag[i2] + 1)] != ctcf_right) {
        local_right = atom->map(tag[i2] + 1);
        delx = x[local_left][0] - x[local_right][0];
        dely = x[local_left][1] - x[local_right][1];
        delz = x[local_left][2] - x[local_right][2];
        rsq = delx*delx + dely*dely + delz*delz;
        //printf("Energy to create: %f\n", 20*rsq);
        if (rsq > distsq_c[local_left] && rsq > distsq_c[local_right]) continue;
        if (rsq < distsq_c[local_left]) {
          distsq_c[local_left] = rsq;
          to_add[local_left] = tag[local_right];
        }
       if (rsq < distsq_c[local_right]) {
          distsq_c[local_right] = rsq;
          to_add[local_right] = tag[local_left];
        }
        //printf("+++ %d\t%d\t%f\t%d\n", to_add[local_left], to_add[local_right], rsq, me);
        delx = x[i1][0] - x[i2][0];
        dely = x[i1][1] - x[i2][1];
        delz = x[i1][2] - x[i2][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq > distsq_b[i1]) { 
          distsq_b[i1] = rsq;
          to_remove[i1]=tag[i2];
        }
        if (rsq > distsq_b[i2]) {
          distsq_b[i2] = rsq;
          to_remove[i2]=tag[i1];
        }
        //printf("--- %d\t%d\t%f\t%d\n", to_remove[i1], to_remove[i2], rsq, me);
      }
      else { // left move, right not
        delx = x[local_left][0] - x[i2][0];
        dely = x[local_left][1] - x[i2][1];
        delz = x[local_left][2] - x[i2][2];
        rsq = delx*delx + dely*dely + delz*delz;
        //printf("Energy to create: %f\n", 20*rsq);
        if (rsq > distsq_c[local_left] && rsq > distsq_c[i2]) continue;
        if (rsq < distsq_c[local_left]) {
          distsq_c[local_left] = rsq;
          to_add[local_left] = tag[i2];
        }
        if (distsq_c[i2] == BIG) {
          distsq_c[i2] = rsq;
          to_add[i2] = tag[local_left];
        }
        else printf ("check me ASAP!!!\n");

        //printf("+++ %d\t%d\t%f\t%d\n", to_add[local_left], to_add[i2], rsq, me);
        delx = x[i1][0] - x[i2][0];
        dely = x[i1][1] - x[i2][1];
        delz = x[i1][2] - x[i2][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq > distsq_b[i1]) { 
          distsq_b[i1] = rsq;
          to_remove[i1]=tag[i2];
        }
        if (rsq > distsq_b[i2]) {
          distsq_b[i2] = rsq;
          to_remove[i2]=tag[i1];
        }
        //printf("--- %d\t%d\t%f\t%d\n", to_remove[i1], to_remove[i2], rsq, me);
      }
    }
    else if (
    num_bond[atom->map(tag[i2] + 1)] - bondcount[atom->map(tag[i2] + 1)] == 2 &&
    bondcount[atom->map(tag[i2] + 1)] == 0 &&
    (atom->type[atom->map(tag[i2] + 1)] != ctcf_right || through_prob > random->uniform())) {
      local_right = atom->map(tag[i2] + 1);
      delx = x[i1][0] - x[local_right][0];
      dely = x[i1][1] - x[local_right][1];
      delz = x[i1][2] - x[local_right][2];
      rsq = delx*delx + dely*dely + delz*delz;
      //printf("Energy to create: %f\n", 20*rsq);
      if (rsq > distsq_c[i1] && rsq > distsq_c[local_right]) continue;
      if (distsq_c[i1] == BIG) {
        distsq_c[i1] = rsq;
        to_add[i1] = tag[local_right];
      }
      else printf("right move left not, problem!\n");
      if (rsq < distsq_c[local_right]) {
        distsq_c[local_right] = rsq;
        to_add[local_right] = tag[i1];
      }
      //printf("+++ %d\t%d\t%f\t%d\n", to_add[i1], to_add[local_right], rsq, me);
      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq > distsq_b[i1]) { 
        distsq_b[i1] = rsq;
        to_remove[i1]=tag[i2];
      }
      if (rsq > distsq_b[i2]) {
        distsq_b[i2] = rsq;
        to_remove[i2]=tag[i1];
      }
      //printf("--- %d\t%d\t%f\t%d\n", to_remove[i1], to_remove[i2], rsq, me);
    }
  }
/*
  // reverse comm of partner info
  commflag=2;
  comm->reverse_comm_fix(this);
  
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("2.5---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("2.5++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }

  for (i = 0; i < nall; i++) {
    if (recreate[i] != 0) 
    {
      recreate[atom->map(recreate[i])] = tag[i];
      printf ("recreate %d isn't zero at %d\n", tag[i], me);
    }
  }
*/
//  commflag = 23;
//  comm->reverse_comm_fix(this,1);

//  for (i = 0; i < nall; i++) {
//    if (recreate[i] != 0) recreate[atom->map(recreate[i])] = tag[i];
//  }
/*
// nullify both atoms because you have already checked them 
  for (i = 0; i < nlocal; i++) {
    if (recreate[i] == 0) continue;
    if (recreate[i] < tag[i]) {
      lb = recreate[i];
      rb = tag[i];
      l_flag=false;
    }
    else {
      lb = tag[i];
      rb = recreate[i];
      l_flag=true;
    }
    //if (l_flag) {
      if (to_remove[atom->map(lb+1)] - rb == 0 || to_remove[atom->map(lb+1)] - rb == -1) {
        to_remove[atom->map(lb+1)] = 0;
        distsq_b[atom->map(lb+1)] = BIG;
        to_remove[atom->map(rb)] = 0;
        distsq_b[atom->map(rb)] = BIG;
      } 
      else if (to_remove[atom->map(lb)] - rb == -1) {
        to_remove[atom->map(lb)] = 0;
        distsq_b[atom->map(lb)] = BIG;
        to_remove[atom->map(rb)] = 0;
        distsq_b[atom->map(rb)] = BIG;
      }
    }
    else {
      if (lb - to_remove[atom->map(rb-1)] == 0 || lb - to_remove[atom->map(rb-1)] == -1) {
        to_remove[atom->map(rb-1)] = 0;
        distsq_b[atom->map(rb-1)] = BIG;
        to_remove[atom->map(lb)] = 0;
        distsq_b[atom->map(lb)] = BIG;
      }
      else if (lb - to_remove[atom->map(rb)] == -1) {
        to_remove[atom->map(rb)] = 0;
        distsq_b[atom->map(rb)] = BIG;
        to_remove[atom->map(lb)] = 0;
        distsq_b[atom->map(lb)] = BIG;
      }
    }
  //}
  commflag=3;
  comm->reverse_comm_fix(this);
  commflag = 92;
  comm->forward_comm_fix(this,1);
  commflag = 91;
  comm->forward_comm_fix(this,1);
  /*
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("2.7---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("2.7++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }
  
//final check starting from to_add
  for (i = 0; i < nlocal; i++) {
    if (to_add[i] != 0) {
      if (tag[i] == to_add[atom->map(to_add[i])]) {
        if (to_add[i] < tag[i]) {
          lb = to_add[i];
          rb = tag[i];
          l_flag=false;
        }
        else {
          lb = tag[i];
          rb = to_add[i];
          l_flag=true;
        }
        if (
          (to_remove[atom->map(lb + 1)] == rb &&      to_remove[atom->map(rb)] == lb + 1) || 
          (to_remove[atom->map(lb + 1)] == rb - 1 &&  to_remove[atom->map(rb - 1)] == lb + 1) || 
          (to_remove[atom->map(lb)] == rb - 1 &&      to_remove[atom->map(rb - 1)] == lb)) {
          continue;
        }
        else {
          to_add[atom->map(to_add[i])] = 0;
          distsq_c[atom->map(to_add[i])] = 0;
          to_add[i] = 0;
          distsq_c[i] = 0;
        }
      }
      else {
        to_add[i] = 0;
        distsq_c[i] = 0;
      }
    }
  }
  /*
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("2.8---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("2.8++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }

//final check starting from to_remove
  for (i = 0; i < nlocal; i++) {
    if (to_remove[i] != 0) {
      if (tag[i] == to_remove[atom->map(to_remove[i])]) {
        if (to_remove[i] < tag[i]) {
          lb = to_remove[i];
          rb = tag[i];
          l_flag=false;
        }
        else {
          lb = tag[i];
          rb = to_remove[i];
          l_flag=true;
        }
        if (
          (to_add[atom->map(lb - 1)] == rb &&     to_add[atom->map(rb)] == lb - 1) ||
          (to_add[atom->map(lb - 1)] == rb + 1 && to_add[atom->map(rb + 1)] == lb - 1) ||
          (to_add[atom->map(lb)] == rb + 1 &&     to_add[atom->map(rb + 1)] == lb)) {
          continue; // pair to_remove has complementary pair to_add
        }
        else {
          to_remove[atom->map(to_remove[i])] = 0;
          distsq_b[atom->map(to_remove[i])] = BIG;
          to_remove[i] = 0;
          distsq_b[i] = BIG;
        }        
      }
      else {
        to_remove[i] = 0;
        distsq_b[i] = BIG;
      }
    }
  }
  commflag=3;
  comm->reverse_comm_fix(this);
  commflag=2;
  comm->reverse_comm_fix(this);
  /*
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("2.9---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("2.9++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }

/*
  for (i = 0; i < nlocal; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i]) {
      //printf("2---- recheck to_remove %d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])]);
      if (to_remove[i] < tag[i]) {
        lb = to_remove[i];
        rb = tag[i];
        l_flag = false;
      }
      else {
        lb = tag[i];
        rb = to_remove[i];
        l_flag = true;
      }
      if (l_flag) {
        if (rb - to_add[atom->map(lb-1)] == 0 || rb - to_add[atom->map(lb-1)] == -1) {
          to_add[atom->map(lb-1)] = 0;
          distsq_c[atom->map(lb-1)] = BIG;
        }
        else if (rb - to_add[atom->map(lb)] == -1) {
          to_add[atom->map(lb)] = 0;
          distsq_c[atom->map(lb)] = BIG;
        }
      }
      else {
        if (to_add[atom->map(rb+1)] - lb == 0 || to_add[atom->map(rb+1)] - lb == -1) {
          to_add[atom->map(rb+1)] = 0;
          distsq_c[atom->map(rb+1)] = BIG;
        }
        else if (to_add[atom->map(rb)] - lb == -1) {
          to_add[atom->map(rb)] = 0;
          distsq_c[atom->map(rb)] = BIG;
        }
      }
        (rb - to_add[atom->map(lb-1)] == 0 || rb - to_add[atom->map(lb-1)] == -1) && 
        (to_add[atom->map(rb+1)] - lb == 0 || to_add[atom->map(rb+1)] - lb == -1)) {
        if (to_add[atom->map(lb-1)] == rb && to_add[atom->map(rb+1)] == lb) {
          if (l_flag) {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(lb-1)])] = 0;
            distsq_c[atom->map(to_add[atom->map(lb-1)])] = BIG;
            to_add[atom->map(lb-1)] = 0;
            distsq_c[atom->map(lb-1)] = BIG;
          }
          else {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(rb+1)])] = 0;
            distsq_c[atom->map(to_add[atom->map(rb+1)])] = BIG;
            to_add[atom->map(rb+1)] = 0;
            distsq_c[atom->map(rb+1)] = BIG;
          }
          continue;
        }
        //printf("2.5---+ recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_add[atom->map(lb-1)], to_add[atom->map(rb+1)]);
        continue;
      }
      else if (
        (rb - to_add[atom->map(lb-1)] == 0 || rb - to_add[atom->map(lb-1)] == -1) && 
        (to_add[atom->map(rb)] - lb == 0 || to_add[atom->map(rb)] - lb == -1)) {
        if (to_add[atom->map(lb-1)] == rb && to_add[atom->map(rb+1)] == lb) {
          if (l_flag) {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(lb-1)])] = 0;
            distsq_c[atom->map(to_add[atom->map(lb-1)])] = BIG;
            to_add[atom->map(lb-1)] = 0;
            distsq_c[atom->map(lb-1)] = BIG;
          }
          else {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(rb)])] = 0;
            distsq_c[atom->map(to_add[atom->map(rb)])] = BIG;
            to_add[atom->map(rb)] = 0;
            distsq_c[atom->map(rb)] = BIG;
          }
          continue;
        }
        //printf("2.5---+ recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_add[atom->map(lb-1)], to_add[atom->map(rb)]);
        continue;
      }
      else if (
        (rb - to_add[atom->map(lb)] == 0 || rb - to_add[atom->map(lb)] == -1) && 
        (to_add[atom->map(rb+1)] - lb == 0 || to_add[atom->map(rb+1)] - lb == -1)) {
          if (to_add[atom->map(lb-1)] == rb && to_add[atom->map(rb+1)] == lb) {
          if (l_flag) {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(lb)])] = 0;
            distsq_c[atom->map(to_add[atom->map(lb)])] = BIG;
            to_add[atom->map(lb)] = 0;
            distsq_c[atom->map(lb)] = BIG;
          }
          else {
            to_remove[atom->map(to_remove[i])] = 0;
            distsq_b[atom->map(to_remove[i])] = 0;
            to_remove[i] = 0;
            distsq_b[i] = 0;
            to_add[atom->map(to_add[atom->map(rb+1)])] = 0;
            distsq_c[atom->map(to_add[atom->map(rb+1)])] = BIG;
            to_add[atom->map(rb+1)] = 0;
            distsq_c[atom->map(rb+1)] = BIG;
          }
          continue;
        }
        //printf("2.5---+ recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_add[atom->map(lb)], to_add[atom->map(rb+1)]);
        continue;
      }
      else {
        to_remove[i] = 0;
        distsq_b[i] = 0;
      }
    }
  }
  
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i]) {
      //printf("2++++ recheck to_add %d\t%d\n", to_add[i], to_add[atom->map(to_add[i])]);
      if (to_add[i] < tag[i]) {
        lb = to_add[i];
        rb = tag[i];
        l_flag = true;
      }
      else {
        lb = tag[i];
        rb = to_add[i];
        l_flag = false;
      }
      if ((to_remove[atom->map(lb+1)] - rb == 0 || to_remove[atom->map(lb+1)] - rb == -1) && 
      (lb - to_remove[atom->map(rb-1)] == 0 || lb - to_remove[atom->map(rb-1)] == -1)) {
        if (to_remove[atom->map(lb+1)] == rb && to_remove[atom->map(rb-1)] == lb) {
          if (l_flag) {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(lb+1)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(lb+1)])] = 0;
            to_remove[atom->map(lb+1)] = 0;
            distsq_b[atom->map(lb+1)] = 0;
          }
          else {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(rb-1)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(rb-1)])] = 0;
            to_remove[atom->map(rb-1)] = 0;
            distsq_b[atom->map(rb-1)] = 0;
          }
          continue;
        }
        //printf ("2.5+++- recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_remove[atom->map(lb+1)], to_remove[atom->map(rb-1)]);
        continue;
      }
      else if ((to_remove[atom->map(lb+1)] - rb == 0 || to_remove[atom->map(lb+1)] - rb == -1) && 
      (lb - to_remove[atom->map(rb)] == 0 || lb - to_remove[atom->map(rb)] == -1)) {
        if (to_remove[atom->map(lb+1)] == rb && to_remove[atom->map(rb)] == lb) {
          if (l_flag) {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(lb+1)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(lb+1)])] = 0;
            to_remove[atom->map(lb+1)] = 0;
            distsq_b[atom->map(lb+1)] = 0;
          }
          else {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(rb)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(rb)])] = 0;
            to_remove[atom->map(rb)] = 0;
            distsq_b[atom->map(rb)] = 0;
          }
          continue;
        }
        //printf ("2.5+++- recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_remove[atom->map(lb+1)], to_remove[atom->map(rb)]);
        continue;
      }
      else if ((to_remove[atom->map(lb)] - rb == 0 || to_remove[atom->map(lb)] - rb == -1) && 
      (lb - to_remove[atom->map(rb-1)] == 0 || lb - to_remove[atom->map(rb-1)] == -1)) {
        if (to_remove[atom->map(lb)] == rb && to_remove[atom->map(rb-1)] == lb) {
          if (l_flag) {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(lb)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(lb)])] = 0;
            to_remove[atom->map(lb)] = 0;
            distsq_b[atom->map(lb)] = 0;
          }
          else {
            to_add[atom->map(to_add[i])] = 0;
            distsq_c[atom->map(to_add[i])] = BIG;
            to_add[i] = 0;
            distsq_c[i] = BIG;
            to_remove[atom->map(to_remove[atom->map(rb-1)])] = 0;
            distsq_b[atom->map(to_remove[atom->map(rb-1)])] = 0;
            to_remove[atom->map(rb-1)] = 0;
            distsq_b[atom->map(rb-1)] = 0;
          }
          continue;
        }
          //printf ("2.5+++- recheck to_add %d\t%d\t%d\t%d\n", lb, rb, to_remove[atom->map(lb)], to_remove[atom->map(rb-1)]);
          continue;
        }
      else {
        to_add[i] = 0;
        distsq_c[i] = BIG;
      }
    }
  }
  
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i]) 
      printf("3---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i]) 
      printf("3++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }
  
*/
  commflag=2;
  comm->reverse_comm_fix(this);
  commflag=3;
  comm->reverse_comm_fix(this);
/*
  for (i = 0; i < nlocal; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i]) 
      printf("3.5---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nlocal; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i]) 
      printf("3.5++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }
*/
  // each atom now knows its winning partner
  commflag = 91;
  comm->forward_comm_fix(this,1);
  /*
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("3.7---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("3.7++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }
  */
  commflag = 92;
  comm->forward_comm_fix(this,1);
  /*
  for (i = 0; i < nall; i++) {
    if (to_remove[i] != 0 && to_remove[atom->map(to_remove[i])] == tag[i] && me==16) 
      printf("4---- recheck to_remove %d\t%d\t%d\n", to_remove[i], to_remove[atom->map(to_remove[i])], me);
  }
  for (i = 0; i < nall; i++) {
    if (to_add[i] != 0 && to_add[atom->map(to_add[i])] == tag[i] && me==16) 
      printf("4++++ recheck to_add %d\t%d\t%d\n", to_add[i], to_add[atom->map(to_add[i])], me);
  }
  */
  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  Now loop over bonds to be removed.  .  .  .  .  .  .  .  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
  
  nbreak=0;
  for (i = 0; i < nlocal; i++) {
    if (to_remove[i] == 0) continue;
    j = atom->map(to_remove[i]);
    j = atom->map(tag[j]);
    if (to_remove[j] != tag[i]) continue;
    if (to_remove[i] < tag[i]) {
      lb = to_remove[i];
      rb = tag[i];
    }
    else {
      lb = tag[i];
      rb = to_remove[i];
    }
    if (
      to_add[atom->map(lb - 1)] == rb &&     to_add[atom->map(rb)] == lb - 1 &&
      to_add[atom->map(lb)] == rb + 1 &&     to_add[atom->map(rb + 1)] == lb) {
        to_add[atom->map(lb - 1)] = rb + 1;
        to_add[atom->map(rb + 1)] = lb - 1;
        to_add[atom->map(lb)] = 0;
        to_add[atom->map(rb)] = 0;
      }
    
    if (
      (to_add[atom->map(lb - 1)] == rb &&     to_add[atom->map(rb)] == lb - 1) ||// &&    to_add[atom->map(rb + 1)] != lb)  ||
      (to_add[atom->map(lb - 1)] == rb + 1 && to_add[atom->map(rb + 1)] == lb - 1) ||
      (to_add[atom->map(lb)] == rb + 1 &&     to_add[atom->map(rb + 1)] == lb)) {// &&    to_add[atom->map(rb)] != lb - 1)) {

    // delete bond from atom I if I stores it
    // atom J will also do this

      for (m = 0; m < num_bond[i]; m++) {
        if (bond_atom[i][m] == to_remove[i]) {
          for (k = m; k < num_bond[i]-1; k++) {
            bond_atom[i][k] = bond_atom[i][k+1];
            bond_type[i][k] = bond_type[i][k+1];
          }
          num_bond[i]--;
          break;
        }
      }

      // remove J from special bond list for atom I
      // atom J will also do this, whatever proc it is on

      slist = special[i];
      n1 = nspecial[i][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == to_remove[i]) break;
      n3 = nspecial[i][2];
      for (; m < n3-1; m++) slist[m] = slist[m+1];
      nspecial[i][0]--;
      nspecial[i][1]--;
      nspecial[i][2]--;

      // store final broken bond partners and count the broken bond once
      //printf ("third assignment remove: %d %d\n", tag[i], tag[j]);
      final_to_remove[i] = tag[j];
      final_to_remove[j] = tag[i];
      //printf("final to remove: %d, %d\n", tag[i], tag[j]);
      if (tag[i] < tag[j]) nbreak++;
      }
  }

  //.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  Now loops over bonds will be formed .  .  .  .  .  .  . 
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
//.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
  ncreate=0;
  for (i = 0; i < nlocal; i++) {
    if (to_add[i] == 0) continue;
    j = atom->map(to_add[i]);
    j = atom->map(tag[j]);
    if (to_add[j] != tag[i]) continue;
    if (to_add[i] < tag[i]) {
      lb = to_add[i];
      rb = tag[i];
    }
    else {
      lb = tag[i];
      rb = to_add[i];
    }
    //if (me == 8 && lb == 168 && rb == 172) printf ("to_add %d\t%d\t to_remove %d\t%d\t%d\t%d\n", 
    //lb, rb, to_remove[atom->map(lb + 1)], to_remove[atom->map(lb)], to_remove[atom->map(rb)], to_remove[atom->map(rb - 1)]);
    if ((
      to_remove[atom->map(lb + 1)] == rb &&
      to_remove[atom->map(rb)] == lb + 1) &&
      to_remove[atom->map(lb)] == rb - 1 &&
      to_remove[atom->map(rb - 1)] == lb) {

      }
    if (
      (to_remove[atom->map(lb + 1)] == rb &&      to_remove[atom->map(rb)] == lb + 1) || 
      (to_remove[atom->map(lb + 1)] == rb - 1 &&  to_remove[atom->map(rb - 1)] == lb + 1) || 
      (to_remove[atom->map(lb)] == rb - 1 &&      to_remove[atom->map(rb - 1)] == lb)) {
    //printf("to add: %d, %d\n", tag[i], tag[j]);

    // if newton_bond is set, only store with I or J
    // if not newton_bond, store bond with both I and J
    // atom J will also do this consistently, whatever proc it is on
    //if (!force->newton_bond || tag[i] < tag[j]) {
      if (num_bond[i] == atom->bond_per_atom){
        printf("problem with bonds %d, %d, %d, proc %d\n", num_bond[i], tag[i], atom->bond_per_atom, me);
        for (int myit = 0; myit < num_bond[i]; myit++)
        {
          printf("bond %d\n", bond_atom[i][myit]);
        }
        printf("bond to add %d\n", tag[j]);
        printf("num bonds type 2 %d\n", bondcount[i]);
        error->one(FLERR,"New bond exceeded bonds per atom in fix extrusion");
      }
      bond_type[i][num_bond[i]] = btype;
      bond_atom[i][num_bond[i]] = tag[j];
      num_bond[i]++;
    //}

    // add a 1-2 neighbor to special bond list for atom I
    // atom J will also do this, whatever proc it is on
    // need to first remove tag[j] from later in list if it appears
    // prevents list from overflowing, will be rebuilt in rebuild_special_one()

      slist = special[i];
      n1 = nspecial[i][0];
      n2 = nspecial[i][1];
      n3 = nspecial[i][2];
      for (m = n1; m < n3; m++)
        if (slist[m] == tag[j]) break;
      if (m < n3) {
        for (n = m; n < n3-1; n++) slist[n] = slist[n+1];
        n3--;
        if (m < n2) n2--;
      }
      if (n3 == atom->maxspecial)
        error->one(FLERR,
                  "New bond exceeded special list size in fix extrusion");
      for (m = n3; m > n1; m--) slist[m] = slist[m-1];
      slist[n1] = tag[j];
      nspecial[i][0] = n1+1;
      nspecial[i][1] = n2+1;
      nspecial[i][2] = n3+1;

      // increment bondcount, convert atom to new type if limit reached
      // atom J will also do this, whatever proc it is on

      bondcount[i]++;

      // store final created bond partners and count the created bond once

      final_to_add[i] = tag[j];
      final_to_add[j] = tag[i];
      if (tag[i] < tag[j]) ncreate++;
      }
    //printf("final to add: %d, %d\n", tag[i], tag[j]);
  }


  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  // trigger reneighboring if any shifts happened
  // this insures neigh lists will immediately reflect the topology changes
  // done if no shifts

  if (!breakcount && !createcount) return;
  else if (breakcount == createcount) next_reneighbor = update->ntimestep;
  else if (breakcount != createcount) {
    printf(" createcount %d, breakcount %d, timestep %d\n", createcount, breakcount, update->ntimestep);
    for (int b=0; b<nlocal; b++)
    {
      //if (abs(tag[bondlist[b][0]] - tag[bondlist[b][1]])>1) printf("_$ %d\t%d\t%d\n", tag[bondlist[b][0]], tag[bondlist[b][1]], me);
      if (final_to_remove[b] != 0) printf("-$ %d\t%d\t%d\t%d\n", tag[b], final_to_remove[atom->map(final_to_remove[b])], final_to_remove[b], me);
      if (final_to_add[b] != 0) printf("+$ %d\t%d\t%d\t%d\n", tag[b], final_to_add[atom->map(final_to_add[b])], final_to_add[b], me);
    }
    error->all(FLERR,"Numbers of created and broken bonds are not equal");
  }
  
  //if (breakcount) next_reneighbor = update->ntimestep;
  //if (!breakcount) return;

  //breakcounttotal += breakcount;

  // communicate final partner and 1-2 special neighbors
  // 1-2 neighs already reflect broken bonds

  //printf ("rank %d chkpnt 7\n", me);
  commflag = 3;
  comm->forward_comm_fix(this);
//printf ("rank %d chkpnt 8\n", me);
  // create list of broken bonds that influence my owned atoms
  //   even if between owned-ghost or ghost-ghost atoms
  // finalpartner is now set for owned and ghost atoms so loop over nall
  // OK if duplicates in broken list due to ghosts duplicating owned atoms
  // check J < 0 to insure a broken bond to unknown atom is included
  //   i.e. bond partner outside of cutoff length

  
  nbreak = 0;
  for (i = 0; i < nall; i++) {
    if (final_to_remove[i] == 0) continue;
      j = atom->map(final_to_remove[i]);
      if (j < 0 || tag[i] < tag[j]) {
        if (nbreak == maxbreak) {
          maxbreak += DELTA;
          memory->grow(broken,maxbreak,2,"extrusion:broken");
        }
        broken[nbreak][0] = tag[i];
        broken[nbreak][1] = final_to_remove[i];
        nbreak++;
      }
  }


//update to_add
//printf ("rank %d chkpnt 9\n", me);
  commflag = 2;
  comm->forward_comm_fix(this);
//printf ("rank %d chkpnt 10\n", me);

  ncreate = 0;
  for (i = 0; i < nall; i++) {
    if (final_to_add[i] == 0) continue;
      j = atom->map(final_to_add[i]);
      if (j < 0 || tag[i] < tag[j]) {
        if (ncreate == maxcreate) {
          maxcreate += DELTA;
          memory->grow(created,maxcreate,2,"extrusion:created");
        }
        created[ncreate][0] = tag[i];
        created[ncreate][1] = final_to_add[i];
        ncreate++;
      }
  }

  // update special neigh lists of all atoms affected by any broken bond
  // also remove angles/dihedrals/impropers broken by broken bonds

  update_topology();

  // DEBUG
  //print_bb();
}

/* ----------------------------------------------------------------------
   insure all atoms 2 hops away from owned atoms are in ghost list
   this allows dihedral 1-2-3-4 to be properly deleted
     and special list of 1 to be properly updated
   if I own atom 1, but not 2,3,4, and bond 3-4 is deleted
     then 2,3 will be ghosts and 3 will store 4 as its finalpartner
------------------------------------------------------------------------- */

void FixExtrusion::check_ghosts()
{
  int i,j,n;
  tagint *slist;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (i = 0; i < nlocal; i++) {
    slist = special[i];
    n = nspecial[i][1];
    for (j = 0; j < n; j++) {
      if (atom->map(slist[j]) < 0) flag = 1;
    }
    //printf ("\n");
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) {
    printf("ATTENTION!   %d\n", flagall);
    error->all(FLERR,"Fix extrusion needs ghost atoms from further away");    
  }
  lastcheck = update->ntimestep;
}

/* ----------------------------------------------------------------------
   double loop over my atoms and broken bonds
   influenced = 1 if atom's topology is affected by any broken bond
     yes if is one of 2 atoms in bond
     yes if both atom IDs appear in atom's special list
     else no
   if influenced:
     check for angles/dihedrals/impropers to break due to specific broken bonds
     rebuild the atom's special list of 1-2,1-3,1-4 neighs
------------------------------------------------------------------------- */

void FixExtrusion::update_topology()
{
  int i,j,k,n,influence,influenced,found;
  tagint id1,id2;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  //printf("NBREAK %d: ",nbreak);
  //for (i = 0; i < nbreak; i++)
  //  printf(" %d %d,",broken[i][0],broken[i][1]);
  //printf("\n");

  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    for (j = 0; j < nbreak; j++) {
      id1 = broken[j][0];
      id2 = broken[j][1];

      influence = 0;
      if (tag[i] == id1 || tag[i] == id2) influence = 1;
      else {
        n = nspecial[i][2];
        found = 0;
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) found++;
        if (found == 2) influence = 1;
      }
      if (!influence) continue;
      influenced = 1;
    }

    if (influenced) rebuild_special_one(i);
  }

  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    for (j = 0; j < ncreate; j++) {
      id1 = created[j][0];
      id2 = created[j][1];

      influence = 0;
      if (tag[i] == id1 || tag[i] == id2) influence = 1;
      else {
        n = nspecial[i][1];
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) {
            influence = 1;
            break;
          }
      }
      if (!influence) continue;
      influenced = 1;
    }
    if (influenced) rebuild_special_one(i);
  }
}

  /*


  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    for (j = 0; j < nbreak; j++) {
      id1 = broken[j][0];
      id2 = broken[j][1];

      influence = 0;
      if (tag[i] == id1 || tag[i] == id2) influence = 1;
      else {
        n = nspecial[i][2]; //interesting what is right: [1] bond_create, [2] bond_break
        found = 0;
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) found++;
        if (found == 2) influence = 1;
        n = nspecial[i][1];
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) {
            influence = 1;
            break;
          }
      }
      if (!influence) continue;
      influenced = 1;
    }

    if (influenced) rebuild_special_one(i);
  }
}
*/

/* ----------------------------------------------------------------------
   re-build special list of atom M
   does not affect 1-2 neighs (already include effects of new bond)
   affects 1-3 and 1-4 neighs due to other atom's augmented 1-2 neighs
------------------------------------------------------------------------- */

void FixExtrusion::rebuild_special_one(int m)
{
  int i,j,n,n1,cn1,cn2,cn3;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // existing 1-2 neighs of atom M

  slist = special[m];
  n1 = nspecial[m][0];
  cn1 = 0;
  for (i = 0; i < n1; i++)
    copy[cn1++] = slist[i];

  // new 1-3 neighs of atom M, based on 1-2 neighs of 1-2 neighs
  // exclude self
  // remove duplicates after adding all possible 1-3 neighs

  cn2 = cn1;
  for (i = 0; i < cn1; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn2++] = slist[j];
  }

  cn2 = dedup(cn1,cn2,copy);
  if (cn2 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix bond/create");
  // new 1-4 neighs of atom M, based on 1-2 neighs of 1-3 neighs
  // exclude self
  // remove duplicates after adding all possible 1-4 neighs

  cn3 = cn2;
  for (i = cn1; i < cn2; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn3++] = slist[j];
  }

  cn3 = dedup(cn2,cn3,copy);
  if (cn3 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix bond/create");
  // store new special list with atom M

  nspecial[m][0] = cn1;
  nspecial[m][1] = cn2;
  nspecial[m][2] = cn3;
  memcpy(special[m],copy,cn3*sizeof(int));
}

/* ----------------------------------------------------------------------
   remove all ID duplicates in copy from Nstart:Nstop-1
   compare to all previous values in copy
   return N decremented by any discarded duplicates
------------------------------------------------------------------------- */

int FixExtrusion::dedup(int nstart, int nstop, tagint *copy)
{
  int i;

  int m = nstart;
  while (m < nstop) {
    for (i = 0; i < m; i++)
      if (copy[i] == copy[m]) {
        copy[m] = copy[nstop-1];
        nstop--;
        break;
      }
    if (i == m) m++;
  }

  return nstop;
}

/* ---------------------------------------------------------------------- */

void FixExtrusion::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixExtrusion::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m,ns;
  m = 0;

  if (commflag == 1) { //fbc cf=1
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(bondcount[j]).d;
    }
    return m;
  }

  if (commflag == 2) { //fbc cf=3
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(final_to_add[j]).d;
      ns = nspecial[j][0];
      buf[m++] = ubuf(ns).d;
      for (k = 0; k < ns; k++)
        buf[m++] = ubuf(special[j][k]).d;
    }
    return m;
  }
  if (commflag == 3) { //fbb cf=2
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(final_to_remove[j]).d;
      ns = nspecial[j][0];
      buf[m++] = ubuf(ns).d;
      for (k = 0; k < ns; k++)
        buf[m++] = ubuf(special[j][k]).d;
    }
    return m;
  }
  
  if (commflag == 91) { //fbb cf=1
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(to_remove[j]).d;
    }
    return m;
  }
  
  if (commflag == 92) { //fbc cf=2
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(to_add[j]).d;
    }
    return m;
  }

  error->all(FLERR,"Wrong commflag in FixExtrusion::pack_forward_comm.");
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixExtrusion::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,m,ns,last;
  m=0;
  last = first + n;
  
  if (commflag == 1) { //fbc cf=1
    for (i = first; i < last; i++)
      bondcount[i] = (int) ubuf(buf[m++]).i;
  } 
  
  if (commflag == 2) { //fbc cf=3
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (i = first; i < last; i++) {
      final_to_add[i] = (tagint) ubuf(buf[m++]).i;
      ns = (int) ubuf(buf[m++]).i;
      nspecial[i][0] = ns;
      for (j = 0; j < ns; j++)
        special[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
  
  if (commflag == 3) { //fbb cf=2
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (i = first; i < last; i++) {
      final_to_remove[i] = (tagint) ubuf(buf[m++]).i;
      ns = (int) ubuf(buf[m++]).i;
      nspecial[i][0] = ns;
      for (j = 0; j < ns; j++)
        special[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
  
  if (commflag == 91) { //fbb cf=1
    for (i = first; i < last; i++) {
      to_remove[i] = (tagint) ubuf(buf[m++]).i;
    }
  }
  
  if (commflag == 92) { //fbc cf=2
    for (i = first; i < last; i++) {
      to_add[i] =    (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixExtrusion::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 1) { //rbc cf=1
    for (i = first; i < last; i++) {
      buf[m++] = ubuf(bondcount[i]).d;
    }
    return m;
  }
  
  if (commflag == 2) { //rbc cf=2
    for (i = first; i < last; i++) {
      buf[m++] = ubuf(to_add[i]).d;
      buf[m++] = distsq_c[i];
    }
  return m;
  }
  
  if (commflag == 3) { //rbb cf=1
    for (i = first; i < last; i++) {
      //printf("to_remove sends %d\t%d\n", ubuf(to_remove[i]).d,to_remove[i]);
      buf[m++] = ubuf(to_remove[i]).d;
      buf[m++] = distsq_b[i];
    }
  return m;
  }

  if (commflag == 25) { //rbc cf=2
    for (i = first; i < last; i++) {
      buf[m++] = ubuf(to_add[i]).d;
      buf[m++] = distsq_c[i];
    }
  return m;
  }
  
  if (commflag == 23) {
    for (i = first; i < last; i++) {
      //printf("to_remove sends %d\t%d\n", ubuf(to_remove[i]).d,to_remove[i]);
      buf[m++] = ubuf(recreate[i]).d;
    }
  return m;
  }
  error->all(FLERR,"Wrong commflag in FixExtrusion::pack_reverse_comm.");
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixExtrusion::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 1) { //rbc cf=1
    for (i = 0; i < n; i++) {
      j = list[i];
      bondcount[j] += (int) ubuf(buf[m++]).i;
    }
  }
  
  if (commflag == 2) { //rbc cf=2
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m+1] < distsq_c[j]) {
        to_add[j] = (tagint) ubuf(buf[m++]).i;
        distsq_c[j] = buf[m++];
        //printf("unpack reverse to add %d, %d, %d, %d\n", to_add[j], distsq_c[j], i, j);
      } 
      else m+=2;
    }
  }
  else if (commflag == 3) { //rbb cf=1
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m+1] > distsq_b[j]) {
        //printf("to_remove recieves %d\n",(tagint) ubuf(buf[m+1]).i);
        to_remove[j] = (tagint) ubuf(buf[m++]).i;
        distsq_b[j] = buf[m++];
      } else m+=2; 
    }
  }
  else if (commflag == 25) { //rbb cf=1
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m+1] < distsq_c[j]) {
        //printf("to_remove recieves %d\n",(tagint) ubuf(buf[m+1]).i);
        to_add[j] = (tagint) ubuf(buf[m++]).i;
        distsq_c[j] = buf[m++];
      } else m+=2; 
    }
  }
  else if (commflag == 23) {
    for (i = 0; i < n; i++) {
      j = list[i];
      //if (buf[m] > 0 && recreate[j] > 0) error->all(FLERR,"Overwrite recreate in FixExtrusion::pack_reverse_comm.");
      if (buf[m] > 0) {
        //printf("to_remove recieves %d\n",(tagint) ubuf(buf[m+1]).i);
        recreate[j] = (tagint) ubuf(buf[m++]).i;
      } 
      else m++;
    }
  }
}


/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixExtrusion::grow_arrays(int nmax)
{
  memory->grow(bondcount,nmax,"bond/create:bondcount");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixExtrusion::copy_arrays(int i, int j, int /*delflag*/)
{
  bondcount[j] = bondcount[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixExtrusion::pack_exchange(int i, double *buf)
{
  buf[0] = bondcount[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixExtrusion::unpack_exchange(int nlocal, double *buf)
{
  bondcount[nlocal] = static_cast<int> (buf[0]);
  return 1;
}

/* ---------------------------------------------------------------------- */


void FixExtrusion::print_bb()
{
  for (int i = 0; i < atom->nlocal; i++) {
    printf("TAG " TAGINT_FORMAT ": %d nbonds: ",atom->tag[i],atom->num_bond[i]);
    for (int j = 0; j < atom->num_bond[i]; j++) {
      printf(" " TAGINT_FORMAT, atom->bond_atom[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d %d %d nspecial: ",atom->tag[i],
           atom->nspecial[i][0],atom->nspecial[i][1],atom->nspecial[i][2]);
    for (int j = 0; j < atom->nspecial[i][2]; j++) {
      printf(" " TAGINT_FORMAT, atom->special[i][j]);
    }
    printf("\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixExtrusion::print_copy(const char *str, tagint m,
                              int n1, int n2, int n3, int *v)
{
  printf("%s " TAGINT_FORMAT ": %d %d %d nspecial: ",str,m,n1,n2,n3);
  for (int j = 0; j < n3; j++) printf(" %d",v[j]);
  printf("\n");
}

/* ---------------------------------------------------------------------- */

double FixExtrusion::compute_vector(int n)
{
  if (n == 0) return (double) breakcount;
  return (double) breakcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixExtrusion::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 2*nmax * sizeof(tagint);
  bytes += nmax * sizeof(double);
  return bytes;
}
