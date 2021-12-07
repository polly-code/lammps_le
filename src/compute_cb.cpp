/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_cb.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCB::ComputeCB(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 6)
    error->all(FLERR, "Illegal compute bond command");

  scalar_flag = 1;
  extscalar = 1;

  bt = utils::inumeric(FLERR, arg[3], false, lmp);
  a1 = utils::inumeric(FLERR, arg[4], false, lmp);
  a2 = utils::inumeric(FLERR, arg[5], false, lmp);
}

/* ---------------------------------------------------------------------- */

double ComputeCB::compute_scalar()
{
  emine = count_bonds(bt, a1, a2);
  MPI_Allreduce(&emine, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  //MPI_Allreduce(emine, vector, nsub, MPI_DOUBLE, MPI_SUM, world);
  return scalar;
}

double ComputeCB::count_bonds(int bt, int a1, int a2)
{
  int i, atom1, atom2;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int m = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++)
  {
    if (!(mask[atom1] & groupbit))
      continue;
    for (i = 0; i < num_bond[atom1]; i++)
    {
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & groupbit))
        continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2])
        continue;
      if (bond_type[atom1][i] == 0)
        continue;
      if (bond_type[atom1][i] == bt && tag[atom1] == a1 && tag[atom2] == a2)
      {
        return 1.0d;
      }
    }
  }
  return 0.0d;
}