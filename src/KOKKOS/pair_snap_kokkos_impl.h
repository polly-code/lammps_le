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

/* ----------------------------------------------------------------------
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL),
                         Evan Weinberg (NVIDIA)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_snap_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "kokkos.h"
#include "sna.h"

#define MAXLINE 1024
#define MAXWORD 3

namespace LAMMPS_NS {

// Outstanding issues with quadratic term
// 1. there seems to a problem with compute_optimized energy calc
// it does not match compute_regular, even when quadratic coeffs = 0

//static double t1 = 0.0;
//static double t2 = 0.0;
//static double t3 = 0.0;
//static double t4 = 0.0;
//static double t5 = 0.0;
//static double t6 = 0.0;
//static double t7 = 0.0;
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::PairSNAPKokkos(LAMMPS *lmp) : PairSNAP(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_cutsq = tdual_fparams("PairSNAPKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;

  host_flag = (execution_space == Host);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::~PairSNAPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->request(this,instance_me);

  neighbor->requests[irequest]->
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) { // still need atomics, even though using a full neigh list
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Must use half neighbor list style with pair snap/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& max_neighs) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }
};

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  int newton_pair = force->newton_pair;
  if (newton_pair == false)
    error->all(FLERR,"PairSNAPKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  /*
  for (int i = 0; i < nlocal; i++) {
    typename t_neigh_list::t_neighs neighs_i = neigh_list.get_neighs(i);
    const int num_neighs = neighs_i.get_num_neighs();
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }*/
  max_neighs = 0;
  Kokkos::parallel_reduce("PairSNAPKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(max_neighs));

  int vector_length_default = 1;
  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 32;//max_neighs;

  if (beta_max < inum) {
    beta_max = inum;
    d_beta = Kokkos::View<F_FLOAT**, DeviceType>("PairSNAPKokkos:beta",ncoeff,inum);
    if (!host_flag)
      d_beta_pack = Kokkos::View<F_FLOAT***, Kokkos::LayoutLeft, DeviceType>("PairSNAPKokkos:beta_pack",32,ncoeff,(inum+32-1)/32);
    d_ninside = Kokkos::View<int*, DeviceType>("PairSNAPKokkos:ninside",inum);
  }

  chunk_size = MIN(chunksize,inum); // "chunksize" variable is set by user
  chunk_offset = 0;

  snaKK.grow_rij(chunk_size,max_neighs);

  EV_FLOAT ev;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    //ComputeNeigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairSNAPComputeNeigh>(chunk_size,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeNeigh> policy_neigh(chunk_size,team_size,vector_length);
      Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
    }

    if (host_flag)
    {
      // Host codepath

      //PreUi
      {
        int vector_length = vector_length_default;
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPPreUiCPU>(chunk_size,team_size,vector_length);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUiCPU> policy_preui_cpu((chunk_size+team_size-1)/team_size,team_size,vector_length);
        Kokkos::parallel_for("PreUiCPU",policy_preui_cpu,*this);
      }

      // ComputeUi
      {
        int vector_length = vector_length_default;
        int team_size = team_size_default;
        // Fused calculation of ulist and accumulation into ulisttot using atomics
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiCPU> policy_ui_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        Kokkos::parallel_for("ComputeUiCPU",policy_ui_cpu,*this);
      }

      {
        // Expand ulisttot -> ulisttot_full
        // Zero out ylist
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPTransformUiCPU> policy_transform_ui_cpu({0,0},{twojmax+1,chunk_size});
        Kokkos::parallel_for("TransformUiCPU",policy_transform_ui_cpu,*this);
      }

      //Compute bispectrum
      if (quadraticflag || eflag) {
        //ComputeZi
        int idxz_max = snaKK.idxz_max;
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeZiCPU> policy_zi_cpu(0,chunk_size*idxz_max);
        Kokkos::parallel_for("ComputeZiCPU",policy_zi_cpu,*this);

        //ComputeBi
        int vector_length = vector_length_default;
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPComputeBiCPU>(chunk_size,team_size,vector_length);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeBiCPU> policy_bi_cpu(chunk_size,team_size,vector_length);
        Kokkos::parallel_for("ComputeBiCPU",policy_bi_cpu,*this);
      }

      //ComputeYi
      {
        //Compute beta = dE_i/dB_i for all i in list
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPBetaCPU> policy_beta(0,chunk_size);
        Kokkos::parallel_for("ComputeBetaCPU",policy_beta,*this);

        //ComputeYi
        int idxz_max = snaKK.idxz_max;
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeYiCPU> policy_yi_cpu(0,chunk_size*idxz_max);
        Kokkos::parallel_for("ComputeYiCPU",policy_yi_cpu,*this);
      } // host flag

      //ComputeDuidrj and Deidrj
      {
        int team_size = team_size_default;
        int vector_length = vector_length_default;

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDuidrjCPU> policy_duidrj_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        snaKK.set_dir(-1); // technically doesn't do anything
        Kokkos::parallel_for("ComputeDuidrjCPU",policy_duidrj_cpu,*this);

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDeidrjCPU> policy_deidrj_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);

        Kokkos::parallel_for("ComputeDeidrjCPU",policy_deidrj_cpu,*this);
      }

    } else { // GPU

#ifdef LMP_KOKKOS_GPU
      //PreUi
      {
        int vector_length = vector_length_default;
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPPreUi>(chunk_size,team_size,vector_length);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUi> policy_preui((chunk_size+team_size-1)/team_size,team_size,vector_length);
        Kokkos::parallel_for("PreUi",policy_preui,*this);
      }

      // ComputeUi w/vector parallelism, shared memory, direct atomicAdd into ulisttot
      {

        int vector_length = 32;
        int team_size = 4; // need to cap b/c of shared memory reqs
        check_team_size_for<TagPairSNAPComputeUi>(chunk_size,team_size,vector_length);

        // scratch size: 2 * team_size * (twojmax+1)^2, to cover all `m1`,`m2` values, div 2 for symmetry
        //   2 is for double buffer

        const int tile_size = (twojmax+1)*(twojmax/2+1);
        typedef Kokkos::View< SNAcomplex*,
                              Kokkos::DefaultExecutionSpace::scratch_memory_space,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged> >
                ScratchViewType;
        int scratch_size = ScratchViewType::shmem_size( 2 * team_size * tile_size );

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUi> policy_ui(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam( scratch_size ));

        Kokkos::parallel_for("ComputeUi",policy_ui,*this);

        //Transform data layout of ulisttot to AoSoA, zero ylist
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPTransformUi> policy_transform_ui({0,0,0},{32,twojmax+1,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("TransformUi",policy_transform_ui,*this);

      }

      //Compute bispectrum in AoSoA data layout, transform Bi
      if (quadraticflag || eflag) {
        //ComputeZi
        int idxz_max = snaKK.idxz_max;
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPComputeZi> policy_compute_zi({0,0,0},{32,idxz_max,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("ComputeZi",policy_compute_zi,*this);

        //ComputeBi
        int idxb_max = snaKK.idxb_max;
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPComputeBi> policy_compute_bi({0,0,0},{32,idxb_max,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("ComputeBi",policy_compute_bi,*this);

        //Transform data layout of blist out of AoSoA
        //We need this b/c `blist` gets used in ComputeForce which doesn't
        //take advantage of AoSoA (which at best would only be beneficial
        //on the margins)
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPTransformBi> policy_transform_bi({0,0,0},{32,idxb_max,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("TransformBi",policy_transform_bi,*this);
      }

      //ComputeYi in AoSoA data layout, transform to AoS for ComputeFusedDeidrj
      //Note zeroing `ylist` is fused into `TransformUi`.
      {
        //Compute beta = dE_i/dB_i for all i in list
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPBeta> policy_beta(0,chunk_size);
        Kokkos::parallel_for("ComputeBeta",policy_beta,*this);

        //ComputeYi
        const int idxz_max = snaKK.idxz_max;
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPComputeYi> policy_compute_yi({0,0,0},{32,idxz_max,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("ComputeYi",policy_compute_yi,*this);

        //Transform data layout of ylist out of AoSoA
        const int idxu_half_max = snaKK.idxu_half_max;
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPTransformYi> policy_transform_yi({0,0,0},{32,idxu_half_max,(chunk_size + 32 - 1) / 32},{32,4,1});
        Kokkos::parallel_for("TransformYi",policy_transform_yi,*this);

      }

      // Fused ComputeDuidrj, ComputeDeidrj
      {
        int vector_length = 32;
        int team_size = 2; // need to cap b/c of shared memory reqs
        check_team_size_for<TagPairSNAPComputeFusedDeidrj>(chunk_size,team_size,vector_length);

        // scratch size: 2 * 2 * team_size * (twojmax+1)*(twojmax/2+1), to cover half `m1`,`m2` values due to symmetry
        // 2 is for double buffer
        const int tile_size = (twojmax+1)*(twojmax/2+1);

        typedef Kokkos::View< SNAcomplex*,
                              Kokkos::DefaultExecutionSpace::scratch_memory_space,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged> >
                ScratchViewType;
        int scratch_size = ScratchViewType::shmem_size( 4 * team_size * tile_size);

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeFusedDeidrj> policy_fused_deidrj(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        policy_fused_deidrj = policy_fused_deidrj.set_scratch_size(0, Kokkos::PerTeam( scratch_size ));

        for (int k = 0; k < 3; k++) {
          snaKK.set_dir(k);
          Kokkos::parallel_for("ComputeFusedDeidrj",policy_fused_deidrj,*this);
        }
      }

#endif // LMP_KOKKOS_GPU

    }

    //ComputeForce
    {
      int team_size = team_size_default;
      int vector_length = vector_length_default;
      if (evflag) {
        if (neighflag == HALF) {
          check_team_size_reduce<TagPairSNAPComputeForce<HALF,1> >(chunk_size,team_size,vector_length);
          typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<HALF,1> > policy_force(chunk_size,team_size,vector_length);
          Kokkos::parallel_reduce(policy_force
            ,*this,ev_tmp);
        } else if (neighflag == HALFTHREAD) {
          check_team_size_reduce<TagPairSNAPComputeForce<HALFTHREAD,1> >(chunk_size,team_size,vector_length);
          typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<HALFTHREAD,1> > policy_force(chunk_size,team_size,vector_length);
          Kokkos::parallel_reduce(policy_force
            ,*this,ev_tmp);
        }
      } else {
        if (neighflag == HALF) {
          check_team_size_for<TagPairSNAPComputeForce<HALF,0> >(chunk_size,team_size,vector_length);
          typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<HALF,0> > policy_force(chunk_size,team_size,vector_length);
          Kokkos::parallel_for(policy_force
            ,*this);
        } else if (neighflag == HALFTHREAD) {
          check_team_size_for<TagPairSNAPComputeForce<HALFTHREAD,0> >(chunk_size,team_size,vector_length);
          typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<HALFTHREAD,0> > policy_force(chunk_size,team_size,vector_length);
          Kokkos::parallel_for(policy_force
            ,*this);
        }
      }
    }
    ev += ev_tmp;
    chunk_offset += chunk_size;

  } // end while

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::allocate()
{
  PairSNAP::allocate();

  int n = atom->ntypes;
  d_map = Kokkos::View<T_INT*, DeviceType>("PairSNAPKokkos::map",n+1);
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairSNAP::init_one(i,j);
  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairSNAP::coeff(narg,arg);

  // Set up element lists

  d_radelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:radelem",nelements);
  d_wjelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:wjelem",nelements);
  d_coeffelem = Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType>("pair:coeffelem",nelements,ncoeffall);

  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_wjelem = Kokkos::create_mirror_view(d_wjelem);
  auto h_coeffelem = Kokkos::create_mirror_view(d_coeffelem);
  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int ielem = 0; ielem < nelements; ielem++) {
    h_radelem(ielem) = radelem[ielem];
    h_wjelem(ielem) = wjelem[ielem];
    for (int jcoeff = 0; jcoeff < ncoeffall; jcoeff++) {
      h_coeffelem(ielem,jcoeff) = coeffelem[ielem][jcoeff];
    }
  }

  for (int i = 1; i <= atom->ntypes; i++) {
    h_map(i) = map[i];
  }

  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_wjelem,h_wjelem);
  Kokkos::deep_copy(d_coeffelem,h_coeffelem);
  Kokkos::deep_copy(d_map,h_map);

  snaKK = SNAKokkos<DeviceType>(rfac0,twojmax,
                  rmin0,switchflag,bzeroflag,chemflag,bnormflag,wselfallflag,nelements);
  snaKK.grow_rij(0,0);
  snaKK.init();
}

/* ----------------------------------------------------------------------
   Begin routines that are called on both CPU and GPU codepaths
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeNeigh>::member_type& team) const {

  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  SNAKokkos<DeviceType> my_sna = snaKK;
  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);
  const int itype = type[i];
  const int ielem = d_map[itype];
  const double radi = d_radelem[ielem];

  const int num_neighs = d_numneigh[i];

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // wj = weights for neighbors of I within cutoff
  // rcutij = cutoffs for neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  int ninside = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,num_neighs),
      [&] (const int jj, int& count) {
    Kokkos::single(Kokkos::PerThread(team), [&] (){
      T_INT j = d_neighbors(i,jj);
      const F_FLOAT dx = x(j,0) - xtmp;
      const F_FLOAT dy = x(j,1) - ytmp;
      const F_FLOAT dz = x(j,2) - ztmp;

      const int jtype = type(j);
      const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;

      if ( rsq < rnd_cutsq(itype,jtype) )
       count++;
    });
  },ninside);

  d_ninside(ii) = ninside;

  if (team.team_rank() == 0)
  Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_neighs),
      [&] (const int jj, int& offset, bool final) {
  //for (int jj = 0; jj < num_neighs; jj++) {
    T_INT j = d_neighbors(i,jj);
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;

    const int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
    const int elem_j = d_map[jtype];

    if ( rsq < rnd_cutsq(itype,jtype) ) {
      if (final) {
#ifdef LMP_KOKKOS_GPU
        if (!host_flag) {
          my_sna.compute_cayley_klein(ii, offset, dx, dy, dz, (radi + d_radelem[elem_j])*rcutfac,
                                      d_wjelem[elem_j]);
        } else {
#endif
          my_sna.rij(ii,offset,0) = dx;
          my_sna.rij(ii,offset,1) = dy;
          my_sna.rij(ii,offset,2) = dz;
          my_sna.wj(ii,offset) = d_wjelem[elem_j];
          my_sna.rcutij(ii,offset) = (radi + d_radelem[elem_j])*rcutfac;
#ifdef LMP_KOKKOS_GPU
        }
#endif
        my_sna.inside(ii,offset) = j;
        if (chemflag)
          my_sna.element(ii,offset) = elem_j;
        else
          my_sna.element(ii,offset) = 0;
      }
      offset++;
    }
  });
}

/* ----------------------------------------------------------------------
   Begin routines that are unique to the GPU codepath. These take advantage
   of AoSoA data layouts and scratch memory for recursive polynomials
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPBeta,const int& ii) const {

  if (ii >= chunk_size) return;

  const int iatom_mod = ii % 32;
  const int iatom_div = ii / 32;

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const int ielem = d_map[itype];
  SNAKokkos<DeviceType> my_sna = snaKK;

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    d_coeffi(d_coeffelem,ielem,Kokkos::ALL);

  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
    d_beta_pack(iatom_mod,icoeff,iatom_div) = d_coeffi[icoeff+1];
  }

  if (quadraticflag) {
    const auto idxb_max = my_sna.idxb_max;
    int k = ncoeff+1;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      const auto idxb = icoeff % idxb_max;
      const auto idx_chem = icoeff / idxb_max;
      double bveci = my_sna.blist(idxb, idx_chem, ii);
      d_beta_pack(iatom_mod,icoeff,iatom_div) += d_coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
        const auto jdxb = jcoeff % idxb_max;
        const auto jdx_chem = jcoeff / idxb_max;
        double bvecj = my_sna.blist(jdxb, jdx_chem, ii);
        d_beta_pack(iatom_mod,icoeff,iatom_div) += d_coeffi[k]*bvecj;
        d_beta_pack(iatom_mod,jcoeff,iatom_div) += d_coeffi[k]*bveci;
        k++;
      }
    }
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPPreUi,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUi>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  const int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  int itype = type(ii);
  int ielem = d_map[itype];

  my_sna.pre_ui(team,ii,ielem);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeUi,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUi>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_ui(team,ii,jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPTransformUi,const int iatom_mod, const int j, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (j > twojmax) return;

  int elem_count = chemflag ? nelements : 1;

  for (int ielem = 0; ielem < elem_count; ielem++) {
    const int jju_half = my_sna.idxu_half_block(j);
    const int jju = my_sna.idxu_block(j);

    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        // Extract top half

        const int idxu_shift = mb * (j + 1) + ma;
        const int idxu_half = jju_half + idxu_shift;
        const int idxu = jju + idxu_shift;

        auto utot_re = my_sna.ulisttot_re(idxu_half, ielem, iatom);
        auto utot_im = my_sna.ulisttot_im(idxu_half, ielem, iatom);

        // Store
        my_sna.ulisttot_pack(iatom_mod, idxu, ielem, iatom_div) = { utot_re, utot_im };

        // Also zero yi
        my_sna.ylist_pack_re(iatom_mod, idxu_half, ielem, iatom_div) = 0.;
        my_sna.ylist_pack_im(iatom_mod, idxu_half, ielem, iatom_div) = 0.;

        // Symmetric term
        const int sign_factor = (((ma+mb)%2==0)?1:-1);
        const int idxu_flip = jju + (j + 1 - mb) * (j + 1) - (ma + 1);

        if (sign_factor == 1) {
          utot_im = -utot_im;
        } else {
          utot_re = -utot_re;
        }

        my_sna.ulisttot_pack(iatom_mod, idxu_flip, ielem, iatom_div) = { utot_re, utot_im };

        // No need to zero symmetrized ylist
      }
    }
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeYi,const int iatom_mod, const int jjz, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (jjz >= my_sna.idxz_max) return;

  my_sna.compute_yi(iatom_mod,jjz,iatom_div,d_beta_pack);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPTransformYi,const int iatom_mod, const int idxu_half, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (idxu_half >= my_sna.idxu_half_max) return;

  int elem_count = chemflag ? nelements : 1;
  for (int ielem = 0; ielem < elem_count; ielem++) {
    const auto y_re = my_sna.ylist_pack_re(iatom_mod, idxu_half, ielem, iatom_div);
    const auto y_im = my_sna.ylist_pack_im(iatom_mod, idxu_half, ielem, iatom_div);

    my_sna.ylist(idxu_half, ielem, iatom) = { y_re, y_im };
  }

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeZi,const int iatom_mod, const int jjz, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (jjz >= my_sna.idxz_max) return;

  my_sna.compute_zi(iatom_mod,jjz,iatom_div);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeBi,const int iatom_mod, const int jjb, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (jjb >= my_sna.idxb_max) return;

  my_sna.compute_bi(iatom_mod,jjb,iatom_div);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPTransformBi,const int iatom_mod, const int idxb, const int iatom_div) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * 32;
  if (iatom >= chunk_size) return;

  if (idxb >= my_sna.idxb_max) return;

  const int ntriples = my_sna.ntriples;

  for (int itriple = 0; itriple < ntriples; itriple++) {

    const auto blocal = my_sna.blist_pack(iatom_mod, idxb, itriple, iatom_div);

    my_sna.blist(idxb, itriple, iatom) = blocal;
  }

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeFusedDeidrj,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeFusedDeidrj>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_fused_deidrj(team,ii,jj);
}

/* ----------------------------------------------------------------------
   Begin routines that are unique to the CPU codepath. These do not take
   advantage of AoSoA data layouts, but that could be a good point of
   future optimization and unification with the above kernels. It's unlikely
   that scratch memory optimizations will ever be useful for the CPU due to
   different arithmetic intensity requirements for the CPU vs GPU.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPBetaCPU,const int& ii) const {

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const int ielem = d_map[itype];
  SNAKokkos<DeviceType> my_sna = snaKK;

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    d_coeffi(d_coeffelem,ielem,Kokkos::ALL);

  for (int icoeff = 0; icoeff < ncoeff; icoeff++)
    d_beta(icoeff,ii) = d_coeffi[icoeff+1];

  if (quadraticflag) {
    const auto idxb_max = my_sna.idxb_max;
    int k = ncoeff+1;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      const auto idxb = icoeff % idxb_max;
      const auto idx_chem = icoeff / idxb_max;
      double bveci = my_sna.blist(idxb,idx_chem,ii);
      d_beta(icoeff,ii) += d_coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
        const auto jdxb = jcoeff % idxb_max;
        const auto jdx_chem = jcoeff / idxb_max;
        double bvecj = my_sna.blist(jdxb,jdx_chem,ii);
        d_beta(icoeff,ii) += d_coeffi[k]*bvecj;
        d_beta(jcoeff,ii) += d_coeffi[k]*bveci;
        k++;
      }
    }
  }
}


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPPreUiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUiCPU>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  const int ii = team.team_rank() + team.team_size() * team.league_rank();
  if (ii >= chunk_size) return;
  int itype = type(ii);
  int ielem = d_map[itype];

  my_sna.pre_ui_cpu(team,ii,ielem);
}



template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeUiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiCPU>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_ui_cpu(team,ii,jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPTransformUiCPU, const int j, const int iatom) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  if (iatom >= chunk_size) return;

  if (j > twojmax) return;

  int elem_count = chemflag ? nelements : 1;

  // De-symmetrize ulisttot
  for (int ielem = 0; ielem < elem_count; ielem++) {

    const int jju_half = my_sna.idxu_half_block(j);
    const int jju = my_sna.idxu_block(j);

    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        // Extract top half

        const int idxu_shift = mb * (j + 1) + ma;
        const int idxu_half = jju_half + idxu_shift;
        const int idxu = jju + idxu_shift;

        // Load ulist
        auto utot = my_sna.ulisttot(idxu_half, ielem, iatom);

        // Store
        my_sna.ulisttot_full(idxu, ielem, iatom) = utot;

        // Zero Yi
        my_sna.ylist(idxu_half, ielem, iatom) = {0., 0.};

        // Symmetric term
        const int sign_factor = (((ma+mb)%2==0)?1:-1);
        const int idxu_flip = jju + (j + 1 - mb) * (j + 1) - (ma + 1);

        if (sign_factor == 1) {
          utot.im = -utot.im;
        } else {
          utot.re = -utot.re;
        }

        my_sna.ulisttot_full(idxu_flip, ielem, iatom) = utot;
      }
    }
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeYiCPU,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_yi_cpu(ii,d_beta);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeZiCPU,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_zi_cpu(ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeBiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeBiCPU>::member_type& team) const {
  int ii = team.league_rank();
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_bi_cpu(team,ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeDuidrjCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDuidrjCPU>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_duidrj_cpu(team,ii,jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeDeidrjCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDeidrjCPU>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_deidrj_cpu(team,ii,jj);
}

/* ----------------------------------------------------------------------
   Also used for both CPU and GPU codepaths. Could maybe benefit from a
   separate GPU/CPU codepath, but this kernel takes so little time it's
   likely not worth it.
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  SNAKokkos<DeviceType> my_sna = snaKK;
  const int ninside = d_ninside(ii);

  Kokkos::parallel_for (Kokkos::TeamThreadRange(team,ninside),
      [&] (const int jj) {
    int j = my_sna.inside(ii,jj);

    F_FLOAT fij[3];
    fij[0] = my_sna.dedr(ii,jj,0);
    fij[1] = my_sna.dedr(ii,jj,1);
    fij[2] = my_sna.dedr(ii,jj,2);

    Kokkos::single(Kokkos::PerThread(team), [&] (){
      a_f(i,0) += fij[0];
      a_f(i,1) += fij[1];
      a_f(i,2) += fij[2];
      a_f(j,0) -= fij[0];
      a_f(j,1) -= fij[1];
      a_f(j,2) -= fij[2];

      // tally global and per-atom virial contribution

      if (EVFLAG) {
        if (vflag_either) {
          v_tally_xyz<NEIGHFLAG>(ev,i,j,
            fij[0],fij[1],fij[2],
            -my_sna.rij(ii,jj,0),-my_sna.rij(ii,jj,1),
            -my_sna.rij(ii,jj,2));
        }
      }

    });
  });

  // tally energy contribution

  if (EVFLAG) {
    if (eflag_either) {

      const int itype = type(i);
      const int ielem = d_map[itype];
      Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        d_coeffi(d_coeffelem,ielem,Kokkos::ALL);

      Kokkos::single(Kokkos::PerTeam(team), [&] () {

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k

        double evdwl = d_coeffi[0];

        // E = beta.B + 0.5*B^t.alpha.B

        const auto idxb_max = snaKK.idxb_max;

        // linear contributions

        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          const auto idxb = icoeff % idxb_max;
          const auto idx_chem = icoeff / idxb_max;
          evdwl += d_coeffi[icoeff+1]*my_sna.blist(idxb,idx_chem,ii);
        }

        // quadratic contributions

        if (quadraticflag) {
          int k = ncoeff+1;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            const auto idxb = icoeff % idxb_max;
            const auto idx_chem = icoeff / idxb_max;
            double bveci = my_sna.blist(idxb,idx_chem,ii);

            evdwl += 0.5*d_coeffi[k++]*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              auto jdxb = jcoeff % idxb_max;
              auto jdx_chem = jcoeff / idxb_max;
              double bvecj = my_sna.blist(jdxb,jdx_chem,ii);

              evdwl += d_coeffi[k++]*bveci*bvecj;
            }
          }
        }

        //ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
        if (eflag_global) ev.evdwl += evdwl;
        if (eflag_atom) d_eatom[i] += evdwl;
      });
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>(), team, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const E_FLOAT v0 = delx*fx;
  const E_FLOAT v1 = dely*fy;
  const E_FLOAT v2 = delz*fz;
  const E_FLOAT v3 = delx*fy;
  const E_FLOAT v4 = delx*fz;
  const E_FLOAT v5 = dely*fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i,0) += 0.5*v0;
    a_vatom(i,1) += 0.5*v1;
    a_vatom(i,2) += 0.5*v2;
    a_vatom(i,3) += 0.5*v3;
    a_vatom(i,4) += 0.5*v4;
    a_vatom(i,5) += 0.5*v5;
    a_vatom(j,0) += 0.5*v0;
    a_vatom(j,1) += 0.5*v1;
    a_vatom(j,2) += 0.5*v2;
    a_vatom(j,3) += 0.5*v3;
    a_vatom(j,4) += 0.5*v4;
    a_vatom(j,5) += 0.5*v5;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);
  bytes += n*n*sizeof(double);
  bytes += (2*ncoeffall)*sizeof(double);
  bytes += (ncoeff*3)*sizeof(double);
  bytes += snaKK.memory_usage();
  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairSNAPKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if(team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType>
template<class TagStyle>
void PairSNAPKokkos<DeviceType>::check_team_size_reduce(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());

  if(team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

}
