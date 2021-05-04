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

PairStyle(snap/kk,PairSNAPKokkos<LMPDeviceType>)
PairStyle(snap/kk/device,PairSNAPKokkos<LMPDeviceType>)
PairStyle(snap/kk/host,PairSNAPKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_SNAP_KOKKOS_H
#define LMP_PAIR_SNAP_KOKKOS_H

#include "pair_snap.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "sna_kokkos.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

// Routines for both the CPU and GPU backend
template<int NEIGHFLAG, int EVFLAG>
struct TagPairSNAPComputeForce{};
struct TagPairSNAPComputeNeigh{};

// GPU backend only
struct TagPairSNAPPreUi{};
struct TagPairSNAPComputeUi{};
struct TagPairSNAPTransformUi{}; // re-order ulisttot from SoA to AoSoA, zero ylist
struct TagPairSNAPComputeZi{};
struct TagPairSNAPBeta{};
struct TagPairSNAPComputeBi{};
struct TagPairSNAPTransformBi{}; // re-order blist from AoSoA to AoS
struct TagPairSNAPComputeYi{};
struct TagPairSNAPTransformYi{}; // re-order ylist from AoSoA to AoS
struct TagPairSNAPComputeFusedDeidrj{};

// CPU backend only
struct TagPairSNAPPreUiCPU{};
struct TagPairSNAPComputeUiCPU{};
struct TagPairSNAPTransformUiCPU{};
struct TagPairSNAPComputeZiCPU{};
struct TagPairSNAPBetaCPU{};
struct TagPairSNAPComputeBiCPU{};
struct TagPairSNAPZeroYiCPU{};
struct TagPairSNAPComputeYiCPU{};
struct TagPairSNAPComputeDuidrjCPU{};
struct TagPairSNAPComputeDeidrjCPU{};

template<class DeviceType>
class PairSNAPKokkos : public PairSNAP {
public:
  enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairSNAPKokkos(class LAMMPS *);
  ~PairSNAPKokkos();

  void coeff(int, char**);
  void init_style();
  double init_one(int, int);
  void compute(int, int);
  double memory_usage();

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template<class TagStyle>
  void check_team_size_reduce(int, int&, int);

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team, EV_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPBetaCPU,const int& ii) const;

  // GPU backend only
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPPreUi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformUi,const int iatom_mod, const int j, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZi,const int iatom_mod, const int idxz, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPBeta, const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBi,const int iatom_mod, const int idxb, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformBi,const int iatom_mod, const int idxb, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYi,const int iatom_mod, const int idxz, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformYi,const int iatom_mod, const int idxu, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrj,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrj>::member_type& team) const;

  // CPU backend only
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUiCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPPreUiCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUiCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformUiCPU, const int j, const int iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZiCPU,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBiCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeBiCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYiCPU,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDuidrjCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDuidrjCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDeidrjCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDeidrjCPU>::member_type& team) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const;

protected:
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  typedef Kokkos::View<F_FLOAT**> t_bvec;
  t_bvec bvec;
  typedef Kokkos::View<F_FLOAT***> t_dbvec;
  t_dbvec dbvec;
  SNAKokkos<DeviceType> snaKK;

  int inum,max_neighs,chunk_size,chunk_offset;
  int host_flag;

  int eflag,vflag;

  void allocate();
  //void read_files(char *, char *);
  /*template<class DeviceType>
inline int equal(double* x,double* y);
  template<class DeviceType>
inline double dist2(double* x,double* y);
  double extra_cutoff();
  void load_balance();
  void set_sna_to_shared(int snaid,int i);
  void build_per_atom_arrays();*/

  int neighflag;

  Kokkos::View<T_INT*, DeviceType> ilistmast;
  Kokkos::View<T_INT*, DeviceType> ghostilist;
  Kokkos::View<T_INT*, DeviceType> ghostnumneigh;
  Kokkos::View<T_INT*, DeviceType> ghostneighs;
  Kokkos::View<T_INT*, DeviceType> ghostfirstneigh;

  Kokkos::View<T_INT**, Kokkos::LayoutRight, DeviceType> i_pairs;
  Kokkos::View<T_INT***, Kokkos::LayoutRight, DeviceType> i_rij;
  Kokkos::View<T_INT**, Kokkos::LayoutRight, DeviceType> i_inside;
  Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType> i_wj;
  Kokkos::View<F_FLOAT***, Kokkos::LayoutRight, DeviceType>i_rcutij;
  Kokkos::View<T_INT*, DeviceType> i_ninside;
  Kokkos::View<F_FLOAT****, Kokkos::LayoutRight, DeviceType> i_uarraytot_r, i_uarraytot_i;
  Kokkos::View<F_FLOAT******, Kokkos::LayoutRight, DeviceType> i_zarray_r, i_zarray_i;

  Kokkos::View<F_FLOAT*, DeviceType> d_radelem;              // element radii
  Kokkos::View<F_FLOAT*, DeviceType> d_wjelem;               // elements weights
  Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType> d_coeffelem;           // element bispectrum coefficients
  Kokkos::View<T_INT*, DeviceType> d_map;                    // mapping from atom types to elements
  Kokkos::View<T_INT*, DeviceType> d_ninside;                // ninside for all atoms in list
  Kokkos::View<F_FLOAT**, DeviceType> d_beta;                // betas for all atoms in list
  Kokkos::View<F_FLOAT***, Kokkos::LayoutLeft, DeviceType> d_beta_pack;          // betas for all atoms in list, GPU
  Kokkos::View<F_FLOAT**, DeviceType> d_bispectrum;          // bispectrum components for all atoms in list

  typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<const F_FLOAT**, DeviceType,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_fparams_rnd;
  t_fparams_rnd rnd_cutsq;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,typename Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_f;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,typename Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_vatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,typename Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_f;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,typename Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_vatom;

  friend void pair_virial_fdotr_compute<PairSNAPKokkos>(PairSNAPKokkos*);

};

}

#endif
#endif
