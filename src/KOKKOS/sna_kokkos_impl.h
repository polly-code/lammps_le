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
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "sna_kokkos.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <type_traits>

namespace LAMMPS_NS {

static const double MY_PI  = 3.14159265358979323846; // pi

template<class DeviceType>
inline
SNAKokkos<DeviceType>::SNAKokkos(double rfac0_in,
         int twojmax_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
         int chem_flag_in, int bnorm_flag_in, int wselfall_flag_in, int nelements_in)
{
  LAMMPS_NS::ExecutionSpace execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  host_flag = (execution_space == LAMMPS_NS::Host);

  wself = 1.0;

  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;
  bzero_flag = bzero_flag_in;

  chem_flag = chem_flag_in;
  if (chem_flag)
    nelements = nelements_in;
  else
    nelements = 1;
  bnorm_flag = bnorm_flag_in;
  wselfall_flag = wselfall_flag_in;

  twojmax = twojmax_in;

  ncoeff = compute_ncoeff();

  nmax = 0;

  build_indexlist();

  int jdimpq = twojmax + 2;
  rootpqarray = t_sna_2d("SNAKokkos::rootpqarray",jdimpq,jdimpq);

  cglist = t_sna_1d("SNAKokkos::cglist",idxcg_max);

  if (bzero_flag) {
    bzero = Kokkos::View<double*, Kokkos::LayoutRight, DeviceType>("sna:bzero",twojmax+1);
    auto h_bzero = Kokkos::create_mirror_view(bzero);

    double www = wself*wself*wself;
    for(int j = 0; j <= twojmax; j++)
      if (bnorm_flag)
        h_bzero[j] = www;
      else
        h_bzero[j] = www*(j+1);
    Kokkos::deep_copy(bzero,h_bzero);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
SNAKokkos<DeviceType>::~SNAKokkos()
{
}

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::build_indexlist()
{
  // index list for cglist

  int jdim = twojmax + 1;
  idxcg_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxcg_block",jdim,jdim,jdim);
  auto h_idxcg_block = Kokkos::create_mirror_view(idxcg_block);

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxcg_block(j1,j2,j) = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;
  Kokkos::deep_copy(idxcg_block,h_idxcg_block);

  // index list for uarray
  // need to include both halves

  idxu_block = Kokkos::View<int*, DeviceType>("SNAKokkos::idxu_block",jdim);
  auto h_idxu_block = Kokkos::create_mirror_view(idxu_block);

  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    h_idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;
  Kokkos::deep_copy(idxu_block,h_idxu_block);

  // index list for half uarray
  idxu_half_block = Kokkos::View<int*, DeviceType>("SNAKokkos::idxu_half_block",jdim);
  auto h_idxu_half_block = Kokkos::create_mirror_view(idxu_half_block);

  int idxu_half_count = 0;
  for(int j = 0; j <= twojmax; j++) {
    h_idxu_half_block[j] = idxu_half_count;
    for(int mb = 0; 2*mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_half_count++;
  }
  idxu_half_max = idxu_half_count;
  Kokkos::deep_copy(idxu_half_block, h_idxu_half_block);

  // index list for "cache" uarray
  // this is the GPU scratch memory requirements
  // applied the CPU structures
  idxu_cache_block = Kokkos::View<int*, DeviceType>("SNAKokkos::idxu_cache_block",jdim);
  auto h_idxu_cache_block = Kokkos::create_mirror_view(idxu_cache_block);

  int idxu_cache_count = 0;
  for(int j = 0; j <= twojmax; j++) {
    h_idxu_cache_block[j] = idxu_cache_count;
    for(int mb = 0; mb < ((j+3)/2); mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_cache_count++;
  }
  idxu_cache_max = idxu_cache_count;
  Kokkos::deep_copy(idxu_cache_block, h_idxu_cache_block);

  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
  idxb = Kokkos::View<int*[3], DeviceType>("SNAKokkos::idxb",idxb_max);
  auto h_idxb = Kokkos::create_mirror_view(idxb);

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          h_idxb(idxb_count,0) = j1;
          h_idxb(idxb_count,1) = j2;
          h_idxb(idxb_count,2) = j;
          idxb_count++;
        }
  Kokkos::deep_copy(idxb,h_idxb);

  // reverse index list for beta and b

  idxb_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxb_block",jdim,jdim,jdim);
  auto h_idxb_block = Kokkos::create_mirror_view(idxb_block);

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          h_idxb_block(j1,j2,j) = idxb_count;
          idxb_count++;
        }
      }
  Kokkos::deep_copy(idxb_block,h_idxb_block);

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  idxz = Kokkos::View<int*[10], DeviceType>("SNAKokkos::idxz",idxz_max);
  auto h_idxz = Kokkos::create_mirror_view(idxz);

  idxz_block = Kokkos::View<int***, DeviceType>("SNAKokkos::idxz_block", jdim,jdim,jdim);
  auto h_idxz_block = Kokkos::create_mirror_view(idxz_block);

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxz_block(j1,j2,j) = idxz_count;

        // find right beta(ii,jjb) entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            h_idxz(idxz_count,0) = j1;
            h_idxz(idxz_count,1) = j2;
            h_idxz(idxz_count,2) = j;
            h_idxz(idxz_count,3) = MAX(0, (2 * ma - j - j2 + j1) / 2);
            h_idxz(idxz_count,4) = (2 * ma - j - (2 * h_idxz(idxz_count,3) - j1) + j2) / 2;
            h_idxz(idxz_count,5) = MAX(0, (2 * mb - j - j2 + j1) / 2);
            h_idxz(idxz_count,6) = (2 * mb - j - (2 * h_idxz(idxz_count,5) - j1) + j2) / 2;
            h_idxz(idxz_count,7) = MIN(j1, (2 * ma - j + j2 + j1) / 2) - h_idxz(idxz_count,3) + 1;
            h_idxz(idxz_count,8) = MIN(j1, (2 * mb - j + j2 + j1) / 2) - h_idxz(idxz_count,5) + 1;

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
            // ylist is "compressed" via symmetry in its
            // contraction with dulist
            const int jju_half = h_idxu_half_block[j] + (j+1)*mb + ma;
            h_idxz(idxz_count,9) = jju_half;

            idxz_count++;
          }
      }
  Kokkos::deep_copy(idxz,h_idxz);
  Kokkos::deep_copy(idxz_block,h_idxz_block);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::grow_rij(int newnatom, int newnmax)
{
  if(newnatom <= natom && newnmax <= nmax) return;
  natom = newnatom;
  nmax = newnmax;

  inside = t_sna_2i(Kokkos::NoInit("sna:inside"),natom,nmax);
  element = t_sna_2i(Kokkos::NoInit("sna:rcutij"),natom,nmax);
  dedr = t_sna_3d(Kokkos::NoInit("sna:dedr"),natom,nmax,3);

#ifdef LMP_KOKKOS_GPU
  if (!host_flag) {

    cayleyklein = t_sna_2ckp(Kokkos::NoInit("sna:cayleyklein"), natom, nmax);
    ulisttot = t_sna_3c_ll(Kokkos::NoInit("sna:ulisttot"),1,1,1); // dummy allocation
    ulisttot_full = t_sna_3c_ll(Kokkos::NoInit("sna:ulisttot"),1,1,1);
    ulisttot_re = t_sna_3d_ll(Kokkos::NoInit("sna:ulisttot_re"),idxu_half_max,nelements,natom);
    ulisttot_im = t_sna_3d_ll(Kokkos::NoInit("sna:ulisttot_im"),idxu_half_max,nelements,natom);
    ulisttot_pack = t_sna_4c_ll(Kokkos::NoInit("sna:ulisttot_pack"),32,idxu_max,nelements,(natom+32-1)/32);
    ulist = t_sna_3c_ll(Kokkos::NoInit("sna:ulist"),1,1,1);
    zlist = t_sna_3c_ll(Kokkos::NoInit("sna:zlist"),1,1,1);
    zlist_pack = t_sna_4c_ll(Kokkos::NoInit("sna:zlist_pack"),32,idxz_max,ndoubles,(natom+32-1)/32);
    blist = t_sna_3d_ll(Kokkos::NoInit("sna:blist"),idxb_max,ntriples,natom);
    blist_pack = t_sna_4d_ll(Kokkos::NoInit("sna:blist_pack"),32,idxb_max,ntriples,(natom+32-1)/32);
    ylist = t_sna_3c_ll(Kokkos::NoInit("sna:ylist"),idxu_half_max,nelements,natom);
    ylist_pack_re = t_sna_4d_ll(Kokkos::NoInit("sna:ylist_pack_re"),32,idxu_half_max,nelements,(natom+32-1)/32);
    ylist_pack_im = t_sna_4d_ll(Kokkos::NoInit("sna:ylist_pack_im"),32,idxu_half_max,nelements,(natom+32-1)/32);
    dulist = t_sna_4c3_ll(Kokkos::NoInit("sna:dulist"),1,1,1);
  } else {
#endif
    rij = t_sna_3d(Kokkos::NoInit("sna:rij"),natom,nmax,3);
    wj = t_sna_2d(Kokkos::NoInit("sna:wj"),natom,nmax);
    rcutij = t_sna_2d(Kokkos::NoInit("sna:rcutij"),natom,nmax);
    ulisttot = t_sna_3c_ll(Kokkos::NoInit("sna:ulisttot"),idxu_half_max,nelements,natom);
    ulisttot_full = t_sna_3c_ll(Kokkos::NoInit("sna:ulisttot_full"),idxu_max,nelements,natom);
    ulisttot_re = t_sna_3d_ll(Kokkos::NoInit("sna:ulisttot_re"),1,1,1);
    ulisttot_im = t_sna_3d_ll(Kokkos::NoInit("sna:ulisttot_im"),1,1,1);
    ulisttot_pack = t_sna_4c_ll(Kokkos::NoInit("sna:ulisttot_pack"),1,1,1,1);
    ulist = t_sna_3c_ll(Kokkos::NoInit("sna:ulist"),idxu_cache_max,natom,nmax);
    zlist = t_sna_3c_ll(Kokkos::NoInit("sna:zlist"),idxz_max,ndoubles,natom);
    zlist_pack = t_sna_4c_ll(Kokkos::NoInit("sna:zlist_pack"),1,1,1,1);
    blist = t_sna_3d_ll(Kokkos::NoInit("sna:blist"),idxb_max,ntriples,natom);
    blist_pack = t_sna_4d_ll(Kokkos::NoInit("sna:blist_pack"),1,1,1,1);
    ylist = t_sna_3c_ll(Kokkos::NoInit("sna:ylist"),idxu_half_max,nelements,natom);
    ylist_pack_re = t_sna_4d_ll(Kokkos::NoInit("sna:ylist_pack_re"),1,1,1,1);
    ylist_pack_im = t_sna_4d_ll(Kokkos::NoInit("sna:ylist_pack_im"),1,1,1,1);
    dulist = t_sna_4c3_ll(Kokkos::NoInit("sna:dulist"),idxu_cache_max,natom,nmax);

#ifdef LMP_KOKKOS_GPU
  }
#endif
}

/* ----------------------------------------------------------------------
 * GPU routines
 * ----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------
   Precompute the Cayley-Klein parameters and the derivatives thereof.
   This routine better exploits parallelism than the GPU ComputeUi and
   ComputeFusedDeidrj, which are one warp per atom-neighbor pair.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_cayley_klein(const int& iatom, const int& jnbor, const double& x, const double& y,
                                         const double& z, const double& rcut, const double& wj_local)
{
  //const double x = rij(iatom,jnbor,0);
  //const double y = rij(iatom,jnbor,1);
  //const double z = rij(iatom,jnbor,2);
  const double rsq = x * x + y * y + z * z;
  const double r = sqrt(rsq);
  //const double rcut = rcutij(iatom, jnbor);
  const double rscale0 = rfac0 * MY_PI / (rcut - rmin0);
  const double theta0 = (r - rmin0) * rscale0;
  double sn, cs;
  sincos(theta0, &sn, &cs);
  const double z0 = r * cs / sn;
  const double dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  //const double wj_local = wj(iatom, jnbor);
  double sfac, dsfac;
  compute_s_dsfac(r, rcut, sfac, dsfac);
  sfac *= wj_local;
  dsfac *= wj_local;

  const double rinv = 1.0 / r;
  const double ux = x * rinv;
  const double uy = y * rinv;
  const double uz = z * rinv;

  const double r0inv = 1.0 / sqrt(r * r + z0 * z0);

  const SNAcomplex a = { z0 * r0inv, -z * r0inv };
  const SNAcomplex b = { r0inv * y, -r0inv * x };

  const double dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  const double dr0invx = dr0invdr * ux;
  const double dr0invy = dr0invdr * uy;
  const double dr0invz = dr0invdr * uz;

  const double dz0x = dz0dr * ux;
  const double dz0y = dz0dr * uy;
  const double dz0z = dz0dr * uz;

  const SNAcomplex dax = { dz0x * r0inv + z0 * dr0invx, -z * dr0invx };
  const SNAcomplex day = { dz0y * r0inv + z0 * dr0invy, -z * dr0invy };
  const SNAcomplex daz = { dz0z * r0inv + z0 * dr0invz, -z * dr0invz - r0inv };

  const SNAcomplex dbx = { y * dr0invx, -x * dr0invx - r0inv };
  const SNAcomplex dby = { y * dr0invy + r0inv, -x * dr0invy };
  const SNAcomplex dbz = { y * dr0invz, -x * dr0invz };

  const double dsfacux = dsfac * ux;
  const double dsfacuy = dsfac * uy;
  const double dsfacuz = dsfac * uz;

  CayleyKleinPack ckp{};
  ckp.a = a;
  ckp.b = b;
  ckp.da[0] = dax;
  ckp.db[0] = dbx;
  ckp.da[1] = day;
  ckp.db[1] = dby;
  ckp.da[2] = daz;
  ckp.db[2] = dbz;
  ckp.sfac = sfac;
  ckp.dsfacu[0] = dsfacux;
  ckp.dsfacu[1] = dsfacuy;
  ckp.dsfacu[2] = dsfacuz;

  // Yes, this breaks the standard mantra of using SoA
  // instead of AoS, but it's net fine because of the
  // one warp per atom/neighbor pair for the recursive
  // polynomials. There's good L1 reuse, anyway.
  cayleyklein(iatom, jnbor) = ckp;
}

/* ----------------------------------------------------------------------
   Initialize ulisttot with self-energy terms.
   Ulisttot uses a "half" data layout which takes
   advantage of the symmetry of the Wigner U matrices.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::pre_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int& iatom, const int& ielem)
{
  for (int jelem = 0; jelem < nelements; jelem++) {
    for (int j = 0; j <= twojmax; j++) {
      const int jju_half = idxu_half_block(j);

      // Only diagonal elements get initialized
      // Top half only: gets symmetrized by TransformUi
      // for (int m = 0; m < (j+1)*(j/2+1); m++)
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, (j+1)*(j/2+1)),
        [&] (const int m) {

        const int jjup = jju_half + m;

        // if m is on the "diagonal", initialize it with the self energy.
        // Otherwise zero it out
        double re_part = 0.;
        if (m % (j+2) == 0 && (!chem_flag || ielem == jelem || wselfall_flag)) { re_part = wself; }

        ulisttot_re(jjup, jelem, iatom) = re_part;
        ulisttot_im(jjup, jelem, iatom) = 0.;
      });
    }
  }

}

/* ----------------------------------------------------------------------
   compute Ui by computing Wigner U-functions for one neighbor and
   accumulating to the total. GPU only.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_ui(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom, const int jnbor)
{

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  // get shared memory offset
  const int max_m_tile = (twojmax+1)*(twojmax/2+1);
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * max_m_tile;

  // shared memory double buffer
  // Each level stores a (j+1) x ((j+3)/2) tile in preparation
  // for the computation on the next level. This is the "cached"
  // data layout, which also gets used on the CPU.
  SNAcomplex* buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;
  SNAcomplex* buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;

  // Each warp is accessing the same pack: take advantage of
  // L1 reuse
  const CayleyKleinPack& ckp = cayleyklein(iatom, jnbor);

  const SNAcomplex a = ckp.a;
  const SNAcomplex b = ckp.b;
  const double sfac = ckp.sfac;

  const int jelem = element(iatom, jnbor);

  // VMK Section 4.8.2

  // All writes go to global memory and shared memory
  // so we can avoid all global memory reads
  Kokkos::single(Kokkos::PerThread(team), [=]() {
    buf1[0] = {1.,0.};
    Kokkos::atomic_add(&(ulisttot_re(0,jelem,iatom)), sfac);
  });

  for (int j = 1; j <= twojmax; j++) {

    const int jju = idxu_half_block[j];

    // fill in left side of matrix layer from previous layer

    // Flatten loop over ma, mb
    // for (int ma = 0; ma <= j; ma++)
    const int n_ma = j+1;
    // for (int mb = 0; 2*mb <= j; mb++)
    const int n_mb = j/2+1;

    const int total_iters = n_ma * n_mb;

    //for (int m = 0; m < total_iters; m++) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, total_iters),
      [&] (const int m) {

      // ma fast, mb slow
      // Equivalent to `int ma = m % n_ma; int mb = m / n_ma;` IF everything's positive.
      const int mb = m / n_ma;
      const int ma = m - mb * n_ma;

      // index into global memory array
      const int jju_index = jju+m;

      // index into shared memory buffer for this level
      const int jju_shared_idx = m;

      // index into shared memory buffer for next level
      const int jjup_shared_idx = jju_shared_idx - mb;

      SNAcomplex u_accum = {0., 0.};

      // VMK recursion relation: grab contribution which is multiplied by b*
      const double rootpq2 = -rootpqarray(ma, j - mb);
      const SNAcomplex u_up2 = rootpq2*((ma > 0)?buf1[jjup_shared_idx-1]:SNAcomplex(0.,0.));

      // u_accum += conj(b) * u_up2
      caconjxpy(b, u_up2, u_accum);

      // VMK recursion relation: grab contribution which is multiplied by a*
      const double rootpq1 = rootpqarray(j - ma, j - mb);
      const SNAcomplex u_up1 = rootpq1*((ma < j)?buf1[jjup_shared_idx]:SNAcomplex(0.,0.));

      // u_accum += conj(a) * u_up1
      caconjxpy(a, u_up1, u_accum);

      // back up into shared memory for next iter
      buf2[jju_shared_idx] = u_accum;

      Kokkos::atomic_add(&(ulisttot_re(jju_index,jelem,iatom)), sfac * u_accum.re);
      Kokkos::atomic_add(&(ulisttot_im(jju_index,jelem,iatom)), sfac * u_accum.im);

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j,mb-j] = (-1)^(ma-mb)*Conj([u[ma,mb))
      // if j is even (-> physical j integer), last element maps to self, skip

      const int sign_factor = (((ma+mb)%2==0)?1:-1);
      const int jju_shared_flip = (j+1-mb)*(j+1)-(ma+1);

      if (sign_factor == 1) {
        u_accum.im = -u_accum.im;
      } else {
        u_accum.re = -u_accum.re;
      }

      if (j%2==1 && mb+1==n_mb) {
        buf2[jju_shared_flip] = u_accum;
      }

      // symmetric part of ulisttot is generated in TransformUi

    });
    // In CUDA backend,
    // ThreadVectorRange has a __syncwarp (appropriately masked for
    // vector lengths < 32) implict at the end

    // swap double buffers
    auto tmp = buf1; buf1 = buf2; buf2 = tmp;


  }
}


/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui,
   AoSoA data layout to take advantage of coalescing, avoiding warp
   divergence. GPU version
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_zi(const int& iatom_mod, const int& jjz, const int& iatom_div)
{

  const int j1 = idxz(jjz, 0);
  const int j2 = idxz(jjz, 1);
  const int j = idxz(jjz, 2);
  const int ma1min = idxz(jjz, 3);
  const int ma2max = idxz(jjz, 4);
  const int mb1min = idxz(jjz, 5);
  const int mb2max = idxz(jjz, 6);
  const int na = idxz(jjz, 7);
  const int nb = idxz(jjz, 8);

  const double* cgblock = cglist.data() + idxcg_block(j1, j2, j);

  int idouble = 0;

  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {
      double ztmp_r = 0.;
      double ztmp_i = 0.;

      int jju1 = idxu_block[j1] + (j1+1)*mb1min;
      int jju2 = idxu_block[j2] + (j2+1)*mb2max;
      int icgb = mb1min*(j2+1) + mb2max;

      #ifdef LMP_KK_DEVICE_COMPILE
      #pragma unroll
      #endif
      for(int ib = 0; ib < nb; ib++) {

        int ma1 = ma1min;
        int ma2 = ma2max;
        int icga = ma1min*(j2+1) + ma2max;

        #ifdef LMP_KK_DEVICE_COMPILE
        #pragma unroll
        #endif
        for(int ia = 0; ia < na; ia++) {
          const SNAcomplex utot1 = ulisttot_pack(iatom_mod, jju1+ma1, elem1, iatom_div);
          const SNAcomplex utot2 = ulisttot_pack(iatom_mod, jju2+ma2, elem2, iatom_div);
          const auto cgcoeff_a = cgblock[icga];
          const auto cgcoeff_b = cgblock[icgb];
          ztmp_r += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.re - utot1.im * utot2.im);
          ztmp_i += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.im + utot1.im * utot2.re);
          ma1++;
          ma2--;
          icga += j2;
        } // end loop over ia

        jju1 += j1 + 1;
        jju2 -= j2 + 1;
        icgb += j2;
      } // end loop over ib

      if (bnorm_flag) {
        ztmp_r /= (j + 1);
        ztmp_i /= (j + 1);
      }

      zlist_pack(iatom_mod,jjz,idouble,iatom_div) = { ztmp_r, ztmp_i };

      idouble++;
    }
  }
}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
   AoSoA data layout to take advantage of coalescing, avoiding warp
   divergence.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_bi(const int& iatom_mod, const int& jjb, const int& iatom_div)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  const int j1 = idxb(jjb,0);
  const int j2 = idxb(jjb,1);
  const int j = idxb(jjb,2);

  const int jjz = idxz_block(j1,j2,j);
  const int jju = idxu_block[j];

  int itriple = 0;
  int idouble = 0;
  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {
      for (int elem3 = 0; elem3 < nelements; elem3++) {

        double sumzu = 0.0;
        double sumzu_temp = 0.0;

        for(int mb = 0; 2*mb < j; mb++) {
          for(int ma = 0; ma <= j; ma++) {
            const int jju_index = jju+mb*(j+1)+ma;
            const int jjz_index = jjz+mb*(j+1)+ma;
            if (2*mb == j) return; // I think we can remove this?
            const auto utot = ulisttot_pack(iatom_mod, jju_index, elem3, iatom_div);
            const auto zloc = zlist_pack(iatom_mod, jjz_index, idouble, iatom_div);
            sumzu_temp += utot.re * zloc.re + utot.im * zloc.im;
          }
        }
        sumzu += sumzu_temp;

        // For j even, special treatment for middle column
        if (j%2 == 0) {
          sumzu_temp = 0.;

          const int mb = j/2;
          for(int ma = 0; ma < mb; ma++) {
            const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
            const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;

            const auto utot = ulisttot_pack(iatom_mod, jju_index, elem3, iatom_div);
            const auto zloc = zlist_pack(iatom_mod, jjz_index, idouble, iatom_div);
            sumzu_temp += utot.re * zloc.re + utot.im * zloc.im;

          }
          sumzu += sumzu_temp;

          const int ma = mb;
          const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
          const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;

          const auto utot = ulisttot_pack(iatom_mod, jju_index, elem3, iatom_div);
          const auto zloc = zlist_pack(iatom_mod, jjz_index, idouble, iatom_div);
          sumzu += 0.5 * (utot.re * zloc.re + utot.im * zloc.im);
        } // end if jeven

        sumzu *= 2.0;
        if (bzero_flag) {
          if (!wselfall_flag) {
            if (elem1 == elem2 && elem1 == elem3) {
              sumzu -= bzero[j];
            }
          } else {
            sumzu -= bzero[j];
          }
        }
        blist_pack(iatom_mod, jjb, itriple, iatom_div) = sumzu;
            //} // end loop over j
          //} // end loop over j1, j2
        itriple++;
      } // end loop over elem3
      idouble++;
    } // end loop over elem2
  } // end loop over elem1
}


/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices.
   AoSoA data layout to take advantage of coalescing, avoiding warp
   divergence. GPU version.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_yi(int iatom_mod, int jjz, int iatom_div,
 const Kokkos::View<F_FLOAT***, Kokkos::LayoutLeft, DeviceType> &beta_pack)
{
  double betaj;

  const int j1 = idxz(jjz, 0);
  const int j2 = idxz(jjz, 1);
  const int j = idxz(jjz, 2);
  const int ma1min = idxz(jjz, 3);
  const int ma2max = idxz(jjz, 4);
  const int mb1min = idxz(jjz, 5);
  const int mb2max = idxz(jjz, 6);
  const int na = idxz(jjz, 7);
  const int nb = idxz(jjz, 8);
  const int jju_half = idxz(jjz, 9);

  const double *cgblock = cglist.data() + idxcg_block(j1,j2,j);
  //int mb = (2 * (mb1min+mb2max) - j1 - j2 + j) / 2;
  //int ma = (2 * (ma1min+ma2max) - j1 - j2 + j) / 2;

  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {

      double ztmp_r = 0.0;
      double ztmp_i = 0.0;

      int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
      int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
      int icgb = mb1min * (j2 + 1) + mb2max;

      #ifdef LMP_KK_DEVICE_COMPILE
      #pragma unroll
      #endif
      for (int ib = 0; ib < nb; ib++) {

        int ma1 = ma1min;
        int ma2 = ma2max;
        int icga = ma1min*(j2+1) + ma2max;

        #ifdef LMP_KK_DEVICE_COMPILE
        #pragma unroll
        #endif
        for (int ia = 0; ia < na; ia++) {
          const SNAcomplex utot1 = ulisttot_pack(iatom_mod,jju1+ma1,elem1,iatom_div);
          const SNAcomplex utot2 = ulisttot_pack(iatom_mod,jju2+ma2,elem2,iatom_div);
          const auto cgcoeff_a = cgblock[icga];
          const auto cgcoeff_b = cgblock[icgb];
          ztmp_r += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.re - utot1.im * utot2.im);
          ztmp_i += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.im + utot1.im * utot2.re);
          ma1++;
          ma2--;
          icga += j2;
        } // end loop over ia

        jju1 += j1 + 1;
        jju2 -= j2 + 1;
        icgb += j2;
      } // end loop over ib

      if (bnorm_flag) {
        ztmp_r /= (j + 1);
        ztmp_i /= (j + 1);
      }

      // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
      // find right y_list[jju] and beta(iatom,jjb) entries
      // multiply and divide by j+1 factors
      // account for multiplicity of 1, 2, or 3

      // pick out right beta value
      for (int elem3 = 0; elem3 < nelements; elem3++) {
        if (j >= j1) {
          const int jjb = idxb_block(j1, j2, j);
          const auto itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
          if (j1 == j) {
            if (j2 == j) betaj = 3 * beta_pack(iatom_mod, itriple, iatom_div);
            else betaj = 2 * beta_pack(iatom_mod, itriple, iatom_div);
          } else betaj = beta_pack(iatom_mod, itriple, iatom_div);
        } else if (j >= j2) {
          const int jjb = idxb_block(j, j2, j1);
          const auto itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
          if (j2 == j) betaj = 2 * beta_pack(iatom_mod, itriple, iatom_div);
          else betaj = beta_pack(iatom_mod, itriple, iatom_div);
        } else {
          const int jjb = idxb_block(j2, j, j1);
          const auto itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
          betaj = beta_pack(iatom_mod, itriple, iatom_div);
        }

        if (!bnorm_flag && j1 > j)
          betaj *= (j1 + 1) / (j + 1.0);


        Kokkos::atomic_add(&(ylist_pack_re(iatom_mod, jju_half, elem3, iatom_div)), betaj*ztmp_r);
        Kokkos::atomic_add(&(ylist_pack_im(iatom_mod, jju_half, elem3, iatom_div)), betaj*ztmp_i);
      } // end loop over elem3
    } // end loop over elem2
  } // end loop over elem1
}

/* ----------------------------------------------------------------------
   Fused calculation of the derivative of Ui w.r.t. atom j
   and accumulation into dEidRj. GPU only.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_fused_deidrj(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom, const int jnbor)
{
  // get shared memory offset
  const int max_m_tile = (twojmax+1)*(twojmax/2+1);
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * max_m_tile;
  const int jelem = element(iatom, jnbor);

  // See notes on data layouts for shared memory caching
  // in `compute_ui`.

  // double buffer for ulist
  SNAcomplex* ulist_buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;
  SNAcomplex* ulist_buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;

  // double buffer for dulist
  SNAcomplex* dulist_buf1 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;
  SNAcomplex* dulist_buf2 = (SNAcomplex*)team.team_shmem( ).get_shmem(team.team_size()*max_m_tile*sizeof(SNAcomplex), 0) + scratch_shift;

  const CayleyKleinPack& ckp = cayleyklein(iatom, jnbor);

  const SNAcomplex a = ckp.a;
  const SNAcomplex b = ckp.b;
  const SNAcomplex da = ckp.da[dir];
  const SNAcomplex db = ckp.db[dir];
  const double sfac = ckp.sfac;
  const double dsfacu = ckp.dsfacu[dir]; // dsfac * u

  // Accumulate the full contribution to dedr on the fly
  const SNAcomplex y_local = ylist(0, jelem, iatom);

  // Symmetry factor of 0.5 b/c 0 element is on diagonal for even j==0
  double dedr_full_sum = 0.5 * dsfacu * y_local.re;

  // single has a warp barrier at the end
  Kokkos::single(Kokkos::PerThread(team), [=]() {

    ulist_buf1[0] = {1., 0.};
    dulist_buf1[0] = {0., 0.};
  });

  for (int j = 1; j <= twojmax; j++) {
    int jju_half = idxu_half_block[j];

    // flatten the loop over ma,mb

    // for (int ma = 0; ma <= j; ma++)
    const int n_ma = j+1;
    // for (int mb = 0; 2*mb <= j; mb++)
    const int n_mb = j/2+1;

    const int total_iters = n_ma * n_mb;

    double dedr_sum = 0.; // j-local sum

    //for (int m = 0; m < total_iters; m++) {
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, total_iters),
      [&] (const int m, double& sum_tmp) {

      // ma fast, mb slow
      const int mb = m / n_ma;
      const int ma = m - mb * n_ma;

      const int jju_half_index = jju_half+m;

      // Load y_local, apply the symmetry scaling factor
      // The "secret" of the shared memory optimization is it eliminates
      // all global memory reads to duidrj in lieu of caching values in
      // shared memory and otherwise always writing, making the kernel
      // ultimately compute bound. We take advantage of that by adding
      // some reads back in.
      auto y_local = ylist(jju_half_index, jelem, iatom);
      if (j % 2 == 0 && 2*mb == j) {
        if (ma == mb) { y_local = 0.5*y_local; }
        else if (ma > mb) { y_local = { 0., 0. }; } // can probably avoid this outright
        // else the ma < mb gets "double counted", cancelling the 0.5.
      }

      // index into shared memory
      const int jju_shared_idx = m;
      const int jjup_shared_idx = jju_shared_idx - mb;

      // Need to compute and accumulate both u and du (mayhaps, we could probably
      // balance some read and compute by reading u each time).
      SNAcomplex u_accum = { 0., 0. };
      SNAcomplex du_accum = { 0., 0. };

      const double rootpq2 = -rootpqarray(ma, j - mb);
      const SNAcomplex u_up2 = rootpq2*((ma > 0) ? ulist_buf1[jjup_shared_idx-1]:SNAcomplex(0.,0.));

      // u_accum += conj(b) * u_up2
      caconjxpy(b, u_up2, u_accum);

      const double rootpq1 = rootpqarray(j - ma, j - mb);
      const SNAcomplex u_up1 = rootpq1*((ma < j) ? ulist_buf1[jjup_shared_idx]:SNAcomplex(0.,0.));

      // u_accum += conj(a) * u_up1
      caconjxpy(a, u_up1, u_accum);

      // Next, spin up du_accum
      const SNAcomplex du_up1 = rootpq1*((ma < j) ? dulist_buf1[jjup_shared_idx] : SNAcomplex(0.,0.));

      // du_accum += conj(da) * u_up1 + conj(a) * du_up1
      caconjxpy(da, u_up1, du_accum);
      caconjxpy(a, du_up1, du_accum);

      const SNAcomplex du_up2 = rootpq2*((ma > 0) ? dulist_buf1[jjup_shared_idx-1] : SNAcomplex(0.,0.));

      // du_accum += conj(db) * u_up2 + conj(b) * du_up2
      caconjxpy(db, u_up2, du_accum);
      caconjxpy(b, du_up2, du_accum);

      // Cache u_accum, du_accum to scratch memory.
      ulist_buf2[jju_shared_idx] = u_accum;
      dulist_buf2[jju_shared_idx] = du_accum;

      // Directly accumulate deidrj into sum_tmp
      const SNAcomplex du_prod = (dsfacu * u_accum) + (sfac * du_accum);
      sum_tmp += du_prod.re * y_local.re + du_prod.im * y_local.im;

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
      if (j%2==1 && mb+1==n_mb) {
        int sign_factor = (((ma+mb)%2==0)?1:-1);

        const int jju_shared_flip = (j+1-mb)*(j+1)-(ma+1);

        if (sign_factor == 1) {
          u_accum.im = -u_accum.im;
          du_accum.im = -du_accum.im;
        } else {
          u_accum.re = -u_accum.re;
          du_accum.re = -du_accum.re;
        }

        // We don't need the second half of the tile for the deidrj accumulation.
        // That's taken care of by the symmetry factor above.
        // We do need it for ortho polynomial generation, though
        ulist_buf2[jju_shared_flip] = u_accum;
        dulist_buf2[jju_shared_flip] = du_accum;
      }

    }, dedr_sum);

    // swap buffers
    auto tmp = ulist_buf1; ulist_buf1 = ulist_buf2; ulist_buf2 = tmp;
    tmp = dulist_buf1; dulist_buf1 = dulist_buf2; dulist_buf2 = tmp;

    // Accumulate dedr. This "should" be in a single, but
    // a Kokkos::single call implies a warp sync, and we may
    // as well avoid that. This does no harm as long as the
    // final assignment is in a single block.
    //Kokkos::single(Kokkos::PerThread(team), [=]() {
    dedr_full_sum += dedr_sum;
    //});
  }

  // Store the accumulated dedr.
  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr(iatom,jnbor,dir) = dedr_full_sum*2.0;
  });
}



/* ----------------------------------------------------------------------
 * CPU routines
 * ----------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
   Ulisttot uses a "half" data layout which takes
   advantage of the symmetry of the Wigner U matrices.
 * ------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::pre_ui_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int& iatom, const int& ielem)
{
  for (int jelem = 0; jelem < nelements; jelem++) {
    for (int j = 0; j <= twojmax; j++) {
      const int jju = idxu_half_block(j);

      // Only diagonal elements get initialized
      // for (int m = 0; m < (j+1)*(j/2+1); m++)
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, (j+1)*(j/2+1)),
        [&] (const int m) {

        const int jjup = jju + m;

        // if m is on the "diagonal", initialize it with the self energy.
        // Otherwise zero it out
        SNAcomplex init = {0., 0.};
        if (m % (j+2) == 0 && (!chem_flag || ielem == jelem || wselfall_flag)) { init = {wself, 0.0}; } //need to map iatom to element

        ulisttot(jjup, jelem, iatom) = init;
      });
    }
  }

}


/* ----------------------------------------------------------------------
   compute Ui by summing over bispectrum components. CPU only.
   See comments above compute_uarray_cpu and add_uarraytot for
   data layout comments.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_ui_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);

  theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  //    theta0 = (r - rmin0) * rscale0;
  z0 = r / tan(theta0);

  compute_uarray_cpu(team, iatom, jnbor, x, y, z, z0, r);
  add_uarraytot(team, iatom, jnbor, r, wj(iatom,jnbor), rcutij(iatom,jnbor), element(iatom, jnbor));

}
/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui, CPU version
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_zi_cpu(const int& iter)
{
  const int iatom = iter / idxz_max;
  const int jjz = iter % idxz_max;

  const int j1 = idxz(jjz, 0);
  const int j2 = idxz(jjz, 1);
  const int j = idxz(jjz, 2);
  const int ma1min = idxz(jjz, 3);
  const int ma2max = idxz(jjz, 4);
  const int mb1min = idxz(jjz, 5);
  const int mb2max = idxz(jjz, 6);
  const int na = idxz(jjz, 7);
  const int nb = idxz(jjz, 8);

  const double *cgblock = cglist.data() + idxcg_block(j1,j2,j);

  int idouble = 0;

  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {
      zlist(jjz, idouble, iatom).re = 0.0;
      zlist(jjz, idouble, iatom).im = 0.0;

      int jju1 = idxu_block[j1] + (j1+1)*mb1min;
      int jju2 = idxu_block[j2] + (j2+1)*mb2max;
      int icgb = mb1min*(j2+1) + mb2max;
      for(int ib = 0; ib < nb; ib++) {

        double suma1_r = 0.0;
        double suma1_i = 0.0;

        int ma1 = ma1min;
        int ma2 = ma2max;
        int icga = ma1min * (j2 + 1) + ma2max;
        for(int ia = 0; ia < na; ia++) {
          suma1_r += cgblock[icga] * (ulisttot_full(jju1+ma1, elem1, iatom).re * ulisttot_full(jju2+ma2, elem2, iatom).re -
                                      ulisttot_full(jju1+ma1, elem1, iatom).im * ulisttot_full(jju2+ma2, elem2, iatom).im);
          suma1_i += cgblock[icga] * (ulisttot_full(jju1+ma1, elem1, iatom).re * ulisttot_full(jju2+ma2, elem2, iatom).im +
                                      ulisttot_full(jju1+ma1, elem1, iatom).im * ulisttot_full(jju2+ma2, elem2, iatom).re);
          ma1++;
          ma2--;
          icga += j2;
        } // end loop over ia

        zlist(jjz, idouble, iatom).re += cgblock[icgb] * suma1_r;
        zlist(jjz, idouble, iatom).im += cgblock[icgb] * suma1_i;

        jju1 += j1 + 1;
        jju2 -= j2 + 1;
        icgb += j2;
      } // end loop over ib

      if (bnorm_flag) {
        zlist(jjz, idouble, iatom).re /= (j+1);
        zlist(jjz, idouble, iatom).im /= (j+1);
      }
      idouble++;
    } // end loop over elem2
  } // end loop over elem1
}


/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi, CPU version
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_bi_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  int itriple = 0;
  int idouble = 0;
  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {
      auto jalloy = idouble; // must be non-const to work around gcc compiler bug
      for (int elem3 = 0; elem3 < nelements; elem3++) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,idxb_max),
          [&] (const int& jjb) {
        //for(int jjb = 0; jjb < idxb_max; jjb++) {
          const int j1 = idxb(jjb, 0);
          const int j2 = idxb(jjb, 1);
          const int j = idxb(jjb, 2);

          int jjz = idxz_block(j1, j2, j);
          int jju = idxu_block[j];
          double sumzu = 0.0;
          double sumzu_temp = 0.0;
          const int bound = (j+2)/2;
          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,(j+1)*bound),
              [&] (const int mbma, double& sum) {
              //for(int mb = 0; 2*mb < j; mb++)
                //for(int ma = 0; ma <= j; ma++) {
              const int ma = mbma % (j + 1);
              const int mb = mbma / (j + 1);
              const int jju_index = jju + mb * (j + 1) + ma;
              const int jjz_index = jjz + mb * (j + 1) + ma;
              if (2*mb == j) return;
              sum +=
                ulisttot_full(jju_index, elem3, iatom).re * zlist(jjz_index, jalloy, iatom).re +
                ulisttot_full(jju_index, elem3, iatom).im * zlist(jjz_index, jalloy, iatom).im;
            },sumzu_temp); // end loop over ma, mb
            sumzu += sumzu_temp;

          // For j even, special treatment for middle column

          if (j%2 == 0) {
            const int mb = j/2;
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, mb),
                [&] (const int ma, double& sum) {
            //for(int ma = 0; ma < mb; ma++) {
              const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
              const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
              sum +=
                ulisttot_full(jju_index, elem3, iatom).re * zlist(jjz_index, jalloy, iatom).re +
                ulisttot_full(jju_index, elem3, iatom).im * zlist(jjz_index, jalloy, iatom).im;
            },sumzu_temp); // end loop over ma
            sumzu += sumzu_temp;

            const int ma = mb;
            const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
            const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;
            sumzu += 0.5*
              (ulisttot_full(jju_index, elem3, iatom).re * zlist(jjz_index, jalloy, iatom).re +
               ulisttot_full(jju_index, elem3, iatom).im * zlist(jjz_index, jalloy, iatom).im);
          } // end if jeven

          Kokkos::single(Kokkos::PerThread(team), [&] () {
            sumzu *= 2.0;

            // apply bzero shift

            if (bzero_flag){
              if (!wselfall_flag) {
                if (elem1 == elem2 && elem1 == elem3) {
                  sumzu -= bzero[j];
                }
              } else {
                sumzu -= bzero[j];
              }
            }

            blist(jjb, itriple, iatom) = sumzu;
          });
        });
          //} // end loop over j
        //} // end loop over j1, j2
        itriple++;
      }
      idouble++;
    } // end loop over elem2
  } // end loop over elem1

}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices,
   CPU version
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_yi_cpu(int iter,
 const Kokkos::View<F_FLOAT**, DeviceType> &beta)
{
  double betaj;
  const int iatom = iter / idxz_max;
  const int jjz = iter % idxz_max;

  const int j1 = idxz(jjz, 0);
  const int j2 = idxz(jjz, 1);
  const int j = idxz(jjz, 2);
  const int ma1min = idxz(jjz, 3);
  const int ma2max = idxz(jjz, 4);
  const int mb1min = idxz(jjz, 5);
  const int mb2max = idxz(jjz, 6);
  const int na = idxz(jjz, 7);
  const int nb = idxz(jjz, 8);
  const int jju_half = idxz(jjz, 9);

  const double *cgblock = cglist.data() + idxcg_block(j1,j2,j);
  //int mb = (2 * (mb1min+mb2max) - j1 - j2 + j) / 2;
  //int ma = (2 * (ma1min+ma2max) - j1 - j2 + j) / 2;

  for (int elem1 = 0; elem1 < nelements; elem1++) {
    for (int elem2 = 0; elem2 < nelements; elem2++) {

      double ztmp_r = 0.0;
      double ztmp_i = 0.0;

      int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
      int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
      int icgb = mb1min * (j2 +1 ) + mb2max;

      for (int ib = 0; ib < nb; ib++) {

        double suma1_r = 0.0;
        double suma1_i = 0.0;

        int ma1 = ma1min;
        int ma2 = ma2max;
        int icga = ma1min*(j2+1) + ma2max;

        for (int ia = 0; ia < na; ia++) {
          suma1_r += cgblock[icga] * (ulisttot_full(jju1+ma1, elem1, iatom).re * ulisttot_full(jju2+ma2, elem2, iatom).re -
                                      ulisttot_full(jju1+ma1, elem1, iatom).im * ulisttot_full(jju2+ma2, elem2, iatom).im);
          suma1_i += cgblock[icga] * (ulisttot_full(jju1+ma1, elem1, iatom).re * ulisttot_full(jju2+ma2, elem2, iatom).im +
                                      ulisttot_full(jju1+ma1, elem1, iatom).im * ulisttot_full(jju2+ma2, elem2, iatom).re);
          ma1++;
          ma2--;
          icga += j2;
        } // end loop over ia

        ztmp_r += cgblock[icgb] * suma1_r;
        ztmp_i += cgblock[icgb] * suma1_i;
        jju1 += j1 + 1;
        jju2 -= j2 + 1;
        icgb += j2;
      } // end loop over ib

      if (bnorm_flag) {
        ztmp_i /= j + 1;
        ztmp_r /= j + 1;
      }

      // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
      // find right y_list[jju] and beta(iatom,jjb) entries
      // multiply and divide by j+1 factors
      // account for multiplicity of 1, 2, or 3

      // pick out right beta value
      for (int elem3 = 0; elem3 < nelements; elem3++) {

        if (j >= j1) {
          const int jjb = idxb_block(j1, j2, j);
          const auto itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
          if (j1 == j) {
            if (j2 == j) betaj = 3 * beta(itriple, iatom);
            else betaj = 2 * beta(itriple, iatom);
          } else betaj = beta(itriple, iatom);
        } else if (j >= j2) {
          const int jjb = idxb_block(j, j2, j1);
          const auto itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
          if (j2 == j) betaj = 2 * beta(itriple, iatom);
          else betaj = beta(itriple, iatom);
        } else {
          const int jjb = idxb_block(j2, j, j1);
          const auto itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
          betaj = beta(itriple, iatom);
        }

        if (!bnorm_flag && j1 > j)
          betaj *= (j1 + 1) / (j + 1.0);

        Kokkos::atomic_add(&(ylist(jju_half, elem3, iatom).re), betaj*ztmp_r);
        Kokkos::atomic_add(&(ylist(jju_half, elem3, iatom).im), betaj*ztmp_i);
      } // end loop over elem3
    } // end loop over elem2
  } // end loop over elem1
}


/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
   see comments above compute_duarray_cpu for comments on the
   data layout
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij(iatom,jnbor,0);
  y = rij(iatom,jnbor,1);
  z = rij(iatom,jnbor,2);
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  compute_duarray_cpu(team, iatom, jnbor, x, y, z, z0, r, dz0dr, wj(iatom,jnbor), rcutij(iatom,jnbor));
}


/* ----------------------------------------------------------------------
   compute dEidRj, CPU path only.
   dulist takes advantage of a `cached` data layout, similar to the
   shared memory layout for the GPU routines, which is efficient for
   compressing the calculation in compute_duarray_cpu. That said,
   dulist only uses the "half" data layout part of that structure.
------------------------------------------------------------------------- */


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_deidrj_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor)
{
  t_scalar3<double> final_sum;
  const int jelem = element(iatom, jnbor);

  //for(int j = 0; j <= twojmax; j++) {
  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,twojmax+1),
      [&] (const int& j, t_scalar3<double>& sum_tmp) {
    int jju_half = idxu_half_block[j];
    int jju_cache = idxu_cache_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
        sum_tmp.x += dulist(jju_cache,iatom,jnbor,0).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,0).im * ylist(jju_half,jelem,iatom).im;
        sum_tmp.y += dulist(jju_cache,iatom,jnbor,1).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,1).im * ylist(jju_half,jelem,iatom).im;
        sum_tmp.z += dulist(jju_cache,iatom,jnbor,2).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,2).im * ylist(jju_half,jelem,iatom).im;
        jju_half++; jju_cache++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        sum_tmp.x += dulist(jju_cache,iatom,jnbor,0).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,0).im * ylist(jju_half,jelem,iatom).im;
        sum_tmp.y += dulist(jju_cache,iatom,jnbor,1).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,1).im * ylist(jju_half,jelem,iatom).im;
        sum_tmp.z += dulist(jju_cache,iatom,jnbor,2).re * ylist(jju_half,jelem,iatom).re +
                     dulist(jju_cache,iatom,jnbor,2).im * ylist(jju_half,jelem,iatom).im;
        jju_half++; jju_cache++;
      }

      //int ma = mb;
      sum_tmp.x += (dulist(jju_cache,iatom,jnbor,0).re * ylist(jju_half,jelem,iatom).re +
                    dulist(jju_cache,iatom,jnbor,0).im * ylist(jju_half,jelem,iatom).im)*0.5;
      sum_tmp.y += (dulist(jju_cache,iatom,jnbor,1).re * ylist(jju_half,jelem,iatom).re +
                    dulist(jju_cache,iatom,jnbor,1).im * ylist(jju_half,jelem,iatom).im)*0.5;
      sum_tmp.z += (dulist(jju_cache,iatom,jnbor,2).re * ylist(jju_half,jelem,iatom).re +
                    dulist(jju_cache,iatom,jnbor,2).im * ylist(jju_half,jelem,iatom).im)*0.5;
    } // end if jeven

  },final_sum); // end loop over j

  Kokkos::single(Kokkos::PerThread(team), [&] () {
    dedr(iatom,jnbor,0) = final_sum.x*2.0;
    dedr(iatom,jnbor,1) = final_sum.y*2.0;
    dedr(iatom,jnbor,2) = final_sum.z*2.0;
  });

}


/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
   ulist is in a "cached" data layout, which is a compressed layout
   which still keeps the recursive calculation simple. On the other hand
   `ulisttot` uses a "half" data layout, which fully takes advantage
   of the symmetry of the Wigner U matrices.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::add_uarraytot(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                                          double r, double wj, double rcut, int jelem)
{
  const double sfac = compute_sfac(r, rcut) * wj;

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,twojmax+1),
      [&] (const int& j) {

    int jju_half = idxu_half_block[j]; // index into ulisttot
    int jju_cache = idxu_cache_block[j]; // index into ulist

    int count = 0;
    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        Kokkos::atomic_add(&(ulisttot(jju_half+count, jelem, iatom).re), sfac * ulist(jju_cache+count, iatom, jnbor).re);
        Kokkos::atomic_add(&(ulisttot(jju_half+count, jelem, iatom).im), sfac * ulist(jju_cache+count, iatom, jnbor).im);
        count++;
      }
    }
  });
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor.
   `ulisttot` uses a "cached" data layout, matching the amount of
   information stored between layers via scratch memory on the GPU path
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_uarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                         double x, double y, double z,
                         double z0, double r)
{
  double r0inv;
  double a_r, b_r, a_i, b_i;
  double rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2

  ulist(0,iatom,jnbor).re = 1.0;
  ulist(0,iatom,jnbor).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    const int jju = idxu_cache_block[j];
    const int jjup = idxu_cache_block[j-1];

    // fill in left side of matrix layer from previous layer

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      ulist(jju_index,iatom,jnbor).re = 0.0;
      ulist(jju_index,iatom,jnbor).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        ulist(jju_index,iatom,jnbor).re +=
          rootpq *
          (a_r * ulist(jjup_index,iatom,jnbor).re +
           a_i * ulist(jjup_index,iatom,jnbor).im);
        ulist(jju_index,iatom,jnbor).im +=
          rootpq *
          (a_r * ulist(jjup_index,iatom,jnbor).im -
           a_i * ulist(jjup_index,iatom,jnbor).re);

        rootpq = rootpqarray(ma + 1,j - mb);
        ulist(jju_index+1,iatom,jnbor).re =
          -rootpq *
          (b_r * ulist(jjup_index,iatom,jnbor).re +
           b_i * ulist(jjup_index,iatom,jnbor).im);
        ulist(jju_index+1,iatom,jnbor).im =
          -rootpq *
          (b_r * ulist(jjup_index,iatom,jnbor).im -
           b_i * ulist(jjup_index,iatom,jnbor).re);
      }

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j,mb-j] = (-1)^(ma-mb)*Conj([u[ma,mb))

      // Only need to add one symmetrized row for convenience
      // Symmetry gets "unfolded" in accumulating ulisttot
      if (j%2==1 && mb==(j/2)) {
        const int mbpar = (mb)%2==0?1:-1;
        int mapar = mbpar;
        for (int ma = 0; ma <= j; ma++) {
          const int jju_index = jju + mb*(j+1) + ma;
          const int jjup_index = jju + (j+1-mb)*(j+1)-(ma+1);
          if (mapar == 1) {
            ulist(jjup_index,iatom,jnbor).re = ulist(jju_index,iatom,jnbor).re;
            ulist(jjup_index,iatom,jnbor).im = -ulist(jju_index,iatom,jnbor).im;
          } else {
            ulist(jjup_index,iatom,jnbor).re = -ulist(jju_index,iatom,jnbor).re;
            ulist(jjup_index,iatom,jnbor).im = ulist(jju_index,iatom,jnbor).im;
          }
          mapar = -mapar;
        }
      }
    });

  }
}

/* ----------------------------------------------------------------------
   compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray_cpu()
   Uses same cached data layout of ulist
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_duarray_cpu(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, int iatom, int jnbor,
                          double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut)
{
double r0inv;
  double a_r, a_i, b_r, b_i;
  double da_r[3], da_i[3], db_r[3], db_i[3];
  double dz0[3], dr0inv[3], dr0invdr;
  double rootpq;

  double rinv = 1.0 / r;
  double ux = x * rinv;
  double uy = y * rinv;
  double uz = z * rinv;

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

  dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  for (int k = 0; k < 3; k++) {
    da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }

  da_i[2] += -r0inv;

  for (int k = 0; k < 3; k++) {
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }

  db_i[0] += -r0inv;
  db_r[1] += r0inv;

  dulist(0,iatom,jnbor,0).re = 0.0;
  dulist(0,iatom,jnbor,1).re = 0.0;
  dulist(0,iatom,jnbor,2).re = 0.0;
  dulist(0,iatom,jnbor,0).im = 0.0;
  dulist(0,iatom,jnbor,1).im = 0.0;
  dulist(0,iatom,jnbor,2).im = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_cache_block[j];
    int jjup = idxu_cache_block[j-1];
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,(j+2)/2),
        [&] (const int& mb) {
    //for (int mb = 0; 2*mb <= j; mb++) {
      const int jju_index = jju+mb+mb*j;
      dulist(jju_index,iatom,jnbor,0).re = 0.0;
      dulist(jju_index,iatom,jnbor,1).re = 0.0;
      dulist(jju_index,iatom,jnbor,2).re = 0.0;
      dulist(jju_index,iatom,jnbor,0).im = 0.0;
      dulist(jju_index,iatom,jnbor,1).im = 0.0;
      dulist(jju_index,iatom,jnbor,2).im = 0.0;

      for (int ma = 0; ma < j; ma++) {
        const int jju_index = jju+mb+mb*j+ma;
        const int jjup_index = jjup+mb*j+ma;
        rootpq = rootpqarray(j - ma,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(jju_index,iatom,jnbor,k).re +=
            rootpq * (da_r[k] * ulist(jjup_index,iatom,jnbor).re +
                      da_i[k] * ulist(jjup_index,iatom,jnbor).im +
                      a_r * dulist(jjup_index,iatom,jnbor,k).re +
                      a_i * dulist(jjup_index,iatom,jnbor,k).im);
          dulist(jju_index,iatom,jnbor,k).im +=
            rootpq * (da_r[k] * ulist(jjup_index,iatom,jnbor).im -
                      da_i[k] * ulist(jjup_index,iatom,jnbor).re +
                      a_r * dulist(jjup_index,iatom,jnbor,k).im -
                      a_i * dulist(jjup_index,iatom,jnbor,k).re);
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        for (int k = 0; k < 3; k++) {
          dulist(jju_index+1,iatom,jnbor,k).re =
            -rootpq * (db_r[k] * ulist(jjup_index,iatom,jnbor).re +
                       db_i[k] * ulist(jjup_index,iatom,jnbor).im +
                       b_r * dulist(jjup_index,iatom,jnbor,k).re +
                       b_i * dulist(jjup_index,iatom,jnbor,k).im);
          dulist(jju_index+1,iatom,jnbor,k).im =
            -rootpq * (db_r[k] * ulist(jjup_index,iatom,jnbor).im -
                       db_i[k] * ulist(jjup_index,iatom,jnbor).re +
                       b_r * dulist(jjup_index,iatom,jnbor,k).im -
                       b_i * dulist(jjup_index,iatom,jnbor,k).re);
        }
      }

      // Only need to add one symmetrized row for convenience
      // Symmetry gets "unfolded" during the dedr accumulation

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

      if (j%2==1 && mb==(j/2)) {
        const int mbpar = (mb)%2==0?1:-1;
        int mapar = mbpar;
        for (int ma = 0; ma <= j; ma++) {
          const int jju_index = jju+mb*(j+1)+ma;
          const int jjup_index = jju+(mb+2)*(j+1)-(ma+1);
          if (mapar == 1) {
            for (int k = 0; k < 3; k++) {
              dulist(jjup_index,iatom,jnbor,k).re = dulist(jju_index,iatom,jnbor,k).re;
              dulist(jjup_index,iatom,jnbor,k).im = -dulist(jju_index,iatom,jnbor,k).im;
            }
          } else {
            for (int k = 0; k < 3; k++) {
              dulist(jjup_index,iatom,jnbor,k).re = -dulist(jju_index,iatom,jnbor,k).re;
              dulist(jjup_index,iatom,jnbor,k).im = dulist(jju_index,iatom,jnbor,k).im;
            }
          }
          mapar = -mapar;
        }
      }
    });
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;

  // Even though we fill out a full "cached" data layout above,
  // we only need the "half" data for the accumulation into dedr.
  // Thus we skip updating any unnecessary data.
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_cache_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist(jju,iatom,jnbor,0).re = dsfac * ulist(jju,iatom,jnbor).re * ux +
                                  sfac * dulist(jju,iatom,jnbor,0).re;
        dulist(jju,iatom,jnbor,0).im = dsfac * ulist(jju,iatom,jnbor).im * ux +
                                  sfac * dulist(jju,iatom,jnbor,0).im;
        dulist(jju,iatom,jnbor,1).re = dsfac * ulist(jju,iatom,jnbor).re * uy +
                                  sfac * dulist(jju,iatom,jnbor,1).re;
        dulist(jju,iatom,jnbor,1).im = dsfac * ulist(jju,iatom,jnbor).im * uy +
                                  sfac * dulist(jju,iatom,jnbor,1).im;
        dulist(jju,iatom,jnbor,2).re = dsfac * ulist(jju,iatom,jnbor).re * uz +
                                  sfac * dulist(jju,iatom,jnbor,2).re;
        dulist(jju,iatom,jnbor,2).im = dsfac * ulist(jju,iatom,jnbor).im * uz +
                                  sfac * dulist(jju,iatom,jnbor,2).im;

        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

template<class DeviceType>
inline
double SNAKokkos<DeviceType>::factorial(int n)
{
  //if (n < 0 || n > nmaxfactorial) {
  //  char str[128];
  //  sprintf(str, "Invalid argument to factorial %d", n);
  //  error->all(FLERR, str);
  //}

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

template<class DeviceType>
const double SNAKokkos<DeviceType>::nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmaxfactorial = 167
};

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

template<class DeviceType>
inline
double SNAKokkos<DeviceType>::deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init_clebsch_gordan()
{
  auto h_cglist = Kokkos::create_mirror_view(cglist);

  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              h_cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(factorial((j1 + aa2) / 2) *
                          factorial((j1 - aa2) / 2) *
                          factorial((j2 + bb2) / 2) *
                          factorial((j2 - bb2) / 2) *
                          factorial((j  + cc2) / 2) *
                          factorial((j  - cc2) / 2) *
                          (j + 1));

            h_cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
  Kokkos::deep_copy(cglist,h_cglist);
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

template<class DeviceType>
inline
void SNAKokkos<DeviceType>::init_rootpqarray()
{
  auto h_rootpqarray = Kokkos::create_mirror_view(rootpqarray);
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      h_rootpqarray(p,q) = sqrt(static_cast<double>(p)/q);
  Kokkos::deep_copy(rootpqarray,h_rootpqarray);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
inline
int SNAKokkos<DeviceType>::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  ndoubles = nelements*nelements;
  ntriples = nelements*nelements*nelements;
  if (chem_flag) ncount *= ntriples;

  return ncount;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double SNAKokkos<DeviceType>::compute_sfac(double r, double rcut)
{
  if (switch_flag == 0) return 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double SNAKokkos<DeviceType>::compute_dsfac(double r, double rcut)
{
  if (switch_flag == 0) return 0.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType>::compute_s_dsfac(const double r, const double rcut, double& sfac, double& dsfac) {
  if (switch_flag == 0) { sfac = 0.; dsfac = 0.; }
  else if (switch_flag == 1) {
    if (r <= rmin0) { sfac = 1.0; dsfac = 0.0; }
    else if (r > rcut) { sfac = 0.; dsfac = 0.; }
    else {
      const double rcutfac = MY_PI / (rcut - rmin0);
      double sn, cs;
      sincos((r - rmin0) * rcutfac, &sn, &cs);
      sfac = 0.5 * (cs + 1.0);
      dsfac = -0.5 * sn * rcutfac;

    }
  } else { sfac = 0.; dsfac = 0.; }
}

/* ---------------------------------------------------------------------- */

// efficient complex FMA (i.e., y += a x)
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::caxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y) {
  y.re += a.re * x.re;
  y.re -= a.im * x.im;
  y.im += a.im * x.re;
  y.im += a.re * x.im;
}

/* ---------------------------------------------------------------------- */

// efficient complex FMA, conjugate of scalar (i.e.) y += (a.re - i a.im) x)
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::caconjxpy(const SNAcomplex& a, const SNAcomplex& x, SNAcomplex& y) {
  y.re += a.re * x.re;
  y.re += a.im * x.im;
  y.im -= a.im * x.re;
  y.im += a.re * x.im;
}

/* ---------------------------------------------------------------------- */

// set direction of batched Duidrj
template<class DeviceType>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType>::set_dir(int dir_) {
  dir = dir_;
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double SNAKokkos<DeviceType>::memory_usage()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  double bytes;

  bytes = 0;

  bytes += jdimpq*jdimpq * sizeof(double);               // pqarray
  bytes += idxcg_max * sizeof(double);                   // cglist

#ifdef LMP_KOKKOS_GPU
  if (!host_flag) {

    auto natom_pad = (natom+32-1)/32;

    bytes += natom * idxu_half_max * nelements * sizeof(double);     // ulisttot_re
    bytes += natom * idxu_half_max * nelements * sizeof(double);     // ulisttot_im
    bytes += natom_pad * idxu_max * nelements * sizeof(double) * 2;  // ulisttot_pack

    bytes += natom_pad * idxz_max * ndoubles * sizeof(double) * 2;   // zlist_pack
    bytes += natom_pad * idxb_max * ntriples * sizeof(double);       // blist_pack

    bytes += natom_pad * idxu_half_max * nelements * sizeof(double); // ylist_pack_re
    bytes += natom_pad * idxu_half_max * nelements * sizeof(double); // ylist_pack_im
    bytes += natom * idxu_half_max * nelements * sizeof(double) * 2; // ylist
  } else {
#endif

    bytes += natom * nmax * idxu_cache_max * sizeof(double) * 2;     // ulist
    bytes += natom * idxu_half_max * nelements * sizeof(double) * 2; // ulisttot

    bytes += natom * idxz_max * ndoubles * sizeof(double) * 2;       // zlist
    bytes += natom * idxb_max * ntriples * sizeof(double);           // blist

    bytes += natom * idxu_half_max * nelements * sizeof(double) * 2; // ylist

    bytes += natom * nmax * idxu_cache_max * 3 * sizeof(double) * 2; // dulist
#ifdef LMP_KOKKOS_GPU
  }
#endif

  bytes += natom * nmax * 3 * sizeof(double);            // dedr

  bytes += jdim * jdim * jdim * sizeof(int);             // idxcg_block
  bytes += jdim * sizeof(int);                           // idxu_block
  bytes += jdim * sizeof(int);                           // idxu_half_block
  bytes += jdim * sizeof(int);                           // idxu_cache_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxz_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxb_block

  bytes += idxz_max * 10 * sizeof(int);                  // idxz
  bytes += idxb_max * 3 * sizeof(int);                   // idxb

  bytes += jdim * sizeof(double);                        // bzero

  bytes += natom * nmax * 3 * sizeof(double);            // rij
  bytes += natom * nmax * sizeof(int);                   // inside
  bytes += natom * nmax * sizeof(double);                // wj
  bytes += natom * nmax * sizeof(double);                // rcutij

  return bytes;
}

} // namespace LAMMPS_NS
