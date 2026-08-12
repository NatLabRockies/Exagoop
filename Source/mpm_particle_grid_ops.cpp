// clang-format off
#include "mpm_specs.H"
#include <mpm_particle_container.H>
#include <interpolants.H>
// clang-format on

/**
 * @brief Computes the total number of rigid-body particles belonging to a given
 * body ID.
 *
 * Performs a parallel reduction over all particles and counts those with:
 *   - phase = 1 (rigid)
 *   - rigid_body_id = body_id
 *
 * @param[in]  body_id     Rigid body identifier to count.
 * @param[out] total_num   Total number of rigid particles for this body.
 *
 * @return None.
 */

void MPMParticleContainer::Calculate_Total_Number_of_rigid_particles(
    int body_id, int &total_num)
{
    total_num = 0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_num = static_cast<int>(
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 1 and
                                 p.idata(intData::rigid_body_id) == body_id)
                             {
                                 return (1);
                             }
                             else
                             {
                                 return (0);
                             }
                         }));

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceIntSum(total_num);
#endif
}

/**
 * @brief Computes the total number of material (non-rigid) particles.
 *
 * Counts all particles with phase = 0 using a parallel reduction.
 *
 * @param[out] total_num  Total number of material particles.
 *
 * @return None.
 */

void MPMParticleContainer::Calculate_Total_Number_of_MaterialParticles(
    int &total_num)
{
    total_num = 0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_num = static_cast<int>(
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 0)
                             {
                                 return (1);
                             }
                             else
                             {
                                 return (0);
                             }
                         }));

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceIntSum(total_num);
#endif
}

/**
 * @brief Computes the total mass of rigid-body particles for a given body ID.
 *
 * Sums p.mass for all particles satisfying:
 *   - phase = 1
 *   - rigid_body_id = body_id
 *
 * @param[in]  body_id      Rigid body identifier.
 * @param[out] total_mass   Total mass of rigid particles for this body.
 *
 * @return None.
 */

void MPMParticleContainer::Calculate_Total_Mass_RigidParticles(int body_id,
                                                               Real &total_mass)
{
    total_mass = 0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_mass =
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 1 and
                                 p.idata(intData::rigid_body_id) == body_id)
                             {
                                 return (p.rdata(realData::mass));
                             }
                             else
                             {
                                 return (0.0);
                             }
                         });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(total_mass);
#endif
}

/**
 * @brief Computes the total mass of all material (phase=0) particles.
 *
 * Performs a parallel reduction summing p.mass for all material points.
 *
 * @param[out] total_mass  Total mass of material particles.
 *
 * @return None.
 */

void MPMParticleContainer::Calculate_Total_Mass_MaterialPoints(Real &total_mass)
{

    total_mass = 0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_mass =
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 0)
                             {
                                 return (p.rdata(realData::mass));
                             }
                             else
                             {
                                 return (0.0);
                             }
                         });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(total_mass);
#endif
}

/**
 * @brief Computes the total volume of all material (phase=0) particles.
 *
 * Performs a parallel reduction summing p.volume for all material points.
 *
 * @param[out] total_vol  Total volume of material particles.
 *
 * @return None.
 */

void MPMParticleContainer::Calculate_Total_Vol_MaterialPoints(Real &total_vol)
{

    total_vol = 0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_vol =
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 0)
                             {
                                 return (p.rdata(realData::volume));
                             }
                             else
                             {
                                 return (0.0);
                             }
                         });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(total_vol);
#endif
}

/**
 * @brief Computes the total volume of rigid-body particles for a given body ID.
 *
 * Sums p.volume for all particles satisfying:
 *   - phase = 1
 *   - rigid_body_id = body_id
 *
 * @param[in] body_id  Rigid body identifier.
 *
 * @return amrex::Real
 *         Total volume of rigid particles belonging to this body.
 */

amrex::Real
MPMParticleContainer::Calculate_Total_Vol_RigidParticles(int body_id)
{

    amrex::Real total_vol = 0.0;

    using PType = typename MPMParticleContainer::SuperParticleType;
    total_vol =
        amrex::ReduceSum(*this,
                         [=] AMREX_GPU_HOST_DEVICE(const PType &p) -> Real
                         {
                             if (p.idata(intData::phase) == 1 and
                                 p.idata(intData::rigid_body_id) == body_id)
                             {
                                 return (p.rdata(realData::volume));
                             }
                             else
                             {
                                 return (0.0);
                             }
                         });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealSum(total_vol);
#endif
    return (total_vol);
}

/**
 * @brief Computes the stencil bounds for interpolation in one dimension.
 *
 * Determines the lower and upper stencil offsets (bmin, bmax) for a given
 * interpolation scheme:
 *
 *   - scheme = 1 → linear (2-point support)
 *   - scheme = 2 or 3 → quadratic/cubic spline (4-point support)
 *
 * Boundary handling:
 *   - Periodic → symmetric stencil [-1, 2]
 *   - Non-periodic → one-sided stencils at domain boundaries
 *
 * @param[in] ivd          Cell index in this dimension.
 * @param[in] lod          Domain lower bound.
 * @param[in] hid          Domain upper bound.
 * @param[in] scheme       Interpolation order (1, 2, or 3).
 * @param[in] is_periodic  Whether the dimension is periodic.
 *
 * @return std::pair<int,int>
 *         (bmin, bmax) stencil bounds.
 */

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE std::pair<int, int>
compute_bounds(int ivd, int lod, int hid, int scheme, bool is_periodic)
{
    int bmin = 0, bmax = 0;

    if (scheme == 1)
    {

        bmin = 0;
        bmax = 2;
    }
    else if (scheme == 3 or scheme == 2)
    {
        if (is_periodic)
        {
            bmin = -1;
            bmax = 3; // symmetric
        }
        else
        {
            if (ivd == lod)
            {
                bmin = 0;
                bmax = 3; // one-sided at lo
            }
            else if (ivd == hid)
            {
                bmin = -1;
                bmax = bmin + 3; // one-sided at hi
            }
            else
            {
                bmin = -1;
                bmax = 3; // interior symmetric
            }
        }
    }
    else
    {
        amrex::Abort("Unsupported interpolation scheme");
    }

    return {bmin, bmax};
}

/**
 * @brief Deposits particle mass, momentum, and internal/external forces onto
 * the nodal grid.
 *
 * For each particle:
 *   1. Determines its cell index and interpolation stencil.
 *   2. Evaluates basis functions and (optionally) basis gradients.
 *   3. Deposits:
 *        - mass                (if update_mass)
 *        - momentum            (if update_vel)
 *        - body forces         (gravity, external loads)
 *        - internal forces     (−V σ ∇N)  if update_forces
 *
 * After deposition, nodal velocities and stresses are normalized by nodal mass.
 *
 * @param[in,out] nodaldata               MultiFab storing nodal quantities.
 * @param[in]     gravity                 Gravity vector.
 * @param[in]     external_loads_present  Flag enabling external slab forces.
 * @param[in]     force_slab_lo           Lower bounds of external force region.
 * @param[in]     force_slab_hi           Upper bounds of external force region.
 * @param[in]     extforce                External force vector.
 * @param[in]     update_mass             Whether to deposit mass.
 * @param[in]     update_vel              Whether to deposit momentum.
 * @param[in]     update_forces           Whether to deposit forces (1=body,
 * 2=stress).
 * @param[in]     mass_tolerance          Threshold for velocity normalization.
 * @param[in]     order_scheme_directional  Per-dimension interpolation order.
 * @param[in]     periodic                Per-dimension periodicity flags.
 *
 * @return None.
 */

void MPMParticleContainer::deposit_onto_grid_momentum(
    MultiFab &nodaldata,
    amrex::Array<Real, AMREX_SPACEDIM> gravity,
    int external_loads_present,
    amrex::Array<Real, AMREX_SPACEDIM> force_slab_lo,
    amrex::Array<Real, AMREX_SPACEDIM> force_slab_hi,
    amrex::Array<Real, AMREX_SPACEDIM> extforce,
    int update_mass,
    int update_vel,
    int update_forces,
    amrex::Real mass_tolerance,
    amrex::GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
    amrex::GpuArray<int, AMREX_SPACEDIM> periodic)
{

    const int lev = 0;
    const Geometry &geom = Geom(lev);
    auto &plev = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const int *loarr = domain.loVect();
    const int *hiarr = domain.hiVect();

    Real grav[] = {AMREX_D_DECL(gravity[XDIR], gravity[YDIR], gravity[ZDIR])};
    Real extpforce[] = {
        AMREX_D_DECL(extforce[XDIR], extforce[YDIR], extforce[ZDIR])};

    int lo[] = {AMREX_D_DECL(loarr[0], loarr[1], loarr[2])};
    int hi[] = {AMREX_D_DECL(hiarr[0], hiarr[1], hiarr[2])};

    int extloads = external_loads_present;

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {

        const Box &nodalbox =
            convert(mfi.validbox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        auto nodal_data_arr = nodaldata.array(mfi);
        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeindex(AMREX_D_DECL(i, j, k));
                if (update_mass)
                {
                    nodal_data_arr(nodeindex, MASS_INDEX) = 0.0;
                }
                if (update_vel)
                {

                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        nodal_data_arr(nodeindex, VELX_INDEX + d) = 0.0;
                    }
                }
                if (update_forces)
                {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        nodal_data_arr(nodeindex, FRCX_INDEX + d) = 0.0;
                    }
                }
                if (update_forces == 2)
                {
                    nodal_data_arr(nodeindex, STRESS_INDEX) = 0.0;
                }
            });
    }

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {

        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        auto &ptile = plev[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto &aos = ptile.GetArrayOfStructs();
        int nt = aos.numRealParticles() + aos.numNeighborParticles();
        auto nodal_data_arr = nodaldata.array(mfi);
        ParticleType *pstruct = aos().dataPtr();

        amrex::ParallelFor(
            nt,
            [=] AMREX_GPU_DEVICE(int ip) noexcept
            {
                ParticleType &p = pstruct[ip];
                if (p.id() < 0)
                    return;
                if (p.idata(intData::phase) != 0)
                    return;

                amrex::Real xp[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    xp[d] = p.pos(d);

                auto iv = getParticleCell(p, plo, dxi, domain);

                // Compute stencil extents per dimension

                int min_idx[AMREX_SPACEDIM], max_idx[AMREX_SPACEDIM];

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    auto bounds = compute_bounds(iv[d], lo[d], hi[d],
                                                 order_scheme_directional[d],
                                                 periodic[d]);
                    min_idx[d] =
                        amrex::max(bounds.first, nodalbox.smallEnd(d) - iv[d]);
                    max_idx[d] = amrex::min(bounds.second,
                                            nodalbox.bigEnd(d) - iv[d] + 1);
                }

                // Nested loops over stencil (specialized per dimension)
                amrex::Real basisvalue = 0.0;
                amrex::Real basisval_grad[AMREX_SPACEDIM] = {
                    AMREX_D_DECL(0.0, 0.0, 0.0)};
                IntVect ivlocal(AMREX_D_DECL(0, 0, 0));

#if (AMREX_SPACEDIM == 3)
                for (int n = min_idx[2]; n < max_idx[2]; n++)
                {
#endif
#if (AMREX_SPACEDIM >= 2)
                    for (int m = min_idx[1]; m < max_idx[1]; m++)
                    {
#endif
                        for (int l = min_idx[0]; l < max_idx[0]; l++)
                        {
                            ivlocal = IntVect(
                                AMREX_D_DECL(iv[0] + l, iv[1] + m, iv[2] + n));

                            IntVect stencil(AMREX_D_DECL(l, m, n));
                            {

                                if (update_forces)
                                {
                                    basisvalue = basisval_and_grad(
                                        stencil, iv, xp, plo, dx,
                                        order_scheme_directional, periodic, lo,
                                        hi, basisval_grad);
                                }
                                else
                                {
                                    basisvalue =
                                        basisval(stencil, iv, xp, plo, dx,
                                                 order_scheme_directional,
                                                 periodic, lo, hi);
                                }

                                if (update_mass)
                                {
                                    amrex::Real mass_contrib =
                                        p.rdata(realData::mass) * basisvalue;
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal, MASS_INDEX),
                                        mass_contrib);
                                }
                                if (update_vel)
                                {

                                    amrex::Real p_contrib[AMREX_SPACEDIM] = {
                                        AMREX_D_DECL(
                                            p.rdata(realData::mass) *
                                                p.rdata(realData::xvel) *
                                                basisvalue,
                                            p.rdata(realData::mass) *
                                                p.rdata(realData::yvel) *
                                                basisvalue,
                                            p.rdata(realData::mass) *
                                                p.rdata(realData::zvel) *
                                                basisvalue)};

                                    for (int dim = 0; dim < AMREX_SPACEDIM;
                                         dim++)
                                    {
                                        amrex::Gpu::Atomic::AddNoRet(
                                            &nodal_data_arr(ivlocal,
                                                            VELX_INDEX + dim),
                                            p_contrib[dim]);
                                    }
                                }

                                if (update_forces)
                                {

                                    amrex::Real stress_tens[AMREX_SPACEDIM *
                                                            AMREX_SPACEDIM];
                                    get_tensor(p, realData::stress,
                                               stress_tens);

                                    amrex::Real bforce_contrib[AMREX_SPACEDIM] =
                                        {AMREX_D_DECL(
                                            p.rdata(realData::mass) *
                                                grav[XDIR] * basisvalue,
                                            p.rdata(realData::mass) *
                                                grav[YDIR] * basisvalue,
                                            p.rdata(realData::mass) *
                                                grav[ZDIR] * basisvalue)};

                                    if (extloads &&
                                        xp[XDIR] > force_slab_lo[XDIR] &&
                                        xp[XDIR] < force_slab_hi[XDIR]
#if (AMREX_SPACEDIM >= 2)
                                        && xp[YDIR] > force_slab_lo[YDIR] &&
                                        xp[YDIR] < force_slab_hi[YDIR]
#endif
#if (AMREX_SPACEDIM == 3)
                                        && xp[ZDIR] > force_slab_lo[ZDIR] &&
                                        xp[ZDIR] < force_slab_hi[ZDIR]
#endif
                                    )
                                    {
                                        for (int dim = 0; dim < AMREX_SPACEDIM;
                                             dim++)
                                        {
                                            bforce_contrib[dim] +=
                                                extpforce[dim] * basisvalue;
                                        }
                                    }

                                    amrex::Real tensvect[AMREX_SPACEDIM];
                                    tensor_vector_pdt(stress_tens,
                                                      basisval_grad, tensvect);

                                    amrex::Real
                                        intforce_contrib[AMREX_SPACEDIM] = {
                                            AMREX_D_DECL(
                                                -p.rdata(realData::volume) *
                                                    tensvect[XDIR],
                                                -p.rdata(realData::volume) *
                                                    tensvect[YDIR],
                                                -p.rdata(realData::volume) *
                                                    tensvect[ZDIR])};

                                    for (int dim = 0; dim < AMREX_SPACEDIM;
                                         dim++)
                                    {
                                        amrex::Gpu::Atomic::AddNoRet(
                                            &nodal_data_arr(ivlocal,
                                                            FRCX_INDEX + dim),
                                            bforce_contrib[dim] +
                                                intforce_contrib[dim]);
                                    }
                                }
                            }
                        }
#if (AMREX_SPACEDIM >= 2)
                    }
#endif
#if (AMREX_SPACEDIM == 3)
                }
#endif
            });
    }

    // Normalize velocities and stresses
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {

        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        auto nodal_data_arr = nodaldata.array(mfi);
        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeindex(AMREX_D_DECL(i, j, k));
                if (update_vel && nodal_data_arr(nodeindex, MASS_INDEX) > 0.0)
                {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        if (nodal_data_arr(nodeindex, MASS_INDEX) >=
                            mass_tolerance)
                        {
                            nodal_data_arr(nodeindex, VELX_INDEX + d) /=
                                nodal_data_arr(nodeindex, MASS_INDEX);
                        }

                        else
                            nodal_data_arr(nodeindex, VELX_INDEX + d) = 0.0;
                    }
                }
                if (update_forces == 2 &&
                    nodal_data_arr(nodeindex, MASS_INDEX) > 0.0)
                {
                    if (nodal_data_arr(nodeindex, MASS_INDEX) >= mass_tolerance)
                        nodal_data_arr(nodeindex, STRESS_INDEX) /=
                            nodal_data_arr(nodeindex, MASS_INDEX);
                    else
                        nodal_data_arr(nodeindex, STRESS_INDEX) = 0.0;
                }
            });
    }
}

/**
 * @brief Deposits thermal quantities (mass*specific_heat, temperature, heat
 * sources) onto the nodal grid.
 *
 * For each particle:
 *   - Deposits mass * specific_heat
 *   - Deposits mass * specific_heat * temperature
 *   - Deposits internal heat flux contributions (−∇·q)
 *   - Deposits external heat sources
 *
 * After deposition, nodal temperature is computed as:
 *      T = (mass_spheat_temp) / (mass_spheat)
 *
 * @param[in,out] nodaldata               Nodal thermal MultiFab.
 * @param[in]     resetnodaldata_to_zero  Whether to zero nodal fields first.
 * @param[in]     update_mass_temp        Whether to deposit thermal mass.
 * @param[in]     update_source           Whether to deposit heat sources.
 * @param[in]     mass_tolerance          Threshold for temperature
 * normalization.
 * @param[in]     order_scheme_directional  Per-dimension interpolation order.
 * @param[in]     periodic                Per-dimension periodicity flags.
 *
 * @return None.
 */

#if USE_TEMP
void MPMParticleContainer::deposit_onto_grid_temperature(
    MultiFab &nodaldata,
    bool resetnodaldata_to_zero,
    bool update_mass_temp,
    bool update_source,
    amrex::Real mass_tolerance,
    GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
    GpuArray<int, AMREX_SPACEDIM> periodic)
{

    const int lev = 0;
    const Geometry &geom = Geom(lev);
    auto &plev = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const int *loarr = domain.loVect();
    const int *hiarr = domain.hiVect();

    int lo[] = {AMREX_D_DECL(loarr[0], loarr[1], loarr[2])};
    int hi[] = {AMREX_D_DECL(hiarr[0], hiarr[1], hiarr[2])};

    // Zero out nodal data
    if (resetnodaldata_to_zero)
    {
        for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
        {

            const Box &nodalbox =
                convert(mfi.validbox(), IntVect(AMREX_D_DECL(1, 1, 1)));
            auto nodal_data_arr = nodaldata.array(mfi);
            amrex::ParallelFor(
                nodalbox,
                [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
                {
                    IntVect nodeindex(AMREX_D_DECL(i, j, k));
                    nodal_data_arr(nodeindex, MASS_SPHEAT) = shunya;
                    nodal_data_arr(nodeindex, MASS_SPHEAT_TEMP) = shunya;
                    nodal_data_arr(nodeindex, SOURCE_TEMP_INDEX) = shunya;
                    if (update_mass_temp)
                        nodal_data_arr(nodeindex, MASS_CONDUCTIVITY) = shunya;
                });
        }
    }

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {

        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        auto &ptile = plev[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto &aos = ptile.GetArrayOfStructs();
        int nt = aos.numRealParticles() + aos.numNeighborParticles();
        auto nodal_data_arr = nodaldata.array(mfi);
        ParticleType *pstruct = aos().dataPtr();

        amrex::ParallelFor(
            nt,
            [=] AMREX_GPU_DEVICE(int ip) noexcept
            {
                ParticleType &p = pstruct[ip];
                if (p.id() < 0)
                    return; // skip invalidated particles (NaN-init guard)
                if (p.idata(intData::phase) != 0)
                    return;

                amrex::Real xp[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    xp[d] = p.pos(d);

                auto iv = getParticleCell(p, plo, dxi, domain);

                // Compute stencil extents per dimension

                int min_idx[AMREX_SPACEDIM], max_idx[AMREX_SPACEDIM];

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    auto bounds = compute_bounds(iv[d], lo[d], hi[d],
                                                 order_scheme_directional[d],
                                                 periodic[d]);

                    min_idx[d] =
                        amrex::max(bounds.first, nodalbox.smallEnd(d) - iv[d]);
                    max_idx[d] = amrex::min(bounds.second,
                                            nodalbox.bigEnd(d) - iv[d] + 1);
                }

                // Nested loops over stencil (specialized per dimension)
                amrex::Real basisvalue = 0.0;
                amrex::Real basisval_grad[AMREX_SPACEDIM] = {
                    AMREX_D_DECL(0.0, 0.0, 0.0)};
                IntVect ivlocal(AMREX_D_DECL(0, 0, 0));

#if (AMREX_SPACEDIM == 3)
                for (int n = min_idx[2]; n < max_idx[2]; n++)
                {
#endif
#if (AMREX_SPACEDIM >= 2)
                    for (int m = min_idx[1]; m < max_idx[1]; m++)
                    {
#endif
                        for (int l = min_idx[0]; l < max_idx[0]; l++)
                        {
                            ivlocal = IntVect(
                                AMREX_D_DECL(iv[0] + l, iv[1] + m, iv[2] + n));
                            IntVect stencil(AMREX_D_DECL(l, m, n));
                            {
                                basisvalue = basisval(stencil, iv, xp, plo, dx,
                                                      order_scheme_directional,
                                                      periodic, lo, hi);

                                if (update_mass_temp)
                                {
                                    amrex::Real mass_spheat_contrib =
                                        p.rdata(realData::mass) *
                                        p.rdata(realData::specific_heat) *
                                        basisvalue;
                                    amrex::Real mass_spheat_temp_contrib =
                                        p.rdata(realData::mass) *
                                        p.rdata(realData::specific_heat) *
                                        p.rdata(realData::temperature) *
                                        basisvalue;
                                    amrex::Real mass_conductivity_contrib =
                                        p.rdata(realData::mass) *
                                        p.rdata(
                                            realData::thermal_conductivity) *
                                        basisvalue;
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal, MASS_SPHEAT),
                                        mass_spheat_contrib);
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal,
                                                        MASS_SPHEAT_TEMP),
                                        mass_spheat_temp_contrib);
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal,
                                                        MASS_CONDUCTIVITY),
                                        mass_conductivity_contrib);
                                }

                                if (update_source)
                                {
                                    for (int dim = 0; dim < AMREX_SPACEDIM;
                                         dim++)
                                    {
                                        basisval_grad[dim] = basisvalder(
                                            dim, stencil, iv, xp, plo, dx,
                                            order_scheme_directional, periodic,
                                            lo, hi);
                                    }

                                    amrex::Real heat_flux_vect[AMREX_SPACEDIM];
                                    amrex::Real net_heatflux = 0.0;
                                    amrex::Real int_source = 0.0;
                                    amrex::Real ext_source = 0.0;
                                    get_vector(p, realData::heat_flux,
                                               heat_flux_vect);

                                    for (int dim = 0; dim < AMREX_SPACEDIM;
                                         dim++)
                                    {
                                        net_heatflux += heat_flux_vect[dim] *
                                                        basisval_grad[dim];
                                    }

                                    int_source = net_heatflux *
                                                 p.rdata(realData::volume);
                                    ext_source =
                                        p.rdata(realData::volume) *
                                        p.rdata(realData::heat_source) *
                                        basisvalue;

                                    amrex::Gpu::Atomic::AddNoRet(
                                        &nodal_data_arr(ivlocal,
                                                        SOURCE_TEMP_INDEX),
                                        int_source + ext_source);
                                }
                            }
                        }
#if (AMREX_SPACEDIM >= 2)
                    }
#endif
#if (AMREX_SPACEDIM == 3)
                }
#endif
            });
    }

    if (update_mass_temp)
    {
        nodaldata.SumBoundary(MASS_SPHEAT, 2, geom.periodicity());
        nodaldata.SumBoundary(MASS_CONDUCTIVITY, 1, geom.periodicity());
    }
    if (update_source)
    {
        nodaldata.SumBoundary(SOURCE_TEMP_INDEX, 1, geom.periodicity());
    }

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        auto nodal_data_arr = nodaldata.array(mfi);
        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeindex(AMREX_D_DECL(i, j, k));
                if (update_mass_temp &&
                    nodal_data_arr(nodeindex, MASS_SPHEAT) > 0.0)
                {

                    if (nodal_data_arr(nodeindex, MASS_SPHEAT) >=
                        mass_tolerance)
                    {
                        nodal_data_arr(nodeindex, TEMPERATURE) =
                            nodal_data_arr(nodeindex, MASS_SPHEAT_TEMP) /
                            nodal_data_arr(nodeindex, MASS_SPHEAT);
                    }

                    else
                    {
                        nodal_data_arr(nodeindex, TEMPERATURE) = 0.0;
                    }
                }
            });
    }
}

#endif

/**
 * @brief Interpolates nodal grid quantities back to particles (G2P transfer).
 *
 * This routine performs the particle update step of the MPM algorithm.
 * For each material particle (phase = 0), it:
 *
 *   1. Determines the particle’s cell index and interpolation stencil.
 *   2. Interpolates **velocity** using:
 *        - linear (order=1)
 *        - quadratic (order=2)
 *        - cubic (order=3)
 *      depending on the directional order scheme.
 *
 *   3. Applies **PIC/FLIP blending**:
 *        vᵖ ← α vᵖ + α Δv_grid + (1−α) v_grid
 *      where α = alpha_pic_flip.
 *
 *   4. If enabled, computes **velocity gradient** ∇v using basis derivatives.
 *
 *   5. Updates the **deformation gradient**:
 *        Fᵖ ← (I + ∇v dt) Fᵖ
 *
 *   6. Computes the **strain‑rate tensor**:
 *        Dᵖ = 0.5 (∇v + ∇vᵀ)
 *
 * @param[in]     nodaldata               MultiFab containing nodal grid values.
 * @param[in]     update_vel              Whether to interpolate velocity.
 * @param[in]     update_strainrate       Whether to compute strain‑rate and F.
 * @param[in]     order_scheme_directional  Per‑dimension interpolation order.
 * @param[in]     periodic                Per‑dimension periodicity flags.
 * @param[in]     alpha_pic_flip          PIC/FLIP blending parameter.
 * @param[in]     dt                      Time step for deformation update.
 *
 * @return None.
 */

void MPMParticleContainer::interpolate_from_grid(
    MultiFab &nodaldata,
    int update_vel,
    int update_strainrate,
    GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
    GpuArray<int, AMREX_SPACEDIM> periodic,
    amrex::Real alpha_pic_flip,
    amrex::Real dt)
{
    const int lev = 0;
    const Geometry &geom = Geom(lev);
    auto &plev = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const int *loarr = domain.loVect();
    const int *hiarr = domain.hiVect();

    int lo[AMREX_SPACEDIM];
    int hi[AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        lo[d] = loarr[d];
        hi[d] = hiarr[d];
    }

    nodaldata.FillBoundary(geom.periodicity());

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto &ptile = plev[index];
        auto &aos = ptile.GetArrayOfStructs();
        const int np = aos.numRealParticles();

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);
        ParticleType *pstruct = aos().dataPtr();

        amrex::ParallelFor(
            np,
            [=] AMREX_GPU_DEVICE(int i) noexcept
            {
                ParticleType &p = pstruct[i];
                if (p.id() < 0)
                    return; // skip invalidated particles (NaN-init guard)
                if (p.idata(intData::phase) != 0)
                    return;

                amrex::Real xp[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    xp[d] = p.pos(d);
                }

                auto iv = getParticleCell(p, plo, dxi, domain);

                // Dimension‑aware min/max indices

                IntVect min_index(AMREX_D_DECL(0, 0, 0));
                IntVect max_index(AMREX_D_DECL(0, 0, 0));

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    auto bounds = compute_bounds(iv[d], lo[d], hi[d],
                                                 order_scheme_directional[d],
                                                 periodic[d]);
                    min_index[d] = bounds.first;
                    max_index[d] = bounds.second;
                }

                if (update_vel)
                {

                    amrex::Real interp_vals[2];
                    for (int dim = 0; dim < AMREX_SPACEDIM; dim++)
                    {

                        if (order_scheme_directional[dim] == 1)
                        {
                            bilin_interp_two(xp, iv, plo, dx, nodal_data_arr,
                                             VELX_INDEX + dim,
                                             DELTA_VELX_INDEX + dim,
                                             interp_vals);
                        }
                        else if (order_scheme_directional[dim] == 2)
                        {
                            quadratic_interp_two(xp, iv, min_index, max_index,
                                                 plo, dx, nodal_data_arr,
                                                 VELX_INDEX + dim,
                                                 DELTA_VELX_INDEX + dim,
                                                 periodic, lo, hi, interp_vals);
                        }
                        else if (order_scheme_directional[dim] == 3)
                        {
                            cubic_interp_two(xp, iv, min_index, max_index, plo,
                                             dx, nodal_data_arr,
                                             VELX_INDEX + dim,
                                             DELTA_VELX_INDEX + dim, periodic,
                                             lo, hi, interp_vals);
                        }

                        p.rdata(realData::xvel_prime + dim) = interp_vals[0];
                        p.rdata(realData::xvel + dim) =
                            alpha_pic_flip * p.rdata(realData::xvel + dim) +
                            alpha_pic_flip * interp_vals[1] +
                            (1 - alpha_pic_flip) * interp_vals[0];
                    }
                }

                // Strain rate update
                if (update_strainrate)
                {
                    amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM] = {};
#if (AMREX_SPACEDIM == 3)
                    for (int n = min_index[2]; n < max_index[2]; ++n)
                    {
#endif
#if (AMREX_SPACEDIM != 1)
                        for (int m = min_index[1]; m < max_index[1]; ++m)
                        {
#endif
                            for (int l = min_index[0]; l < max_index[0]; ++l)
                            {
                                IntVect stencil(AMREX_D_DECL(l, m, n));
                                amrex::Real basisval_grad[AMREX_SPACEDIM];
                                basisval_and_grad(stencil, iv, xp, plo, dx,
                                                  order_scheme_directional,
                                                  periodic, lo, hi,
                                                  basisval_grad);

                                IntVect nodeindex(AMREX_D_DECL(
                                    iv[0] + l, iv[1] + m, iv[2] + n));
                                for (int d1 = 0; d1 < AMREX_SPACEDIM; ++d1)
                                {
                                    amrex::Real v_d1 = nodal_data_arr(
                                        nodeindex, VELX_INDEX + d1);
                                    for (int d2 = 0; d2 < AMREX_SPACEDIM; ++d2)
                                    {
                                        gradvp[d1][d2] +=
                                            v_d1 * basisval_grad[d2];
                                    }
                                }
                            }
#if (AMREX_SPACEDIM != 1)
                        }
#endif
#if (AMREX_SPACEDIM == 3)
                    }
#endif

                    // Calculate deformation gradient tensor. F_p^{t+dt} at time
                    // t+dt
                    get_deformation_gradient_tensor(
                        p, realData::deformation_gradient, gradvp, dt);

                    // Calculate strain rate tensor. D_p^{t+dt} at time t+dt
                    int ind = 0;
                    for (int d1 = 0; d1 < AMREX_SPACEDIM; ++d1)
                    {
                        for (int d2 = d1; d2 < AMREX_SPACEDIM; ++d2)
                        {
                            p.rdata(realData::strainrate + ind) =
                                0.5 * (gradvp[d1][d2] + gradvp[d2][d1]);
                            ind++;
                        }
                    }
                }
            });
    }
}

#if USE_TEMP
/**
 * @brief Interpolates temperature and heat‑flux from nodal grid to particles.
 *
 * For each material particle (phase = 0), this routine:
 *
 *   1. Determines the interpolation stencil based on the directional order.
 *
 *   2. If update_temperature is enabled:
 *        Tᵖ ← Tᵖ + ΔT_grid(xᵖ)
 *      using linear, quadratic, or cubic interpolation.
 *
 *   3. If update_heatflux is enabled:
 *        Computes ∇T using basis derivatives:
 *            ∇T = Σ (T_node ∇N)
 *        Then computes particle heat flux:
 *            qᵖ = −k ∇T
 *
 * @param[in]     nodaldata               MultiFab containing nodal temperature.
 * @param[in]     update_temperature      Whether to interpolate ΔT.
 * @param[in]     update_heatflux         Whether to compute heat‑flux.
 * @param[in]     order_scheme_directional  Per‑dimension interpolation order.
 * @param[in]     periodic                Per‑dimension periodicity flags.
 * @param[in]     alpha_pic_flip          PIC/FLIP blending parameter (unused
 * here).
 *
 * @return None.
 */

void MPMParticleContainer::interpolate_from_grid_temperature(
    MultiFab &nodaldata,
    bool update_temperature,
    bool update_heatflux,
    GpuArray<int, AMREX_SPACEDIM> order_scheme_directional,
    GpuArray<int, AMREX_SPACEDIM> periodic,
    [[maybe_unused]] amrex::Real alpha_pic_flip)
{
    const int lev = 0;
    const Geometry &geom = Geom(lev);
    auto &plev = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const int *loarr = domain.loVect();
    const int *hiarr = domain.hiVect();

    int lo[AMREX_SPACEDIM];
    int hi[AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        lo[d] = loarr[d];
        hi[d] = hiarr[d];
    }

    nodaldata.FillBoundary(geom.periodicity());

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto &ptile = plev[index];
        auto &aos = ptile.GetArrayOfStructs();
        const int np = aos.numRealParticles();

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);
        ParticleType *pstruct = aos().dataPtr();

        amrex::ParallelFor(
            np,
            [=] AMREX_GPU_DEVICE(int i) noexcept
            {
                ParticleType &p = pstruct[i];
                if (p.id() < 0)
                    return;
                if (p.idata(intData::phase) != 0)
                    return;

                amrex::Real xp[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    xp[d] = p.pos(d);
                }

                auto iv = getParticleCell(p, plo, dxi, domain);

                // Dimension‑aware min/max indices

                IntVect min_index(AMREX_D_DECL(0, 0, 0));
                IntVect max_index(AMREX_D_DECL(0, 0, 0));

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    auto bounds = compute_bounds(iv[d], lo[d], hi[d],
                                                 order_scheme_directional[d],
                                                 periodic[d]);
                    min_index[d] = bounds.first;
                    max_index[d] = bounds.second;
                }

                if (update_temperature)
                {
                    if (order_scheme_directional[0] == 1)
                    {
                        p.rdata(realData::temperature) += bilin_interp(
                            xp, iv, plo, dx, nodal_data_arr, DELTA_TEMPERATURE);
                    }
                    else if (order_scheme_directional[0] == 2)
                    {
                        p.rdata(realData::temperature) += quadratic_interp(
                            xp, iv, min_index, max_index, plo, dx,
                            nodal_data_arr, DELTA_TEMPERATURE, periodic, lo,
                            hi);
                    }
                    else if (order_scheme_directional[0] == 3)
                    {
                        p.rdata(realData::temperature) +=
                            cubic_interp(xp, iv, min_index, max_index, plo, dx,
                                         nodal_data_arr, DELTA_TEMPERATURE,
                                         periodic, lo, hi);
                    }
                }

                // Strain rate update
                if (update_heatflux)
                {
                    amrex::Real gradT[AMREX_SPACEDIM] = {};
#if (AMREX_SPACEDIM == 3)
                    for (int n = min_index[2]; n < max_index[2]; ++n)
                    {
#endif
#if (AMREX_SPACEDIM != 1)
                        for (int m = min_index[1]; m < max_index[1]; ++m)
                        {
#endif
                            for (int l = min_index[0]; l < max_index[0]; ++l)
                            {
                                IntVect stencil(AMREX_D_DECL(l, m, n));
                                amrex::Real basisval_grad[AMREX_SPACEDIM];
                                basisval_and_grad(stencil, iv, xp, plo, dx,
                                                  order_scheme_directional,
                                                  periodic, lo, hi,
                                                  basisval_grad);

                                IntVect nodeindex(AMREX_D_DECL(
                                    iv[0] + l, iv[1] + m, iv[2] + n));
                                amrex::Real T_node =
                                    nodal_data_arr(nodeindex, TEMPERATURE);
                                for (int d1 = 0; d1 < AMREX_SPACEDIM; ++d1)
                                {
                                    gradT[d1] += T_node * basisval_grad[d1];
                                }
                            }
#if (AMREX_SPACEDIM != 1)
                        }
#endif
#if (AMREX_SPACEDIM == 3)
                    }
#endif
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        p.rdata(realData::heat_flux + d) =
                            -p.rdata(realData::thermal_conductivity) * gradT[d];
                    }
                }
            });
    }
}
#endif

/**
 * @brief Computes nodal normals from particle mass‑weighted basis gradients.
 *
 * This function is currently disabled (commented out), but its intended
 * behavior is:
 *
 *   1. Zero nodal normal components.
 *   2. For each particle:
 *        - Compute ∇N at the particle location.
 *        - Accumulate mass‑weighted gradients into nodal normals.
 *   3. Normalize nodal normals at each node.
 *
 * This is typically used for contact algorithms or boundary detection.
 *
 * @note The implementation is commented out and not executed.
 *
 * @return None.
 */

void MPMParticleContainer::calculate_nodal_normal(
    MultiFab & /*nodaldata*/,
    amrex::Real /*mass_tolerance*/,
    GpuArray<int, AMREX_SPACEDIM> /*order_scheme_directional*/,
    GpuArray<int, AMREX_SPACEDIM> /*periodic*/)
{
    // To be implemented
}
