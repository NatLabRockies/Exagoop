// clang-format off
#include <nodal_data_ops.H>
#include <mpm_specs.H>
#include <mpm_eb.H>
#include <mpm_kernels.H>
#include <AMReX_iMultiFab.H>
#include <cmath>
// clang-format on

using namespace amrex;

/**
 * @brief Writes nodal grid data to a single‑level plotfile.
 *
 * Converts nodal (node‑centered) data to cell‑centered form using
 * average_node_to_cellcenter(), then writes the resulting MultiFab
 * to a plotfile directory.
 *
 * @param[in] fname        Output plotfile name.
 * @param[in] nodaldata    Nodal MultiFab to be written.
 * @param[in] fieldnames   Names of each component in nodaldata.
 * @param[in] geom         Geometry describing the domain.
 * @param[in] ba           BoxArray for the output MultiFab.
 * @param[in] dm           DistributionMapping for parallel layout.
 * @param[in] time         Simulation time for the plotfile.
 *
 * @return None.
 */

void write_grid_file(std::string fname,
                     MultiFab &nodaldata,
                     Vector<std::string> fieldnames,
                     Geometry geom,
                     BoxArray ba,
                     DistributionMapping dm,
                     Real time)
{
    MultiFab plotmf(ba, dm, nodaldata.nComp(), 0);
    average_node_to_cellcenter(plotmf, 0, nodaldata, 0, nodaldata.nComp());
    WriteSingleLevelPlotfile(fname, plotmf, fieldnames, geom, time, 0);
}

/**
 * @brief Stores the current nodal mass and velocity for later PIC/FLIP updates.
 *
 * For each node with positive mass:
 *   - MASS_OLD_INDEX ← MASS_INDEX
 *   - DELTA_VELX_INDEX[d] ← VELX_INDEX[d]
 *
 * This is used to compute Δv during the next P2G/G2P cycle.
 *
 * @param[in,out] nodaldata  Nodal MultiFab containing mass and velocity fields.
 *
 * @return None.
 */

void backup_current_velocity(MultiFab &nodaldata)
{
    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (nodal_data_arr(i, j, k, MASS_INDEX) > shunya)
                {
                    nodal_data_arr(i, j, k, MASS_OLD_INDEX) =
                        nodal_data_arr(i, j, k, MASS_INDEX);
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                    {
                        nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d) =
                            nodal_data_arr(i, j, k, VELX_INDEX + d);
                    }
                }
            });
    }
}

#if USE_TEMP
/**
 * @brief Stores the current nodal temperature for later ΔT updates.
 *
 * For each node with positive thermal mass (MASS_SPHEAT):
 *   DELTA_TEMPERATURE ← TEMPERATURE
 *
 * @param[in,out] nodaldata  Nodal MultiFab containing thermal fields.
 *
 * @return None.
 */

void backup_current_temperature(MultiFab &nodaldata)
{

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);
        amrex::ParallelFor(nodalbox,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                           {
                               if (nodal_data_arr(i, j, k, MASS_SPHEAT) >
                                   shunya)
                               {
                                   nodal_data_arr(i, j, k, DELTA_TEMPERATURE) =
                                       nodal_data_arr(i, j, k, TEMPERATURE);
                               }
                           });
    }
}
#endif

#if USE_EB
/**
 * @brief Applies level-set momentum BCs to nodal velocities.
 * Per-body parameters (mom_bc_type, wall_mu, wall_vel) are stored in each
 * LevelSetBody and set by init_eb() from the input file.
 *
 * @param[in,out] nodaldata  Nodal MultiFab containing velocity fields.
 * @param[in]     geom       Coarse-level geometry.
 * @param[in]     dt         Time step (unused; reserved for future use).
 *
 * @return None.
 */
void nodal_levelset_bcs(MultiFab &nodaldata,
                        const Geometry &geom,
                        amrex::Real & /*dt*/)
{
    const auto plo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

    const int num_bodies = static_cast<int>(mpm_ebtools::ls_bodies.size());

    for (int b = 0; b < num_bodies; ++b)
    {
        const LevelSetBody &body = mpm_ebtools::ls_bodies[b];
        const int lsref = body.ls_refinement;

        const int bc_int = body.mom_bc_int();
        const amrex::Real wmu = body.wall_mu;
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wvel = body.wall_vel;

        MultiFab lsphi_coarse(nodaldata.boxArray(), nodaldata.DistributionMap(),
                              1,  // ncomp
                              1); // nghost — must be >= 1 for interpolation
        lsphi_coarse.setVal(1.0);
        amrex::average_down_nodal(*body.lsphi, lsphi_coarse,
                                  amrex::IntVect(lsref));
        lsphi_coarse.FillBoundary(geom.periodicity());

        for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
        {
            Box nodalbox =
                convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

            Array4<Real> nodal_data_arr = nodaldata.array(mfi);
            Array4<Real> lsarr = lsphi_coarse.array(mfi);

            amrex::ParallelFor(
                nodalbox,
                [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
                {
                    IntVect nodeid(AMREX_D_DECL(i, j, k));

                    amrex::Real xp[AMREX_SPACEDIM] = {AMREX_D_DECL(
                        plo[XDIR] + i * dx[XDIR], plo[YDIR] + j * dx[YDIR],
                        plo[ZDIR] + k * dx[ZDIR])};

                    amrex::Real lsval =
                        get_levelset_value(lsarr, plo, dx, xp, /*lsref=*/1);

                    if (lsval >= 0.0 ||
                        nodal_data_arr(nodeid, MASS_INDEX) <= shunya)
                        return;

                    amrex::Real normaldir[AMREX_SPACEDIM] = {
                        AMREX_D_DECL(1.0, 0.0, 0.0)};
                    get_levelset_grad(lsarr, plo, dx, xp, /*lsref=*/1,
                                      normaldir);

                    amrex::Real gradmag = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        gradmag += normaldir[d] * normaldir[d];
                    gradmag = std::sqrt(gradmag);

                    if (!(gradmag >= 1.0e-10))
                        return;

                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        normaldir[d] /= gradmag;

                    amrex::Real relvel_in[AMREX_SPACEDIM];
                    amrex::Real relvel_out[AMREX_SPACEDIM];
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                    {
                        relvel_in[d] =
                            nodal_data_arr(nodeid, VELX_INDEX + d) - wvel[d];
                        relvel_out[d] = relvel_in[d];
                    }

                    amrex::Real veln = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        veln += relvel_in[d] * normaldir[d];

                    if (veln > 0.0)
                        return;

                    applybc(relvel_in, relvel_out, wmu, normaldir, bc_int);

                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        nodal_data_arr(nodeid, VELX_INDEX + d) =
                            relvel_out[d] + wvel[d];
                });
        }
    }
}

#if USE_TEMP
/**
 * @brief Applies level-set temperature BCs to nodal temperatures.
 * When dirichlet_only is true only "isothermal" bodies are processed;
 * flux and convection bodies are skipped (matches the predictor-pass
 * convention used by Apply_Nodal_BCs_Temperature).
 *
 * @param[in,out] nodaldata      Nodal MultiFab containing temperature fields.
 * @param[in]     geom           Coarse-level geometry.
 * @param[in]     dirichlet_only If true, skip non-Dirichlet body types.
 *
 * @return None.
 */
void nodal_levelset_bcs_temperature(MultiFab &nodaldata,
                                    const Geometry &geom,
                                    bool dirichlet_only)
{
    const auto plo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

    const int num_bodies = static_cast<int>(mpm_ebtools::ls_bodies.size());

    for (int b = 0; b < num_bodies; ++b)
    {
        const LevelSetBody &body = mpm_ebtools::ls_bodies[b];

        const int bc_int = body.temp_bc_int();

        if (bc_int == 0)
            continue;

        if (dirichlet_only && bc_int != 1)
            continue;

        const int lsref = body.ls_refinement;

        const amrex::Real T_wall_v = body.T_wall;
        const amrex::Real heat_flux_v = body.heat_flux;
        const amrex::Real h_conv_v = body.h_conv;
        const amrex::Real T_inf_v = body.T_inf;

        MultiFab lsphi_coarse(nodaldata.boxArray(), nodaldata.DistributionMap(),
                              1, 1);
        amrex::average_down_nodal(*body.lsphi, lsphi_coarse,
                                  amrex::IntVect(lsref));
        lsphi_coarse.FillBoundary(geom.periodicity());

        for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
        {
            Box nodalbox =
                convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

            Array4<Real> arr = nodaldata.array(mfi);
            Array4<Real> lsarr = lsphi_coarse.array(mfi);

            const Box lsbox = lsphi_coarse[mfi].box();

            amrex::ParallelFor(
                nodalbox,
                [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
                {
                    IntVect nodeid(AMREX_D_DECL(i, j, k));

                    amrex::Real xp[AMREX_SPACEDIM] = {AMREX_D_DECL(
                        plo[XDIR] + i * dx[XDIR], plo[YDIR] + j * dx[YDIR],
                        plo[ZDIR] + k * dx[ZDIR])};

                    amrex::Real lsval =
                        get_levelset_value(lsarr, plo, dx, xp, /*lsref=*/1);

                    if (lsval >= 0.0 || arr(nodeid, MASS_SPHEAT) <= shunya)
                        return;

                    amrex::Real normaldir[AMREX_SPACEDIM] = {
                        AMREX_D_DECL(1.0, 0.0, 0.0)};
                    get_levelset_grad(lsarr, plo, dx, xp, /*lsref=*/1,
                                      normaldir);

                    amrex::Real gradmag = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        gradmag += normaldir[d] * normaldir[d];
                    gradmag = std::sqrt(gradmag);

                    if (gradmag < 1.0e-10)
                        return;

                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                        normaldir[d] /= gradmag;

                    if (bc_int == 1)
                    {
                        arr(nodeid, TEMPERATURE) = T_wall_v;
                        return;
                    }

                    int dom_dir = 0;
                    for (int d = 1; d < AMREX_SPACEDIM; d++)
                        if (std::abs(normaldir[d]) >
                            std::abs(normaldir[dom_dir]))
                            dom_dir = d;

                    IntVect nb = nodeid;
                    nb[dom_dir] += (normaldir[dom_dir] > 0.0) ? 1 : -1;

                    if (!lsbox.contains(nb))
                        return;

                    amrex::Real mk = arr(nodeid, MASS_CONDUCTIVITY);
                    amrex::Real m = arr(nodeid, MASS_INDEX);
                    amrex::Real k_node = (m > shunya) ? mk / m : eka;

                    amrex::Real T_nb = arr(nb, TEMPERATURE);

                    if (bc_int == 3)
                    {
                        arr(nodeid, TEMPERATURE) =
                            T_nb + heat_flux_v * dx[dom_dir] / k_node;
                    }
                    else if (bc_int == 4)
                    {
                        amrex::Real Bi = h_conv_v * dx[dom_dir] / k_node;
                        arr(nodeid, TEMPERATURE) =
                            (T_nb + Bi * T_inf_v) / (eka + Bi);
                    }
                });
        }
    }
}
#endif
#endif

/**
 * @brief Computes nodal Δv = v_new − v_old for PIC/FLIP blending.
 *
 * Requires that backup_current_velocity() has already stored v_old in
 * DELTA_VELX_INDEX. This routine overwrites DELTA_VELX_INDEX with:
 *
 *      Δv[d] = v[d] − v_old[d]
 *
 * Only nodes with positive mass are updated.
 *
 * @param[in,out] nodaldata  Nodal MultiFab containing velocity fields.
 *
 * @return None.
 */

void store_delta_velocity(MultiFab &nodaldata)
{

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (nodal_data_arr(i, j, k, MASS_INDEX) > shunya)
                {
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                    {
                        nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d) =
                            nodal_data_arr(i, j, k, VELX_INDEX + d) -
                            nodal_data_arr(i, j, k, DELTA_VELX_INDEX + d);
                    }
                }
            });
    }
}

#if USE_TEMP
/**
 * @brief Computes nodal ΔT = T_new − T_old for thermal PIC/FLIP updates.
 *
 * Requires that backup_current_temperature() has already stored T_old in
 * DELTA_TEMPERATURE. This routine overwrites DELTA_TEMPERATURE with:
 *
 *      ΔT = T − T_old
 *
 * Only nodes with positive MASS_SPHEAT are updated.
 *
 * @param[in,out] nodaldata  Nodal MultiFab containing thermal fields.
 *
 * @return None.
 */
void store_delta_temperature(MultiFab &nodaldata)
{

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (nodal_data_arr(i, j, k, MASS_SPHEAT) > shunya)
                {
                    nodal_data_arr(i, j, k, DELTA_TEMPERATURE) =
                        nodal_data_arr(i, j, k, TEMPERATURE) -
                        nodal_data_arr(i, j, k, DELTA_TEMPERATURE);
                }
            });
    }
}
#endif

/**
 * @brief Advances nodal velocities using nodal forces (explicit time
 * integration).
 *
 * For each node with mass ≥ mass_tolerance:
 *
 *      v[d] ← v[d] + (F[d] / m) * dt
 *
 * Otherwise, nodal velocity is zeroed.
 *
 * @param[in,out] nodaldata      Nodal MultiFab containing mass, force,
 * velocity.
 * @param[in]     dt             Time step.
 * @param[in]     mass_tolerance Minimum mass required to update velocity.
 *
 * @return None.
 */

void Nodal_Time_Update_Momentum(MultiFab &nodaldata,
                                const amrex::Real &dt,
                                const amrex::Real &mass_tolerance)
{
    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (nodal_data_arr(i, j, k, MASS_INDEX) >= mass_tolerance)
                {
                    for (int d = 0; d < AMREX_SPACEDIM; d++)
                    {
                        nodal_data_arr(i, j, k, VELX_INDEX + d) +=
                            nodal_data_arr(i, j, k, FRCX_INDEX + d) /
                            nodal_data_arr(i, j, k, MASS_INDEX) * dt;
                    }
                }
                else
                {
                    nodal_data_arr(i, j, k, VELX_INDEX) = 0.0;
                    nodal_data_arr(i, j, k, VELY_INDEX) = 0.0;
                    nodal_data_arr(i, j, k, VELZ_INDEX) = 0.0;
                }
            });
    }
}

#if USE_TEMP
/**
 * @brief Advances nodal temperature using nodal heat sources.
 *
 * For each node with MASS_SPHEAT ≥ mass_tolerance:
 *
 *      T ← T + (source / MASS_SPHEAT) * dt
 *
 * Otherwise, temperature is set to zero.
 *
 * @param[in,out] nodaldata      Nodal MultiFab containing thermal fields.
 * @param[in]     dt             Time step.
 * @param[in]     mass_tolerance Minimum thermal mass required to update T.
 *
 * @return None.
 */
void Nodal_Time_Update_Temperature(MultiFab &nodaldata,
                                   const amrex::Real &dt,
                                   const amrex::Real &mass_tolerance)
{

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (nodal_data_arr(i, j, k, MASS_SPHEAT) >= mass_tolerance)
                {
                    nodal_data_arr(i, j, k, TEMPERATURE) +=
                        nodal_data_arr(i, j, k, SOURCE_TEMP_INDEX) /
                        nodal_data_arr(i, j, k, MASS_SPHEAT) * dt;
                }
                else
                {
                    nodal_data_arr(i, j, k, TEMPERATURE) = 0.0;
                }
            });
    }
}
#endif

/**
 * @brief Precomputes per‑node shape‑function index categories for boundary
 * handling.
 *
 * For each node, assigns an integer code (0–4) in each dimension indicating:
 *   - 0: domain lower boundary
 *   - 1: near lower boundary
 *   - 2: interior
 *   - 3: near upper boundary
 *   - 4: domain upper boundary
 *
 * These indices are used to select one‑sided or symmetric interpolation
 * stencils.
 *
 * @param[out] shapefunctionindex  iMultiFab storing index codes.
 * @param[in]  geom                Geometry describing the domain.
 *
 * @return None.
 */

void initialise_shape_function_indices(iMultiFab &shapefunctionindex,
                                       const amrex::Geometry geom)
{
    const int *domloarr = geom.Domain().loVect();
    const int *domhiarr = geom.Domain().hiVect();

    GpuArray<int, AMREX_SPACEDIM> lo = {
        AMREX_D_DECL(domloarr[0], domloarr[1], domloarr[2])};
    GpuArray<int, AMREX_SPACEDIM> hi = {
        AMREX_D_DECL(domhiarr[0], domhiarr[1], domhiarr[2])};

    for (MFIter mfi(shapefunctionindex); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<int> shapefunctionindex_arr = shapefunctionindex.array(mfi);

        amrex::ParallelFor(nodalbox,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                           {
                               if (i == lo[0])
                               {
                                   shapefunctionindex_arr(i, j, k, 0) = 0;
                               }
                               else if (i == lo[0] + 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 0) = 1;
                               }
                               else if (i == hi[0])
                               {
                                   shapefunctionindex_arr(i, j, k, 0) = 4;
                               }
                               else if (i == hi[0] - 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 0) = 3;
                               }
                               else
                               {
                                   shapefunctionindex_arr(i, j, k, 0) = 2;
                               }

                               if (j == lo[1])
                               {
                                   shapefunctionindex_arr(i, j, k, 1) = 0;
                               }
                               else if (j == lo[1] + 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 1) = 1;
                               }
                               else if (j == hi[1])
                               {
                                   shapefunctionindex_arr(i, j, k, 1) = 4;
                               }
                               else if (j == hi[1] - 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 1) = 3;
                               }
                               else
                               {
                                   shapefunctionindex_arr(i, j, k, 1) = 2;
                               }

                               if (k == lo[2])
                               {
                                   shapefunctionindex_arr(i, j, k, 2) = 0;
                               }
                               else if (k == lo[2] + 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 2) = 1;
                               }
                               else if (k == hi[2])
                               {
                                   shapefunctionindex_arr(i, j, k, 2) = 4;
                               }
                               else if (k == hi[2] - 1)
                               {
                                   shapefunctionindex_arr(i, j, k, 2) = 3;
                               }
                               else
                               {
                                   shapefunctionindex_arr(i, j, k, 2) = 2;
                               }
                           });
    }
}

/**
 * @brief Applies velocity boundary conditions at nodal locations.
 *
 * For each node lying on a domain boundary:
 *   1. Subtracts wall velocity to obtain relative velocity.
 *   2. Applies applybc() using the appropriate wall friction μ and BC type.
 *   3. Restores wall velocity to obtain the final nodal velocity.
 *
 * Supports:
 *   - No‑slip
 *   - Slip
 *   - Partial slip
 *   - Periodic
 *
 * @param[in]     geom             Geometry describing the domain.
 * @param[in,out] nodaldata        Nodal MultiFab containing velocity fields.
 * @param[in]     bcloarr          BC types at lower faces.
 * @param[in]     bchiarr          BC types at upper faces.
 * @param[in]     wall_mu_loarr    Friction coefficients at lower faces.
 * @param[in]     wall_mu_hiarr    Friction coefficients at upper faces.
 * @param[in]     wall_vel_loarr   Wall velocities at lower faces.
 * @param[in]     wall_vel_hiarr   Wall velocities at upper faces.
 * @param[in]     dt               Time step (unused).
 *
 * @return None.
 */

void nodal_bcs(const amrex::Geometry geom,
               MultiFab &nodaldata,
               int bcloarr[AMREX_SPACEDIM],
               int bchiarr[AMREX_SPACEDIM],
               amrex::Real wall_mu_loarr[AMREX_SPACEDIM],
               amrex::Real wall_mu_hiarr[AMREX_SPACEDIM],
               amrex::Real wall_vel_loarr[AMREX_SPACEDIM * AMREX_SPACEDIM],
               amrex::Real wall_vel_hiarr[AMREX_SPACEDIM * AMREX_SPACEDIM],
               const amrex::Real & /*dt*/)
{
    const int *domloarr = geom.Domain().loVect();
    const int *domhiarr = geom.Domain().hiVect();

    GpuArray<int, AMREX_SPACEDIM> domlo;
    GpuArray<int, AMREX_SPACEDIM> domhi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        domlo[d] = domloarr[d];
        domhi[d] = domhiarr[d];
    }

    GpuArray<int, AMREX_SPACEDIM> bclo;
    GpuArray<int, AMREX_SPACEDIM> bchi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        bclo[d] = bcloarr[d];
        bchi[d] = bchiarr[d];
    }

    GpuArray<amrex::Real, AMREX_SPACEDIM> wall_mu_lo;
    GpuArray<amrex::Real, AMREX_SPACEDIM> wall_mu_hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        wall_mu_lo[d] = wall_mu_loarr[d];
        wall_mu_hi[d] = wall_mu_hiarr[d];
    }

    GpuArray<amrex::Real, AMREX_SPACEDIM * AMREX_SPACEDIM> wall_vel_lo;
    GpuArray<amrex::Real, AMREX_SPACEDIM * AMREX_SPACEDIM> wall_vel_hi;
    for (int d = 0; d < AMREX_SPACEDIM * AMREX_SPACEDIM; ++d)
    {
        wall_vel_lo[d] = wall_vel_loarr[d];
        wall_vel_hi[d] = wall_vel_hiarr[d];
    }

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));

        Array4<amrex::Real> nodal_data_arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeid(AMREX_D_DECL(i, j, k));

                amrex::Real relvel_in[AMREX_SPACEDIM],
                    relvel_out[AMREX_SPACEDIM];
                amrex::Real wallvel[AMREX_SPACEDIM] = {0.0};

                // Initialize relvel arrays
                for (int d = 0; d < AMREX_SPACEDIM; d++)
                {
                    relvel_in[d] = nodal_data_arr(nodeid, VELX_INDEX + d);
                    relvel_out[d] = relvel_in[d];
                }

                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    if (nodeid[dir] == domlo[dir])
                    {
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        {
                            wallvel[d] = wall_vel_lo[dir * AMREX_SPACEDIM + d];
                            relvel_in[d] -= wallvel[d];
                        }
                        amrex::Real normaldir[AMREX_SPACEDIM] = {0.0};
                        normaldir[dir] = 1.0;
                        applybc(relvel_in, relvel_out, wall_mu_lo[dir],
                                normaldir, bclo[dir]);
                        for (int c = 0; c < AMREX_SPACEDIM; ++c)
                            nodal_data_arr(nodeid, VELX_INDEX + c) =
                                relvel_out[c] + wallvel[c];
                    }
                    else if (nodeid[dir] == domhi[dir] + 1)
                    {
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        {
                            wallvel[d] = wall_vel_hi[dir * AMREX_SPACEDIM + d];
                            relvel_in[d] -= wallvel[d];
                        }
                        amrex::Real normaldir[AMREX_SPACEDIM] = {0.0};
                        normaldir[dir] = -1.0;
                        applybc(relvel_in, relvel_out, wall_mu_hi[dir],
                                normaldir, bchi[dir]);
                        for (int c = 0; c < AMREX_SPACEDIM; ++c)
                            nodal_data_arr(nodeid, VELX_INDEX + c) =
                                relvel_out[c] + wallvel[c];
                    }
                }
            });
    }
}

#if USE_TEMP
/**
 * @brief Applies Dirichlet temperature boundary conditions at nodal locations.
 *
 * For each dimension d:
 *   - At nodeid[d] = domlo[d], set T = wall_temp_lo[d]
 *   - At nodeid[d] = domhi[d] + 1, set T = wall_temp_hi[d]
 *
 * @param[in]     geom                   Geometry describing the domain.
 * @param[in,out] nodaldata              Nodal MultiFab containing temperature.
 * @param[in]     bcloarr                BC types at lower faces (unused).
 * @param[in]     bchiarr                BC types at upper faces (unused).
 * @param[in]     dirichlet_temperature_lo  Temperature at lower faces.
 * @param[in]     dirichlet_temperature_hi  Temperature at upper faces.
 *
 * @return None.
 */

void nodal_bcs_temperature(const amrex::Geometry geom,
                           MultiFab &nodaldata,
                           int bcloarr[AMREX_SPACEDIM],
                           int bchiarr[AMREX_SPACEDIM],
                           amrex::Real T_wall_lo[AMREX_SPACEDIM],
                           amrex::Real T_wall_hi[AMREX_SPACEDIM],
                           amrex::Real flux_lo[AMREX_SPACEDIM],
                           amrex::Real flux_hi[AMREX_SPACEDIM],
                           amrex::Real h_lo[AMREX_SPACEDIM],
                           amrex::Real h_hi[AMREX_SPACEDIM],
                           amrex::Real Tinf_lo[AMREX_SPACEDIM],
                           amrex::Real Tinf_hi[AMREX_SPACEDIM])
{
    const auto dx = geom.CellSize();
    const int *domloarr = geom.Domain().loVect();
    const int *domhiarr = geom.Domain().hiVect();

    GpuArray<int, AMREX_SPACEDIM> domlo;
    GpuArray<int, AMREX_SPACEDIM> domhi;
    GpuArray<int, AMREX_SPACEDIM> bclo;
    GpuArray<int, AMREX_SPACEDIM> bchi;
    GpuArray<amrex::Real, AMREX_SPACEDIM> T_wall_lo_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> T_wall_hi_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> flux_lo_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> flux_hi_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> h_lo_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> h_hi_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> Tinf_lo_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> Tinf_hi_g;
    GpuArray<amrex::Real, AMREX_SPACEDIM> dx_g;

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        domlo[d] = domloarr[d];
        domhi[d] = domhiarr[d];
        bclo[d] = bcloarr[d];
        bchi[d] = bchiarr[d];
        T_wall_lo_g[d] = T_wall_lo[d];
        T_wall_hi_g[d] = T_wall_hi[d];
        flux_lo_g[d] = flux_lo[d];
        flux_hi_g[d] = flux_hi[d];
        h_lo_g[d] = h_lo[d];
        h_hi_g[d] = h_hi[d];
        Tinf_lo_g[d] = Tinf_lo[d];
        Tinf_hi_g[d] = Tinf_hi[d];
        dx_g[d] = dx[d];
    }

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        Array4<amrex::Real> arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                (void)j;
                (void)k;
                IntVect nodeid(AMREX_D_DECL(i, j, k));

                bool bc_applied = false;
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    bool is_lo = (nodeid[d] == domlo[d]);
                    bool is_hi = (nodeid[d] == domhi[d] + 1);

                    if (is_lo || is_hi)
                    {
                        int bc_type = is_lo ? bclo[d] : bchi[d];
                        int sign = is_lo ? 1 : -1;

                        if (bc_type == BC_TEMP_PERIODIC)
                        {
                            // Periodic. Do nothing
                            bc_applied = true;
                        }
                        else if (bc_type == BC_TEMP_ISOTHERMAL)
                        {
                            amrex::Real Tw =
                                is_lo ? T_wall_lo_g[d] : T_wall_hi_g[d];
                            arr(nodeid, TEMPERATURE) = Tw;
                            bc_applied = true;
                        }
                        else if (bc_type == BC_TEMP_ADIABATIC)
                        {
                            if (!bc_applied)
                            {
                                IntVect nb = nodeid;
                                nb[d] += sign;
                                arr(nodeid, TEMPERATURE) = arr(nb, TEMPERATURE);
                                bc_applied = true;
                            }
                        }
                        else if (bc_type == BC_TEMP_USERFLUX)
                        {
                            IntVect nb = nodeid;
                            nb[d] += sign;
                            amrex::Real q = is_lo ? flux_lo_g[d] : flux_hi_g[d];
                            amrex::Real mk = arr(nodeid, MASS_CONDUCTIVITY);
                            amrex::Real m = arr(nodeid, MASS_INDEX);
                            amrex::Real k_node = (m > shunya) ? mk / m : eka;
                            arr(nodeid, TEMPERATURE) =
                                arr(nb, TEMPERATURE) + q * dx_g[d] / k_node;
                            bc_applied = true;
                        }
                        else if (bc_type == BC_TEMP_CONVECTION)
                        {
                            IntVect nb = nodeid;
                            nb[d] += sign;
                            amrex::Real hc = is_lo ? h_lo_g[d] : h_hi_g[d];
                            amrex::Real Tinf =
                                is_lo ? Tinf_lo_g[d] : Tinf_hi_g[d];
                            amrex::Real Bi = hc * dx_g[d];
                            arr(nodeid, TEMPERATURE) =
                                (arr(nb, TEMPERATURE) + Bi * Tinf) / (1.0 + Bi);
                            bc_applied = true;
                        }
                    }
                }
            });
    }
}
#endif

/**
 * @brief Sets all nodal data components to zero (shunya).
 *
 * Resets the MultiFab including ghost cells.
 *
 * @param[in,out] nodaldata            Nodal MultiFab to reset.
 * @param[in]     ng_cells_nodaldata   Number of ghost cells to include.
 *
 * @return None.
 */

void Reset_Nodaldata_to_Zero(amrex::MultiFab &nodaldata, int ng_cells_nodaldata)
{
    nodaldata.setVal(shunya, ng_cells_nodaldata);
}

void compute_udf_wall_vel_at_nodes(const amrex::Geometry &geom,
                                   MPMspecs &specs,
                                   amrex::Real t)
{
#if (AMREX_SPACEDIM >= 2)
    const auto dx = geom.CellSizeArray();
#endif
    const auto plo_g = geom.ProbLoArray();
    const auto phi_g = geom.ProbHiArray();

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
#if (AMREX_SPACEDIM >= 2)
        int p0 = (d == 0) ? 1 : 0;
#endif
#if (AMREX_SPACEDIM == 3)
        int p1 = (d == 2) ? 1 : 2;
#endif

        if (specs.udf_func_ptr_lo[d] != nullptr)
        {
#if (AMREX_SPACEDIM == 1)
            int n_nodes = 1;
#elif (AMREX_SPACEDIM == 2)
            int n_nodes = specs.ncells[p0] + 1;
#else
            int n_nodes = (specs.ncells[p0] + 1) * (specs.ncells[p1] + 1);
#endif
            amrex::Vector<amrex::Real> hv(n_nodes * AMREX_SPACEDIM, 0.0);

#if (AMREX_SPACEDIM == 1)
            {
                double vel[3] = {0.0, 0.0, 0.0};
                specs.udf_func_ptr_lo[d]((double)plo_g[0], 0.0, 0.0, (double)t,
                                         vel);
                for (int c = 0; c < AMREX_SPACEDIM; ++c)
                    hv[c] = (amrex::Real)vel[c];
            }
#elif (AMREX_SPACEDIM == 2)
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                double cx = (d == 0) ? (double)plo_g[0]
                                     : (double)(plo_g[0] + j * dx[0]);
                double cy = (d == 1) ? (double)plo_g[1]
                                     : (double)(plo_g[1] + j * dx[1]);
                double vel[3] = {0.0, 0.0, 0.0};
                specs.udf_func_ptr_lo[d](cx, cy, 0.0, (double)t, vel);
                for (int c = 0; c < AMREX_SPACEDIM; ++c)
                    hv[j * AMREX_SPACEDIM + c] = (amrex::Real)vel[c];
            }
#else
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                for (int k = 0; k <= specs.ncells[p1]; ++k)
                {
                    double coords[3] = {0.0, 0.0, 0.0};
                    coords[d] = (double)plo_g[d];
                    coords[p0] = (double)(plo_g[p0] + j * dx[p0]);
                    coords[p1] = (double)(plo_g[p1] + k * dx[p1]);
                    double vel[3] = {0.0, 0.0, 0.0};
                    specs.udf_func_ptr_lo[d](coords[0], coords[1], coords[2],
                                             (double)t, vel);
                    int flat =
                        (j * (specs.ncells[p1] + 1) + k) * AMREX_SPACEDIM;
                    for (int c = 0; c < AMREX_SPACEDIM; ++c)
                        hv[flat + c] = (amrex::Real)vel[c];
                }
            }
#endif
            specs.udf_wall_vel_nodes_lo[d].resize(hv.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, hv.begin(), hv.end(),
                             specs.udf_wall_vel_nodes_lo[d].begin());
        }

        if (specs.udf_func_ptr_hi[d] != nullptr)
        {
#if (AMREX_SPACEDIM == 1)
            int n_nodes = 1;
#elif (AMREX_SPACEDIM == 2)
            int n_nodes = specs.ncells[p0] + 1;
#else
            int n_nodes = (specs.ncells[p0] + 1) * (specs.ncells[p1] + 1);
#endif
            amrex::Vector<amrex::Real> hv(n_nodes * AMREX_SPACEDIM, 0.0);

#if (AMREX_SPACEDIM == 1)
            {
                double vel[3] = {0.0, 0.0, 0.0};
                specs.udf_func_ptr_hi[d]((double)phi_g[0], 0.0, 0.0, (double)t,
                                         vel);
                for (int c = 0; c < AMREX_SPACEDIM; ++c)
                    hv[c] = (amrex::Real)vel[c];
            }
#elif (AMREX_SPACEDIM == 2)
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                double cx = (d == 0) ? (double)phi_g[0]
                                     : (double)(plo_g[0] + j * dx[0]);
                double cy = (d == 1) ? (double)phi_g[1]
                                     : (double)(plo_g[1] + j * dx[1]);
                double vel[3] = {0.0, 0.0, 0.0};
                specs.udf_func_ptr_hi[d](cx, cy, 0.0, (double)t, vel);
                for (int c = 0; c < AMREX_SPACEDIM; ++c)
                    hv[j * AMREX_SPACEDIM + c] = (amrex::Real)vel[c];
            }
#else
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                for (int k = 0; k <= specs.ncells[p1]; ++k)
                {
                    double coords[3] = {0.0, 0.0, 0.0};
                    coords[d] = (double)phi_g[d];
                    coords[p0] = (double)(plo_g[p0] + j * dx[p0]);
                    coords[p1] = (double)(plo_g[p1] + k * dx[p1]);
                    double vel[3] = {0.0, 0.0, 0.0};
                    specs.udf_func_ptr_hi[d](coords[0], coords[1], coords[2],
                                             (double)t, vel);
                    int flat =
                        (j * (specs.ncells[p1] + 1) + k) * AMREX_SPACEDIM;
                    for (int c = 0; c < AMREX_SPACEDIM; ++c)
                        hv[flat + c] = (amrex::Real)vel[c];
                }
            }
#endif
            specs.udf_wall_vel_nodes_hi[d].resize(hv.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, hv.begin(), hv.end(),
                             specs.udf_wall_vel_nodes_hi[d].begin());
        }
    }
    amrex::Gpu::streamSynchronize();
}

void apply_udf_nodal_bcs(const amrex::Geometry &geom,
                         amrex::MultiFab &nodaldata,
                         MPMspecs &specs)
{
    bool any_udf = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        if (specs.udf_func_ptr_lo[d] || specs.udf_func_ptr_hi[d])
        {
            any_udf = true;
            break;
        }
    }
    if (!any_udf)
        return;

    const int *domloarr = geom.Domain().loVect();
    const int *domhiarr = geom.Domain().hiVect();

    amrex::GpuArray<int, AMREX_SPACEDIM> domlo, domhi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        domlo[d] = domloarr[d];
        domhi[d] = domhiarr[d];
    }

    amrex::GpuArray<int, AMREX_SPACEDIM> bclo_arr, bchi_arr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        bclo_arr[d] = specs.bclo[d];
        bchi_arr[d] = specs.bchi[d];
    }

    amrex::GpuArray<const amrex::Real *, AMREX_SPACEDIM> udf_lo_ptrs,
        udf_hi_ptrs;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        udf_lo_ptrs[d] = specs.udf_func_ptr_lo[d]
                             ? specs.udf_wall_vel_nodes_lo[d].data()
                             : nullptr;
        udf_hi_ptrs[d] = specs.udf_func_ptr_hi[d]
                             ? specs.udf_wall_vel_nodes_hi[d].data()
                             : nullptr;
    }

    amrex::GpuArray<int, AMREX_SPACEDIM> ncells_arr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        ncells_arr[d] = specs.ncells[d];

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        Array4<amrex::Real> arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeid(AMREX_D_DECL(i, j, k));

                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    const amrex::Real *udf_ptr = nullptr;
                    int bc_type = -1;

                    if (nodeid[dir] == domlo[dir] &&
                        udf_lo_ptrs[dir] != nullptr)
                    {
                        udf_ptr = udf_lo_ptrs[dir];
                        bc_type = bclo_arr[dir];
                    }
                    else if (nodeid[dir] == domhi[dir] + 1 &&
                             udf_hi_ptrs[dir] != nullptr)
                    {
                        udf_ptr = udf_hi_ptrs[dir];
                        bc_type = bchi_arr[dir];
                    }

                    if (udf_ptr == nullptr)
                        continue;

#if (AMREX_SPACEDIM >= 2)
                    int p0 = (dir == 0) ? 1 : 0;
#endif
#if (AMREX_SPACEDIM == 3)
                    int p1 = (dir == 2) ? 1 : 2;
#endif

                    int flat = 0;
#if (AMREX_SPACEDIM == 1)
                    flat = 0;
#elif (AMREX_SPACEDIM == 2)
                    flat = nodeid[p0] * AMREX_SPACEDIM;
#else
                    flat = (nodeid[p0] * (ncells_arr[p1] + 1) + nodeid[p1]) *
                           AMREX_SPACEDIM;
#endif

                    amrex::Real W[AMREX_SPACEDIM];
                    for (int c = 0; c < AMREX_SPACEDIM; ++c)
                        W[c] = udf_ptr[flat + c];

                    if (bc_type == BC_NOSLIPWALL)
                    {
                        for (int c = 0; c < AMREX_SPACEDIM; ++c)
                            arr(nodeid, VELX_INDEX + c) = W[c];
                    }
                    else
                    {
                        arr(nodeid, VELX_INDEX + dir) = W[dir];
                    }
                }
            });
    }
}

#if USE_TEMP
void compute_udf_temp_at_nodes(const amrex::Geometry &geom,
                               MPMspecs &specs,
                               amrex::Real t)
{
#if (AMREX_SPACEDIM >= 2)
    const auto dx = geom.CellSizeArray();
#endif
    const auto plo_g = geom.ProbLoArray();
    const auto phi_g = geom.ProbHiArray();

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
#if (AMREX_SPACEDIM >= 2)
        int p0 = (d == 0) ? 1 : 0;
#endif
#if (AMREX_SPACEDIM == 3)
        int p1 = (d == 2) ? 1 : 2;
#endif

        if (specs.udf_temp_func_ptr_lo[d] != nullptr)
        {
#if (AMREX_SPACEDIM == 1)
            int n_nodes = 1;
#elif (AMREX_SPACEDIM == 2)
            int n_nodes = specs.ncells[p0] + 1;
#else
            int n_nodes = (specs.ncells[p0] + 1) * (specs.ncells[p1] + 1);
#endif
            amrex::Vector<amrex::Real> hv(n_nodes * 2, 0.0);

#if (AMREX_SPACEDIM == 1)
            {
                double out[2] = {0.0, 0.0};
                specs.udf_temp_func_ptr_lo[d]((double)plo_g[0], 0.0, 0.0,
                                              (double)t, out);
                hv[0] = (amrex::Real)out[0];
                hv[1] = (amrex::Real)out[1];
            }
#elif (AMREX_SPACEDIM == 2)
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                double cx = (d == 0) ? (double)plo_g[0]
                                     : (double)(plo_g[0] + j * dx[0]);
                double cy = (d == 1) ? (double)plo_g[1]
                                     : (double)(plo_g[1] + j * dx[1]);
                double out[2] = {0.0, 0.0};
                specs.udf_temp_func_ptr_lo[d](cx, cy, 0.0, (double)t, out);
                hv[j * 2] = (amrex::Real)out[0];
                hv[j * 2 + 1] = (amrex::Real)out[1];
            }
#else
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                for (int k = 0; k <= specs.ncells[p1]; ++k)
                {
                    double coords[3] = {0.0, 0.0, 0.0};
                    coords[d] = (double)plo_g[d];
                    coords[p0] = (double)(plo_g[p0] + j * dx[p0]);
                    coords[p1] = (double)(plo_g[p1] + k * dx[p1]);
                    double out[2] = {0.0, 0.0};
                    specs.udf_temp_func_ptr_lo[d](coords[0], coords[1],
                                                  coords[2], (double)t, out);
                    int flat = (j * (specs.ncells[p1] + 1) + k) * 2;
                    hv[flat] = (amrex::Real)out[0];
                    hv[flat + 1] = (amrex::Real)out[1];
                }
            }
#endif
            specs.udf_temp_nodes_lo[d].resize(hv.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, hv.begin(), hv.end(),
                             specs.udf_temp_nodes_lo[d].begin());
        }

        if (specs.udf_temp_func_ptr_hi[d] != nullptr)
        {
#if (AMREX_SPACEDIM == 1)
            int n_nodes = 1;
#elif (AMREX_SPACEDIM == 2)
            int n_nodes = specs.ncells[p0] + 1;
#else
            int n_nodes = (specs.ncells[p0] + 1) * (specs.ncells[p1] + 1);
#endif
            amrex::Vector<amrex::Real> hv(n_nodes * 2, 0.0);

#if (AMREX_SPACEDIM == 1)
            {
                double out[2] = {0.0, 0.0};
                specs.udf_temp_func_ptr_hi[d]((double)phi_g[0], 0.0, 0.0,
                                              (double)t, out);
                hv[0] = (amrex::Real)out[0];
                hv[1] = (amrex::Real)out[1];
            }
#elif (AMREX_SPACEDIM == 2)
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                double cx = (d == 0) ? (double)phi_g[0]
                                     : (double)(plo_g[0] + j * dx[0]);
                double cy = (d == 1) ? (double)phi_g[1]
                                     : (double)(plo_g[1] + j * dx[1]);
                double out[2] = {0.0, 0.0};
                specs.udf_temp_func_ptr_hi[d](cx, cy, 0.0, (double)t, out);
                hv[j * 2] = (amrex::Real)out[0];
                hv[j * 2 + 1] = (amrex::Real)out[1];
            }
#else
            for (int j = 0; j <= specs.ncells[p0]; ++j)
            {
                for (int k = 0; k <= specs.ncells[p1]; ++k)
                {
                    double coords[3] = {0.0, 0.0, 0.0};
                    coords[d] = (double)phi_g[d];
                    coords[p0] = (double)(plo_g[p0] + j * dx[p0]);
                    coords[p1] = (double)(plo_g[p1] + k * dx[p1]);
                    double out[2] = {0.0, 0.0};
                    specs.udf_temp_func_ptr_hi[d](coords[0], coords[1],
                                                  coords[2], (double)t, out);
                    int flat = (j * (specs.ncells[p1] + 1) + k) * 2;
                    hv[flat] = (amrex::Real)out[0];
                    hv[flat + 1] = (amrex::Real)out[1];
                }
            }
#endif
            specs.udf_temp_nodes_hi[d].resize(hv.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, hv.begin(), hv.end(),
                             specs.udf_temp_nodes_hi[d].begin());
        }
    }
    amrex::Gpu::streamSynchronize();
}

void apply_udf_nodal_bcs_temperature(const amrex::Geometry &geom,
                                     amrex::MultiFab &nodaldata,
                                     MPMspecs &specs)
{
    bool any_udf = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        if (specs.udf_temp_func_ptr_lo[d] || specs.udf_temp_func_ptr_hi[d])
        {
            any_udf = true;
            break;
        }
    }
    if (!any_udf)
        return;

    const auto dx = geom.CellSize();
    const int *domloarr = geom.Domain().loVect();
    const int *domhiarr = geom.Domain().hiVect();

    amrex::GpuArray<int, AMREX_SPACEDIM> domlo, domhi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        domlo[d] = domloarr[d];
        domhi[d] = domhiarr[d];
    }

    amrex::GpuArray<int, AMREX_SPACEDIM> bclo_arr, bchi_arr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        bclo_arr[d] = specs.bclo_temp[d];
        bchi_arr[d] = specs.bchi_temp[d];
    }

    amrex::GpuArray<const amrex::Real *, AMREX_SPACEDIM> udf_lo_ptrs,
        udf_hi_ptrs;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        udf_lo_ptrs[d] = specs.udf_temp_func_ptr_lo[d]
                             ? specs.udf_temp_nodes_lo[d].data()
                             : nullptr;
        udf_hi_ptrs[d] = specs.udf_temp_func_ptr_hi[d]
                             ? specs.udf_temp_nodes_hi[d].data()
                             : nullptr;
    }

    amrex::GpuArray<int, AMREX_SPACEDIM> ncells_arr;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        ncells_arr[d] = specs.ncells[d];

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_g;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        dx_g[d] = dx[d];

    for (MFIter mfi(nodaldata); mfi.isValid(); ++mfi)
    {
        Box nodalbox = convert(mfi.tilebox(), IntVect(AMREX_D_DECL(1, 1, 1)));
        Array4<amrex::Real> arr = nodaldata.array(mfi);

        amrex::ParallelFor(
            nodalbox,
            [=] AMREX_GPU_DEVICE(AMREX_D_DECL(int i, int j, int k)) noexcept
            {
                IntVect nodeid(AMREX_D_DECL(i, j, k));

                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    const amrex::Real *udf_ptr = nullptr;
                    int bc_type = -1;
                    int sign = 0;

                    if (nodeid[dir] == domlo[dir] &&
                        udf_lo_ptrs[dir] != nullptr)
                    {
                        udf_ptr = udf_lo_ptrs[dir];
                        bc_type = bclo_arr[dir];
                        sign = 1;
                    }
                    else if (nodeid[dir] == domhi[dir] + 1 &&
                             udf_hi_ptrs[dir] != nullptr)
                    {
                        udf_ptr = udf_hi_ptrs[dir];
                        bc_type = bchi_arr[dir];
                        sign = -1;
                    }

                    if (udf_ptr == nullptr)
                        continue;

#if (AMREX_SPACEDIM >= 2)
                    int p0 = (dir == 0) ? 1 : 0;
#endif
#if (AMREX_SPACEDIM == 3)
                    int p1 = (dir == 2) ? 1 : 2;
#endif

                    int flat = 0;
#if (AMREX_SPACEDIM == 1)
                    flat = 0;
#elif (AMREX_SPACEDIM == 2)
                    flat = nodeid[p0] * 2;
#else
                    flat = (nodeid[p0] * (ncells_arr[p1] + 1) + nodeid[p1]) * 2;
#endif

                    amrex::Real val0 = udf_ptr[flat];
                    amrex::Real val1 = udf_ptr[flat + 1];

                    if (bc_type == BC_TEMP_ISOTHERMAL)
                    {
                        arr(nodeid, TEMPERATURE) = val0;
                    }
                    else if (bc_type == BC_TEMP_USERFLUX)
                    {
                        IntVect nb = nodeid;
                        nb[dir] += sign;
                        amrex::Real mk = arr(nodeid, MASS_CONDUCTIVITY);
                        amrex::Real m = arr(nodeid, MASS_INDEX);
                        amrex::Real k_node = (m > shunya) ? mk / m : eka;
                        arr(nodeid, TEMPERATURE) =
                            arr(nb, TEMPERATURE) + val0 * dx_g[dir] / k_node;
                    }
                    else if (bc_type == BC_TEMP_CONVECTION)
                    {
                        IntVect nb = nodeid;
                        nb[dir] += sign;
                        amrex::Real hc = val0;
                        amrex::Real Tinf = val1;
                        amrex::Real Bi = hc * dx_g[dir];
                        arr(nodeid, TEMPERATURE) =
                            (arr(nb, TEMPERATURE) + Bi * Tinf) / (eka + Bi);
                    }
                }
            });
    }
}
#endif
