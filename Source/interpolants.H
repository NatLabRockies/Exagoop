#ifndef PARTICLEBASIS_H_
#define PARTICLEBASIS_H_
// clang-format off
#include<constants.H>
#include <AMReX.H>
#include<mpm_particle_container.H>
// clang-format on

struct CubicSplineCoeff
{
    amrex::Real a, b, c, d;
};

struct CubicSplineDerCoeff
{
    amrex::Real a, b, c; // f(zi) = a*zi^2 + b*zi + c
};

struct QuadraticSplineCoeff
{
    amrex::Real a, b, c;
};

struct QuadraticSplineDerCoeff
{
    amrex::Real a, b;
};

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_GPU_CONSTANT
#endif
    constexpr int interval_map_quadbspline[5][6] = {
        {-1, -1, -1, -1, -1, -1}, // type 0 filler
        {0, 0, 1, 2, 3, 3},       // type 1
        {0, 1, 1, 2, 2, -1},      // type 2
        {0, 0, 1, 1, 2, 2},       // type 3
        {0, 0, 1, 1, 2, -1}       // type 4
};
static constexpr int interval_map_quadbspline_host[5][6] = {
    {-1, -1, -1, -1, -1, -1}, // type 0 filler
    {0, 0, 1, 2, 3, 3},       // type 1
    {0, 1, 1, 2, 2, -1},      // type 2
    {0, 0, 1, 1, 2, 2},       // type 3
    {0, 0, 1, 1, 2, -1}       // type 4
};
#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_GPU_CONSTANT
#endif
    constexpr CubicSplineCoeff cubic_spline_coeffs[5][4] = {
        {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},

        // shapefunctiontype = 1 (Boundary node)
        {
            {onebysix, 1.0, 2.0, fourbythree},  // [-2,-1]
            {-onebysix, 0.0, 1.0, 1.0},         // [-1,0]
            {onebysix, 0.0, -1.0, 1.0},         // [0,1]
            {-onebysix, 1.0, -2.0, fourbythree} // [1,2]
        },

        // shapefunctiontype = 2 (Near boundary node)
        {{0.0, 0.0, 0.0, 0.0},                 // [-2,-1]
         {-onebythree, -1.0, 0.0, twobythree}, // [-1,0]
         {0.5, -1.0, 0.0, twobythree},         // [0,1]
         {-onebysix, 1.0, -2.0, fourbythree}},

        // shapefunctiontype = 3 (Interior node)
        {{onebysix, 1.0, 2.0, fourbythree},
         {-0.5, -1.0, 0.0, twobythree},
         {0.5, -1.0, 0.0, twobythree},
         {-onebysix, 1.0, -2.0, fourbythree}},

        // shapefunctiontype = 4 (Near boundary node, high boundary)
        {{onebysix, 1.0, 2.0, fourbythree},
         {-0.5, -1.0, 0.0, twobythree},
         {onebythree, -1.0, 0.0, twobythree},
         {0.0, 0.0, 0.0, 0.0}}};

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_GPU_CONSTANT
#endif
    constexpr CubicSplineDerCoeff cubic_splineder_coeffs[5][4] = {
        {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},

        // shapefunctiontype = 1 (Boundary node)
        {{0.5, 2.0, 2.0},    // [-2,-1] : 0.5*zi^2 + 2*zi + 2
         {-0.5, 0.0, 1.0},   // [-1,0]  : -0.5*zi^2 + 1
         {0.5, 0.0, -1.0},   // [0,1]   : 0.5*zi^2 - 1
         {-0.5, 2.0, -2.0}}, // [1,2]   : -0.5*zi^2 + 2*zi - 2

        // shapefunctiontype = 2 (Near boundary node)
        {
            {0, 0, 0},
            {-1.0, -2.0, 0.0}, // [-1,0]  : -zi^2 - 2*zi
            {1.5, -2.0, 0.0},  // [0,1]   : 1.5*zi^2 - 2*zi
            {-0.5, 2.0, -2.0}  // [1,2]   : -0.5*zi^2 + 2*zi - 2
        },                     // filler

        // shapefunctiontype = 3 (Interior node)
        {{0.5, 2.0, 2.0},    // [-2,-1] : 0.5*zi^2 + 2*zi + 2
         {-1.5, -2.0, 0.0},  // [-1,0]  : -1.5*zi^2 - 2*zi
         {1.5, -2.0, 0.0},   // [0,1]   : 1.5*zi^2 - 2*zi
         {-0.5, 2.0, -2.0}}, // [1,2]   : -0.5*zi^2 + 2*zi - 2

        // shapefunctiontype = 4 (Alt interior node)
        {{0.5, 2.0, 2.0},   // [-2,-1] : 0.5*zi^2 + 2*zi + 2
         {-1.5, -2.0, 0.0}, // [-1,0]  : -1.5*zi^2 - 2*zi
         {1.0, -2.0, 0.0},  // [0,1]   : 1.0*zi^2 - 2*zi
         {0, 0, 0}}         // filler
};

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_GPU_CONSTANT
#endif
    constexpr QuadraticSplineCoeff quad_spline_coeffs[5][4] = {
        {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},

        // shapefunctiontype = 1 (Boundary node, low/high end)
        {{twobythree, 2.0, 1.5},   // [-1.5,-0.5] : (2/3)zi² + 2zi + 1.5
         {-fourbythree, 0.0, 1.0}, // [-0.5,0]    : 1 - (4/3)zi²
         {-fourbythree, 0.0, 1.0}, // [0,0.5]     : 1 - (4/3)zi²
         {twobythree, -2.0, 1.5}}, // [0.5,1.5]   : (2/3)zi² - 2zi + 1.5

        // shapefunctiontype = 2
        {{4.0 / 3.0, 8.0 / 3.0, 4.0 / 3.0},
         {-7.0 / 6.0, 1.0 / 6.0, 17.0 / 24.0},
         {0.5, -1.5, 1.125}, // was 0.5625; correct is 9/8 = 1.125
         {0, 0, 0}},

        // shapefunctiontype = 3 (Interior node)
        {{0.5, 1.5, 1.125},  // [-1.5,-0.5] : (2zi+3)²/8 = 0.5zi² + 1.5zi + 9/8
         {-1.0, 0.0, 0.75},  // [-0.5,0.5]  : (-8zi²+6)/8 = -zi² + 0.75
         {0.5, -1.5, 1.125}, // [0.5,1.5]   : (2zi-3)²/8 = 0.5zi² -1.5zi + 9/8
         {0, 0, 0}},         // filler

        // shapefunctiontype = 4 (Near boundary node, high end)
        {{0.5, 1.5, 1.125},                     // [-1.5,-0.5]
         {-7.0 / 6.0, -1.0 / 6.0, 17.0 / 24.0}, // [-0.5,0.5]
         {4.0 / 3.0, -8.0 / 3.0, 4.0 / 3.0},    // [0.5,1.0]  ← corrected
         {0, 0, 0}}};

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_GPU_CONSTANT
#endif
    constexpr QuadraticSplineDerCoeff quad_splineder_coeffs[5][4] = {
        {{0, 0}, {0, 0}, {0, 0}, {0, 0}},

        // shapefunctiontype = 1 (Boundary node)
        {{4.0 / 3.0, 2.0},   // [-1.5,-0.5] : 4/3*zi + 2
         {-8.0 / 3.0, 0.0},  // [-0.5,0]    : -8/3*zi
         {-8.0 / 3.0, 0.0},  // [0,0.5]     : -8/3*zi
         {4.0 / 3.0, -2.0}}, // [0.5,1.5]   : 4/3*zi - 2

        // shapefunctiontype = 2 (Near boundary node)
        {{8.0 / 3.0, 8.0 / 3.0},   // [-1.0,-0.5] : 8/3*(zi+1) = 8/3*zi + 8/3
         {-14.0 / 6.0, 1.0 / 6.0}, // [-0.5,0.5]  : -14/6*zi + 1/6
         {1.0, -1.5},              // [0.5,1.5]   : zi - 1.5
         {0, 0}},                  // filler

        // shapefunctiontype = 3 (Interior node)
        {{2.0 / 2.0, 3.0 / 2.0},  // [-1.5,-0.5] : (2zi+3)/2 = zi + 1.5
         {-2.0, 0.0},             // [-0.5,0.5]  : -2zi
         {2.0 / 2.0, -3.0 / 2.0}, // [0.5,1.5]   : (2zi-3)/2 = zi - 1.5
         {0, 0}},                 // filler

        // shapefunctiontype = 4 (Alt interior node)
        {{1.0, 1.5},                // [-1.5,-0.5] : zi + 1.5
         {-14.0 / 6.0, -1.0 / 6.0}, // [-0.5,0.5]  : -14/6*zi - 1/6
         {8.0 / 3.0, -8.0 / 3.0},   // [0.5,1.0]   : 8/3*(zi-1) = 8/3*zi - 8/3
         {0, 0}}                    // filler
};

/**
 * @brief Computes the interval index for a quadratic B‑spline basis function.
 *
 * Determines which polynomial segment should be used for evaluating a
 * quadratic B‑spline basis function, based on the shape-function type and
 * the shifted/scaled coordinate zi. Returns -1 if zi lies outside the
 * supported range.
 *
 * @param[in] shapefunctiontype  Type of spline (1=boundary, 2=near-boundary,
 *                               3=interior, 4=alt interior)
 * @param[in] zi                 Local coordinate in spline space
 *
 * @return Interval index in [0,3], or -1 if zi is outside support.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
interval_idx_quadratic(int shapefunctiontype, amrex::Real zi)
{
    amrex::Real shift = (shapefunctiontype == 2) ? 1.0 : 1.5;
    int raw = static_cast<int>(std::floor((zi + shift) * 2));
    if (raw < 0 || raw >= 6)
        return -1;
#ifdef __CUDA_ARCH__
    return interval_map_quadbspline[shapefunctiontype][raw];
#else
    return interval_map_quadbspline_host[shapefunctiontype][raw];
#endif
}

/**
 * @brief Evaluates a 1D cubic B‑spline basis function.
 *
 * Selects the correct cubic polynomial segment based on floor(zi) and
 * evaluates it using precomputed coefficients for the given shape-function
 * type. Supports boundary, near-boundary, and interior spline variants.
 *
 * @param[in] shapefunctiontype  Spline type controlling which coefficient set
 * to use
 * @param[in] zi                 Local coordinate in spline space
 *
 * @return Value of the cubic B‑spline basis function at zi.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cubicspline_1d(int shapefunctiontype, amrex::Real zi)
{
    amrex::Real value = 0.0;
    amrex::Real zi_cubed = zi * zi * zi;
    amrex::Real zi_squared = zi * zi;

    int interval = static_cast<int>(std::floor(zi)) + 2; // maps [-2,-1]→0, etc.
    const auto &c = cubic_spline_coeffs[shapefunctiontype][interval];
    value = c.a * zi_cubed + c.b * zi_squared + c.c * zi + c.d;
    return value;
}

/**
 * @brief Evaluates the derivative of a 1D cubic B‑spline basis function.
 *
 * Uses precomputed derivative coefficients for each interval and shape-function
 * type. Computes dN/dzi for the appropriate cubic polynomial segment.
 *
 * @param[in] shapefunctiontype  Spline type controlling which derivative
 * coefficients to use
 * @param[in] zi                 Local coordinate in spline space
 *
 * @return Derivative of the cubic B‑spline basis function at zi.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cubicspline_1d_der(int shapefunctiontype, amrex::Real zi)
{
    amrex::Real value = 0.0;
    amrex::Real zi_squared = zi * zi;
    int interval = static_cast<int>(std::floor(zi)) + 2; // maps [-2,-1]→0, etc.
    const auto &c = cubic_splineder_coeffs[shapefunctiontype][interval];
    value = c.a * zi_squared + c.b * zi + c.c;
    return value;
}

/**
 * @brief Evaluates a 1D quadratic B‑spline basis function.
 *
 * Supports two evaluation modes:
 *   - Coefficient-based evaluation
 *
 * Handles boundary, near-boundary, and interior spline types.
 *
 * @param[in] shapefunctiontype  Spline type controlling which polynomial to use
 * @param[in] zi                 Local coordinate in spline space
 *
 * @return Value of the quadratic B‑spline basis function at zi.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
quadraticspline_1d(int shapefunctiontype, amrex::Real zi)
{
    amrex::Real value = 0.0;
    amrex::Real zi_squared = zi * zi;
    int interval = interval_idx_quadratic(shapefunctiontype, zi);
    if (interval < 0)
        return 0.0;
    const auto &c = quad_spline_coeffs[shapefunctiontype][interval];
    value = c.a * zi_squared + c.b * zi + c.c;
    return value;
}

/**
 * @brief Evaluates the derivative of a 1D quadratic B‑spline basis function.
 *
 * Uses precomputed derivative coefficients for each interval and spline type.
 * Computes dN/dzi for the appropriate quadratic polynomial segment.
 *
 * @param[in] shapefunctiontype  Spline type controlling which derivative
 * coefficients to use
 * @param[in] zi                 Local coordinate in spline space
 *
 * @return Derivative of the quadratic B‑spline basis function at zi.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
quadraticspline_1d_der(int shapefunctiontype, amrex::Real zi)
{
    amrex::Real value = 0.0;
    int interval =
        interval_idx_quadratic(shapefunctiontype, zi); // arithmetic mapping
    if (interval < 0)
        return 0.0;

    const auto &c = quad_splineder_coeffs[shapefunctiontype][interval];
    value = c.a * zi + c.b;
    return value;
}

/**
 * @brief Evaluates a normalized spherical Gaussian kernel.
 *
 * Computes exp(-‖xi - xp‖² / r0²) multiplied by the correct normalization
 * constant for 1D, 2D, or 3D. Used for smooth particle kernels, density
 * estimation, or level-set approximations.
 *
 * @param[in] xi   Evaluation point in physical space
 * @param[in] xp   Kernel center (particle position)
 * @param[in] r0   Kernel radius (standard deviation-like parameter)
 *
 * @return Gaussian kernel value at xi.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
spherical_gaussian(amrex::Real xi[AMREX_SPACEDIM],
                   amrex::Real xp[AMREX_SPACEDIM],
                   amrex::Real r0)
{
    amrex::Real r2 = 0.0;
#if (AMREX_SPACEDIM == 1)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]);
    // Normalization factor in 1D: (π^(-0.5)) * r0^(-1)
    return std::exp(-r2 / (r0 * r0)) * std::pow(PI, -0.5) * std::pow(r0, -1.0);

#elif (AMREX_SPACEDIM == 2)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]);
    // Normalization factor in 2D: (π^(-1)) * r0^(-2)
    return std::exp(-r2 / (r0 * r0)) * std::pow(PI, -1.0) * std::pow(r0, -2.0);

#elif (AMREX_SPACEDIM == 3)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]) +
         (xi[2] - xp[2]) * (xi[2] - xp[2]);
    // Normalization factor in 3D: (π^(-1.5)) * r0^(-3)
    return std::exp(-r2 / (r0 * r0)) * std::pow(PI, -1.5) * std::pow(r0, -3.0);
#endif
}

/**
 * @brief Evaluates a compact-support box kernel.
 *
 * Returns 1.0 if the squared distance between xi and xp is less than r0²,
 * otherwise returns 0.0. Useful for simple particle-to-grid transfer or
 * indicator functions.
 *
 * @param[in] xi   Evaluation point
 * @param[in] xp   Kernel center
 * @param[in] r0   Support radius
 *
 * @return 1.0 if inside support, 0.0 otherwise.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
box_kernel(amrex::Real xi[AMREX_SPACEDIM],
           amrex::Real xp[AMREX_SPACEDIM],
           amrex::Real r0)
{
    amrex::Real r2 = 0.0;
#if (AMREX_SPACEDIM == 1)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]);

#elif (AMREX_SPACEDIM == 2)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]);

#elif (AMREX_SPACEDIM == 3)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]) +
         (xi[2] - xp[2]) * (xi[2] - xp[2]);
#endif

    return (r2 < r0 * r0) ? 1.0 : 0.0;
}

/**
 * @brief Computes a signed-distance-like level-set value.
 *
 * Returns:
 *   - minval if inside radius r0
 *   - sqrt(r²) outside
 *   - clamped to maxval beyond 3*r0
 *
 * Used for implicit geometry, particle-based level sets, or interface tracking.
 *
 * @param[in] xi      Evaluation point
 * @param[in] xp      Center point
 * @param[in] r0      Radius threshold
 * @param[in] minval  Value returned inside r0
 * @param[in] maxval  Maximum clamped value
 *
 * @return Level-set-like scalar value.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
levelset(amrex::Real xi[AMREX_SPACEDIM],
         amrex::Real xp[AMREX_SPACEDIM],
         amrex::Real r0,
         amrex::Real minval = -1.0,
         amrex::Real maxval = 1.0)
{
    amrex::Real r2 = 0.0;

#if (AMREX_SPACEDIM == 1)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]);

#elif (AMREX_SPACEDIM == 2)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]);

#elif (AMREX_SPACEDIM == 3)
    r2 = (xi[0] - xp[0]) * (xi[0] - xp[0]) + (xi[1] - xp[1]) * (xi[1] - xp[1]) +
         (xi[2] - xp[2]) * (xi[2] - xp[2]);
#endif

    amrex::Real f1 = (r2 < r0 * r0) ? minval : std::sqrt(r2);
    amrex::Real f2 = (f1 > 3.0 * r0) ? maxval : f1;

    return f2;
}

/**
 * @brief Evaluates a 1D hat (tent) basis function.
 *
 * Produces a triangular function centered at xp with width hatsize and
 * peak value proportional to 1/hatsize. Used for simple interpolation or
 * particle-to-grid transfer.
 *
 * @param[in] xi       Evaluation coordinate
 * @param[in] xp       Center coordinate
 * @param[in] hatsize  Width of the hat function
 *
 * @return Hat function value at xi.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
hat1d(amrex::Real xi, amrex::Real xp, amrex::Real hatsize)
{
    amrex::Real zeta = xi - xp;
    amrex::Real halfhatsize = hatsize * HALF_CONST;
    amrex::Real funcval = shunya;
    amrex::Real maxval = eka / halfhatsize;
    amrex::Real slope = four / (hatsize * hatsize);

    if (amrex::Math::abs(zeta) < halfhatsize)
    {
        if (zeta > 0.0)
        {
            funcval = -slope * zeta + maxval;
        }
        else
        {
            funcval = slope * zeta + maxval;
        }
    }
    return (funcval);
}

/**
 * @brief Evaluates an N-dimensional hat function as a tensor product of 1D
 * hats.
 *
 * Supports 1D, 2D, and 3D. Each dimension uses hat1d() with its own width.
 *
 * @param[in] xi       Evaluation point
 * @param[in] xp       Center point
 * @param[in] hatsize  Hat width per dimension
 *
 * @return Tensor-product hat basis value.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
hatND(amrex::Real xi[AMREX_SPACEDIM],
      amrex::Real xp[AMREX_SPACEDIM],
      amrex::Real hatsize[AMREX_SPACEDIM])
{
    amrex::Real funcval = eka;

#if (AMREX_SPACEDIM == 1)
    funcval *= hat1d(xi[0], xp[0], hatsize[0]);

#elif (AMREX_SPACEDIM == 2)
    funcval *= hat1d(xi[0], xp[0], hatsize[0]);
    funcval *= hat1d(xi[1], xp[1], hatsize[1]);

#elif (AMREX_SPACEDIM == 3)
    funcval *= hat1d(xi[0], xp[0], hatsize[0]);
    funcval *= hat1d(xi[1], xp[1], hatsize[1]);
    funcval *= hat1d(xi[2], xp[2], hatsize[2]);
#endif

    return funcval;
}

/**
 * @brief Performs linear, bilinear, or trilinear interpolation on nodal data.
 *
 * Computes the interpolated value at particle position xp using nodal values
 * stored in nodaldata. Dimensionality is determined by AMREX_SPACEDIM.
 *
 * @param[in] xp        Particle position
 * @param[in] iv        Cell index containing xp
 * @param[in] plo       Physical lower bounds of the domain
 * @param[in] dx        Grid spacing
 * @param[in] nodaldata Nodal field values
 * @param[in] comp      Component index to interpolate
 *
 * @return Interpolated nodal value at xp.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
bilin_interp(amrex::Real xp[AMREX_SPACEDIM],
             amrex::IntVect iv,
             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
             amrex::Array4<amrex::Real> nodaldata,
             int comp)
{
    amrex::Real value = shunya;

#if (AMREX_SPACEDIM == 1)
    int i = iv[XDIR];
    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    value = (1.0 - alpha) * nodaldata(i, 0, 0, comp) +
            alpha * nodaldata(i + 1, 0, 0, comp);

#elif (AMREX_SPACEDIM == 2)
    int i = iv[XDIR];
    int j = iv[YDIR];
    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    amrex::Real beta = (xp[1] - (plo[1] + j * dx[1])) / dx[1];

    value = (1.0 - alpha) * (1.0 - beta) * nodaldata(i, j, 0, comp);
    value += alpha * (1.0 - beta) * nodaldata(i + 1, j, 0, comp);
    value += (1.0 - alpha) * beta * nodaldata(i, j + 1, 0, comp);
    value += alpha * beta * nodaldata(i + 1, j + 1, 0, comp);

#elif (AMREX_SPACEDIM == 3)
    int i = iv[XDIR];
    int j = iv[YDIR];
    int k = iv[ZDIR];

    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    amrex::Real beta = (xp[1] - (plo[1] + j * dx[1])) / dx[1];
    amrex::Real gama = (xp[2] - (plo[2] + k * dx[2])) / dx[2];

    value =
        (1.0 - alpha) * (1.0 - beta) * (1.0 - gama) * nodaldata(i, j, k, comp);
    value += alpha * (1.0 - beta) * (1.0 - gama) * nodaldata(i + 1, j, k, comp);
    value += (1.0 - alpha) * beta * (1.0 - gama) * nodaldata(i, j + 1, k, comp);
    value += alpha * beta * (1.0 - gama) * nodaldata(i + 1, j + 1, k, comp);

    value += (1.0 - alpha) * (1.0 - beta) * gama * nodaldata(i, j, k + 1, comp);
    value += alpha * (1.0 - beta) * gama * nodaldata(i + 1, j, k + 1, comp);
    value += (1.0 - alpha) * beta * gama * nodaldata(i, j + 1, k + 1, comp);
    value += alpha * beta * gama * nodaldata(i + 1, j + 1, k + 1, comp);
#endif

    return value;
}

/**
 * @brief Interpolates two nodal field components at a particle position in a
 * single traversal, reusing the basis weights for both.
 *
 * Structurally identical to bilin_interp but accumulates into vals[0] (comp0)
 * and vals[1] (comp1) using the same weight for each node, eliminating the
 * redundant weight computation when two fields are sampled at the same
 * location.
 *
 * @param[in]  xp        Particle position
 * @param[in]  iv        Base cell index
 * @param[in]  plo       Domain lower bounds
 * @param[in]  dx        Grid spacing
 * @param[in]  nodaldata Nodal field array
 * @param[in]  comp0     First component index
 * @param[in]  comp1     Second component index
 * @param[out] vals      vals[0] = interpolated comp0, vals[1] = interpolated
 * comp1
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bilin_interp_two(amrex::Real xp[AMREX_SPACEDIM],
                 amrex::IntVect iv,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                 amrex::Array4<amrex::Real> nodaldata,
                 int comp0,
                 int comp1,
                 amrex::Real vals[2])
{
    vals[0] = shunya;
    vals[1] = shunya;

#if (AMREX_SPACEDIM == 1)
    int i = iv[XDIR];
    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    amrex::Real w0 = 1.0 - alpha;
    amrex::Real w1 = alpha;
    vals[0] =
        w0 * nodaldata(i, 0, 0, comp0) + w1 * nodaldata(i + 1, 0, 0, comp0);
    vals[1] =
        w0 * nodaldata(i, 0, 0, comp1) + w1 * nodaldata(i + 1, 0, 0, comp1);

#elif (AMREX_SPACEDIM == 2)
    int i = iv[XDIR];
    int j = iv[YDIR];
    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    amrex::Real beta = (xp[1] - (plo[1] + j * dx[1])) / dx[1];
    amrex::Real w00 = (1.0 - alpha) * (1.0 - beta);
    amrex::Real w10 = alpha * (1.0 - beta);
    amrex::Real w01 = (1.0 - alpha) * beta;
    amrex::Real w11 = alpha * beta;
    vals[0] = w00 * nodaldata(i, j, 0, comp0) +
              w10 * nodaldata(i + 1, j, 0, comp0) +
              w01 * nodaldata(i, j + 1, 0, comp0) +
              w11 * nodaldata(i + 1, j + 1, 0, comp0);
    vals[1] = w00 * nodaldata(i, j, 0, comp1) +
              w10 * nodaldata(i + 1, j, 0, comp1) +
              w01 * nodaldata(i, j + 1, 0, comp1) +
              w11 * nodaldata(i + 1, j + 1, 0, comp1);

#elif (AMREX_SPACEDIM == 3)
    int i = iv[XDIR];
    int j = iv[YDIR];
    int k = iv[ZDIR];
    amrex::Real alpha = (xp[0] - (plo[0] + i * dx[0])) / dx[0];
    amrex::Real beta = (xp[1] - (plo[1] + j * dx[1])) / dx[1];
    amrex::Real gama = (xp[2] - (plo[2] + k * dx[2])) / dx[2];
    amrex::Real w000 = (1.0 - alpha) * (1.0 - beta) * (1.0 - gama);
    amrex::Real w100 = alpha * (1.0 - beta) * (1.0 - gama);
    amrex::Real w010 = (1.0 - alpha) * beta * (1.0 - gama);
    amrex::Real w110 = alpha * beta * (1.0 - gama);
    amrex::Real w001 = (1.0 - alpha) * (1.0 - beta) * gama;
    amrex::Real w101 = alpha * (1.0 - beta) * gama;
    amrex::Real w011 = (1.0 - alpha) * beta * gama;
    amrex::Real w111 = alpha * beta * gama;
    vals[0] = w000 * nodaldata(i, j, k, comp0) +
              w100 * nodaldata(i + 1, j, k, comp0) +
              w010 * nodaldata(i, j + 1, k, comp0) +
              w110 * nodaldata(i + 1, j + 1, k, comp0) +
              w001 * nodaldata(i, j, k + 1, comp0) +
              w101 * nodaldata(i + 1, j, k + 1, comp0) +
              w011 * nodaldata(i, j + 1, k + 1, comp0) +
              w111 * nodaldata(i + 1, j + 1, k + 1, comp0);
    vals[1] = w000 * nodaldata(i, j, k, comp1) +
              w100 * nodaldata(i + 1, j, k, comp1) +
              w010 * nodaldata(i, j + 1, k, comp1) +
              w110 * nodaldata(i + 1, j + 1, k, comp1) +
              w001 * nodaldata(i, j, k + 1, comp1) +
              w101 * nodaldata(i + 1, j, k + 1, comp1) +
              w011 * nodaldata(i, j + 1, k + 1, comp1) +
              w111 * nodaldata(i + 1, j + 1, k + 1, comp1);
#endif
}

/**
 * @brief Performs quadratic B‑spline interpolation in 1D, 2D, or 3D.
 *
 * Evaluates tensor-product quadratic B‑splines over the stencil defined by
 * min_index and max_index. Boundary-aware shape-function types are selected
 * automatically based on lo/hi and periodicity.
 *
 * @param[in] xp         Particle position
 * @param[in] iv         Base cell index
 * @param[in] min_index  Minimum stencil offset
 * @param[in] max_index  Maximum stencil offset
 * @param[in] plo        Domain lower bounds
 * @param[in] dx         Grid spacing
 * @param[in] nodaldata  Nodal field values
 * @param[in] comp       Component index
 * @param[in] lo         Lower physical boundary indices
 * @param[in] hi         Upper physical boundary indices
 *
 * @return Quadratic B‑spline interpolated value.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
quadratic_interp(amrex::Real xp[AMREX_SPACEDIM],
                 IntVect iv,
                 IntVect min_index,
                 IntVect max_index,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                 amrex::Array4<amrex::Real> nodaldata,
                 int comp,
                 const amrex::GpuArray<int, AMREX_SPACEDIM> &periodic,
                 const int *lo,
                 const int *hi)
{
    amrex::Real value = shunya;

#if (AMREX_SPACEDIM == 1)
    amrex::Real rx;
    amrex::Real lval;
    int shapetypex;

    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];

    for (int l = lmin; l < lmax; l++)
    {
        if (i + l == lo[0])
            shapetypex = 1;
        else if (i + l == lo[0] + 1)
            shapetypex = 2;
        else if (i + l == hi[0] + 1)
            shapetypex = 1;
        else if (i + l == hi[0])
            shapetypex = 4;
        else
            shapetypex = 3;
        rx = (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
        lval = quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
        value += lval * nodaldata(i + l, 0, 0, comp);
    }

#elif (AMREX_SPACEDIM == 2)
    amrex::Real rx, ry;
    amrex::Real lval, mval;
    int shapetypex, shapetypey;

    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];

    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    for (int m = mmin; m < mmax; m++)
    {
        for (int l = lmin; l < lmax; l++)
        {
            if (i + l == lo[0])
                shapetypex = 1;
            else if (i + l == lo[0] + 1)
                shapetypex = 2;
            else if (i + l == hi[0] + 1)
                shapetypex = 1;
            else if (i + l == hi[0])
                shapetypex = 4;
            else
                shapetypex = 3;

            if (j + m == lo[1])
                shapetypey = 1;
            else if (j + m == lo[1] + 1)
                shapetypey = 2;
            else if (j + m == hi[1] + 1)
                shapetypey = 1;
            else if (j + m == hi[1])
                shapetypey = 4;
            else
                shapetypey = 3;

            rx = (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
            ry = (xp[YDIR] - (plo[YDIR] + (j + m) * dx[YDIR])) / dx[YDIR];

            lval = quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
            mval = quadraticspline_1d(periodic[1] ? 3 : shapetypey, ry);

            value += lval * mval * nodaldata(i + l, j + m, 0, comp);
        }
    }

#elif (AMREX_SPACEDIM == 3)
    amrex::Real rx, ry, rz;
    amrex::Real lval, mval, nval;
    int shapetypex, shapetypey, shapetypez;

    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];

    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];

    int k = iv[ZDIR];
    int nmin = min_index[ZDIR];
    int nmax = max_index[ZDIR];

    for (int n = nmin; n < nmax; n++)
    {
        for (int m = mmin; m < mmax; m++)
        {
            for (int l = lmin; l < lmax; l++)
            {
                if (i + l == lo[0])
                    shapetypex = 1;
                else if (i + l == lo[0] + 1)
                    shapetypex = 2;
                else if (i + l == hi[0] + 1)
                    shapetypex = 1;
                else if (i + l == hi[0])
                    shapetypex = 4;
                else
                    shapetypex = 3;

                if (j + m == lo[1])
                    shapetypey = 1;
                else if (j + m == lo[1] + 1)
                    shapetypey = 2;
                else if (j + m == hi[1] + 1)
                    shapetypey = 1;
                else if (j + m == hi[1])
                    shapetypey = 4;
                else
                    shapetypey = 3;

                if (k + n == lo[2])
                    shapetypez = 1;
                else if (k + n == lo[2] + 1)
                    shapetypez = 2;
                else if (k + n == hi[2] + 1)
                    shapetypez = 1;
                else if (k + n == hi[2])
                    shapetypez = 4;
                else
                    shapetypez = 3;

                rx = (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
                ry = (xp[YDIR] - (plo[YDIR] + (j + m) * dx[YDIR])) / dx[YDIR];
                rz = (xp[ZDIR] - (plo[ZDIR] + (k + n) * dx[ZDIR])) / dx[ZDIR];

                lval = quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
                mval = quadraticspline_1d(periodic[1] ? 3 : shapetypey, ry);
                nval = quadraticspline_1d(periodic[2] ? 3 : shapetypez, rz);
                value +=
                    lval * mval * nval * nodaldata(i + l, j + m, k + n, comp);
            }
        }
    }
#endif

    return (value);
}

/**
 * @brief Quadratic B‑spline interpolation for two components in one stencil
 * pass.
 *
 * Computes the tensor-product quadratic B‑spline weight for each stencil node
 * once and accumulates it into vals[0] (comp0) and vals[1] (comp1).
 *
 * @param[in]  xp        Particle position
 * @param[in]  iv        Base cell index
 * @param[in]  min_index Minimum stencil offset per dimension
 * @param[in]  max_index Maximum stencil offset per dimension
 * @param[in]  plo       Domain lower bounds
 * @param[in]  dx        Grid spacing
 * @param[in]  nodaldata Nodal field array
 * @param[in]  comp0     First component index
 * @param[in]  comp1     Second component index
 * @param[in]  lo        Lower physical boundary indices
 * @param[in]  hi        Upper physical boundary indices
 * @param[out] vals      vals[0] = interpolated comp0, vals[1] = interpolated
 * comp1
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
quadratic_interp_two(amrex::Real xp[AMREX_SPACEDIM],
                     IntVect iv,
                     IntVect min_index,
                     IntVect max_index,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                     amrex::Array4<amrex::Real> nodaldata,
                     int comp0,
                     int comp1,
                     const amrex::GpuArray<int, AMREX_SPACEDIM> &periodic,
                     const int *lo,
                     const int *hi,
                     amrex::Real vals[2])
{
    vals[0] = shunya;
    vals[1] = shunya;

#if (AMREX_SPACEDIM == 1)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    for (int l = lmin; l < lmax; l++)
    {
        int shapetypex;
        if (i + l == lo[0] || i + l == hi[0] + 1)
            shapetypex = 1;
        else if (i + l == lo[0] + 1)
            shapetypex = 2;
        else if (i + l == hi[0])
            shapetypex = 4;
        else
            shapetypex = 3;
        amrex::Real rx =
            (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
        amrex::Real lval = quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
        vals[0] += lval * nodaldata(i + l, 0, 0, comp0);
        vals[1] += lval * nodaldata(i + l, 0, 0, comp1);
    }

#elif (AMREX_SPACEDIM == 2)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    for (int m = mmin; m < mmax; m++)
    {
        int shapetypey;
        if (j + m == lo[1] || j + m == hi[1] + 1)
            shapetypey = 1;
        else if (j + m == lo[1] + 1)
            shapetypey = 2;
        else if (j + m == hi[1])
            shapetypey = 4;
        else
            shapetypey = 3;
        amrex::Real ry =
            (xp[YDIR] - (plo[YDIR] + (j + m) * dx[YDIR])) / dx[YDIR];
        amrex::Real mval = quadraticspline_1d(periodic[1] ? 3 : shapetypey, ry);
        for (int l = lmin; l < lmax; l++)
        {
            int shapetypex;
            if (i + l == lo[0] || i + l == hi[0] + 1)
                shapetypex = 1;
            else if (i + l == lo[0] + 1)
                shapetypex = 2;
            else if (i + l == hi[0])
                shapetypex = 4;
            else
                shapetypex = 3;
            amrex::Real rx =
                (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
            amrex::Real lval =
                quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
            amrex::Real weight = lval * mval;
            vals[0] += weight * nodaldata(i + l, j + m, 0, comp0);
            vals[1] += weight * nodaldata(i + l, j + m, 0, comp1);
        }
    }

#elif (AMREX_SPACEDIM == 3)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    int k = iv[ZDIR];
    int nmin = min_index[ZDIR];
    int nmax = max_index[ZDIR];
    for (int n = nmin; n < nmax; n++)
    {
        int shapetypez;
        if (k + n == lo[2] || k + n == hi[2] + 1)
            shapetypez = 1;
        else if (k + n == lo[2] + 1)
            shapetypez = 2;
        else if (k + n == hi[2])
            shapetypez = 4;
        else
            shapetypez = 3;
        amrex::Real rz =
            (xp[ZDIR] - (plo[ZDIR] + (k + n) * dx[ZDIR])) / dx[ZDIR];
        amrex::Real nval = quadraticspline_1d(periodic[2] ? 3 : shapetypez, rz);
        for (int m = mmin; m < mmax; m++)
        {
            int shapetypey;
            if (j + m == lo[1] || j + m == hi[1] + 1)
                shapetypey = 1;
            else if (j + m == lo[1] + 1)
                shapetypey = 2;
            else if (j + m == hi[1])
                shapetypey = 4;
            else
                shapetypey = 3;
            amrex::Real ry =
                (xp[YDIR] - (plo[YDIR] + (j + m) * dx[YDIR])) / dx[YDIR];
            amrex::Real mval =
                quadraticspline_1d(periodic[1] ? 3 : shapetypey, ry);
            for (int l = lmin; l < lmax; l++)
            {
                int shapetypex;
                if (i + l == lo[0] || i + l == hi[0] + 1)
                    shapetypex = 1;
                else if (i + l == lo[0] + 1)
                    shapetypex = 2;
                else if (i + l == hi[0])
                    shapetypex = 4;
                else
                    shapetypex = 3;
                amrex::Real rx =
                    (xp[XDIR] - (plo[XDIR] + (i + l) * dx[XDIR])) / dx[XDIR];
                amrex::Real lval =
                    quadraticspline_1d(periodic[0] ? 3 : shapetypex, rx);
                amrex::Real weight = lval * mval * nval;
                vals[0] += weight * nodaldata(i + l, j + m, k + n, comp0);
                vals[1] += weight * nodaldata(i + l, j + m, k + n, comp1);
            }
        }
    }
#endif
}

/**
 * @brief Performs cubic B‑spline interpolation in 1D, 2D, or 3D.
 *
 * Evaluates tensor-product cubic B‑splines using precomputed coefficients and
 * boundary-aware shape-function types. Supports arbitrary stencils.
 *
 * @param[in] xp         Particle position
 * @param[in] iv         Base cell index
 * @param[in] min_index  Minimum stencil offset
 * @param[in] max_index  Maximum stencil offset
 * @param[in] plo        Domain lower bounds
 * @param[in] dx         Grid spacing
 * @param[in] nodaldata  Nodal field values
 * @param[in] comp       Component index
 * @param[in] lo         Lower physical boundary indices
 * @param[in] hi         Upper physical boundary indices
 *
 * @return Cubic B‑spline interpolated value.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
cubic_interp(amrex::Real xp[AMREX_SPACEDIM],
             IntVect iv,
             IntVect min_index,
             IntVect max_index,
             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
             amrex::Array4<amrex::Real> nodaldata,
             int comp,
             const amrex::GpuArray<int, AMREX_SPACEDIM> &periodic,
             const int *lo,
             const int *hi)
{
    amrex::Real value = shunya;

#if (AMREX_SPACEDIM == 1)

    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];

    amrex::Real inv_dx = 1.0 / dx[0];
    for (int l = lmin; l < lmax; ++l)
    {
        int shapetypex;
        if (i + l == lo[0])
            shapetypex = 1;
        else if (i + l == lo[0] + 1)
            shapetypex = 2;
        else if (i + l == hi[0] + 1)
            shapetypex = 1;
        else if (i + l == hi[0])
            shapetypex = 4;
        else
            shapetypex = 3;

        amrex::Real cell_center = plo[0] + (i + l) * dx[0];

        amrex::Real rx = (xp[0] - cell_center) * inv_dx;
        amrex::Real lval = cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
        value += lval * nodaldata(i + l, 0, 0, comp);
    }

#elif (AMREX_SPACEDIM == 2)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];

    amrex::Real inv_dx = 1.0 / dx[0];
    amrex::Real inv_dy = 1.0 / dx[1];

    for (int m = mmin; m < mmax; ++m)
    {
        for (int l = lmin; l < lmax; ++l)
        {
            int shapetypex, shapetypey;
            // X direction
            if (i + l == lo[0])
                shapetypex = 1;
            else if (i + l == lo[0] + 1)
                shapetypex = 2;
            else if (i + l == hi[0] + 1)
                shapetypex = 1;
            else if (i + l == hi[0])
                shapetypex = 4;
            else
                shapetypex = 3;
            // Y direction
            if (j + m == lo[1])
                shapetypey = 1;
            else if (j + m == lo[1] + 1)
                shapetypey = 2;
            else if (j + m == hi[1] + 1)
                shapetypey = 1;
            else if (j + m == hi[1])
                shapetypey = 4;
            else
                shapetypey = 3;

            amrex::Real cell_centerx = plo[0] + (i + l) * dx[0];
            amrex::Real cell_centery = plo[1] + (j + m) * dx[1];

            amrex::Real rx = (xp[0] - (cell_centerx)) * inv_dx;
            amrex::Real ry = (xp[1] - (cell_centery)) * inv_dy;

            amrex::Real lval = cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
            amrex::Real mval = cubicspline_1d(periodic[1] ? 3 : shapetypey, ry);

            value += lval * mval * nodaldata(i + l, j + m, 0, comp);
        }
    }

#elif (AMREX_SPACEDIM == 3)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    int k = iv[ZDIR];
    int nmin = min_index[ZDIR];
    int nmax = max_index[ZDIR];

    amrex::Real inv_dx = 1.0 / dx[0];
    amrex::Real inv_dy = 1.0 / dx[1];
    amrex::Real inv_dz = 1.0 / dx[2];

    for (int n = nmin; n < nmax; ++n)
    {
        for (int m = mmin; m < mmax; ++m)
        {
            for (int l = lmin; l < lmax; ++l)
            {
                int shapetypex, shapetypey, shapetypez;
                // X direction
                if (i + l == lo[0])
                    shapetypex = 1;
                else if (i + l == lo[0] + 1)
                    shapetypex = 2;
                else if (i + l == hi[0] + 1)
                    shapetypex = 1;
                else if (i + l == hi[0])
                    shapetypex = 4;
                else
                    shapetypex = 3;
                // Y direction
                if (j + m == lo[1])
                    shapetypey = 1;
                else if (j + m == lo[1] + 1)
                    shapetypey = 2;
                else if (j + m == hi[1] + 1)
                    shapetypey = 1;
                else if (j + m == hi[1])
                    shapetypey = 4;
                else
                    shapetypey = 3;
                // Z direction
                if (k + n == lo[2])
                    shapetypez = 1;
                else if (k + n == lo[2] + 1)
                    shapetypez = 2;
                else if (k + n == hi[2] + 1)
                    shapetypez = 1;
                else if (k + n == hi[2])
                    shapetypez = 4;
                else
                    shapetypez = 3;

                amrex::Real cell_centerx = plo[0] + (i + l) * dx[0];
                amrex::Real cell_centery = plo[1] + (j + m) * dx[1];
                amrex::Real cell_centerz = plo[2] + (k + n) * dx[2];

                amrex::Real rx = (xp[0] - (cell_centerx)) * inv_dx;
                amrex::Real ry = (xp[1] - (cell_centery)) * inv_dy;
                amrex::Real rz = (xp[2] - (cell_centerz)) * inv_dz;

                amrex::Real lval =
                    cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
                amrex::Real mval =
                    cubicspline_1d(periodic[1] ? 3 : shapetypey, ry);
                amrex::Real nval =
                    cubicspline_1d(periodic[2] ? 3 : shapetypez, rz);

                value +=
                    lval * mval * nval * nodaldata(i + l, j + m, k + n, comp);
            }
        }
    }
#endif

    return value;
}

/**
 * @brief Cubic B‑spline interpolation for two components in one stencil pass.
 *
 * Computes the tensor-product cubic B‑spline weight for each stencil node
 * once and accumulates it into vals[0] (comp0) and vals[1] (comp1).
 *
 * @param[in]  xp        Particle position
 * @param[in]  iv        Base cell index
 * @param[in]  min_index Minimum stencil offset per dimension
 * @param[in]  max_index Maximum stencil offset per dimension
 * @param[in]  plo       Domain lower bounds
 * @param[in]  dx        Grid spacing
 * @param[in]  nodaldata Nodal field array
 * @param[in]  comp0     First component index
 * @param[in]  comp1     Second component index
 * @param[in]  lo        Lower physical boundary indices
 * @param[in]  hi        Upper physical boundary indices
 * @param[out] vals      vals[0] = interpolated comp0, vals[1] = interpolated
 * comp1
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cubic_interp_two(amrex::Real xp[AMREX_SPACEDIM],
                 IntVect iv,
                 IntVect min_index,
                 IntVect max_index,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                 amrex::Array4<amrex::Real> nodaldata,
                 int comp0,
                 int comp1,
                 const amrex::GpuArray<int, AMREX_SPACEDIM> &periodic,
                 const int *lo,
                 const int *hi,
                 amrex::Real vals[2])
{
    vals[0] = shunya;
    vals[1] = shunya;

#if (AMREX_SPACEDIM == 1)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    amrex::Real inv_dx = 1.0 / dx[0];
    for (int l = lmin; l < lmax; ++l)
    {
        int shapetypex;
        if (i + l == lo[0] || i + l == hi[0] + 1)
            shapetypex = 1;
        else if (i + l == lo[0] + 1)
            shapetypex = 2;
        else if (i + l == hi[0])
            shapetypex = 4;
        else
            shapetypex = 3;
        amrex::Real rx = (xp[0] - (plo[0] + (i + l) * dx[0])) * inv_dx;
        amrex::Real lval = cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
        vals[0] += lval * nodaldata(i + l, 0, 0, comp0);
        vals[1] += lval * nodaldata(i + l, 0, 0, comp1);
    }

#elif (AMREX_SPACEDIM == 2)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    amrex::Real inv_dx = 1.0 / dx[0];
    amrex::Real inv_dy = 1.0 / dx[1];
    for (int m = mmin; m < mmax; ++m)
    {
        int shapetypey;
        if (j + m == lo[1] || j + m == hi[1] + 1)
            shapetypey = 1;
        else if (j + m == lo[1] + 1)
            shapetypey = 2;
        else if (j + m == hi[1])
            shapetypey = 4;
        else
            shapetypey = 3;
        amrex::Real ry = (xp[1] - (plo[1] + (j + m) * dx[1])) * inv_dy;
        amrex::Real mval = cubicspline_1d(periodic[1] ? 3 : shapetypey, ry);
        for (int l = lmin; l < lmax; ++l)
        {
            int shapetypex;
            if (i + l == lo[0] || i + l == hi[0] + 1)
                shapetypex = 1;
            else if (i + l == lo[0] + 1)
                shapetypex = 2;
            else if (i + l == hi[0])
                shapetypex = 4;
            else
                shapetypex = 3;
            amrex::Real rx = (xp[0] - (plo[0] + (i + l) * dx[0])) * inv_dx;
            amrex::Real lval = cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
            amrex::Real weight = lval * mval;
            vals[0] += weight * nodaldata(i + l, j + m, 0, comp0);
            vals[1] += weight * nodaldata(i + l, j + m, 0, comp1);
        }
    }

#elif (AMREX_SPACEDIM == 3)
    int i = iv[XDIR];
    int lmin = min_index[XDIR];
    int lmax = max_index[XDIR];
    int j = iv[YDIR];
    int mmin = min_index[YDIR];
    int mmax = max_index[YDIR];
    int k = iv[ZDIR];
    int nmin = min_index[ZDIR];
    int nmax = max_index[ZDIR];
    amrex::Real inv_dx = 1.0 / dx[0];
    amrex::Real inv_dy = 1.0 / dx[1];
    amrex::Real inv_dz = 1.0 / dx[2];
    for (int n = nmin; n < nmax; ++n)
    {
        int shapetypez;
        if (k + n == lo[2] || k + n == hi[2] + 1)
            shapetypez = 1;
        else if (k + n == lo[2] + 1)
            shapetypez = 2;
        else if (k + n == hi[2])
            shapetypez = 4;
        else
            shapetypez = 3;
        amrex::Real rz = (xp[2] - (plo[2] + (k + n) * dx[2])) * inv_dz;
        amrex::Real nval = cubicspline_1d(periodic[2] ? 3 : shapetypez, rz);
        for (int m = mmin; m < mmax; ++m)
        {
            int shapetypey;
            if (j + m == lo[1] || j + m == hi[1] + 1)
                shapetypey = 1;
            else if (j + m == lo[1] + 1)
                shapetypey = 2;
            else if (j + m == hi[1])
                shapetypey = 4;
            else
                shapetypey = 3;
            amrex::Real ry = (xp[1] - (plo[1] + (j + m) * dx[1])) * inv_dy;
            amrex::Real mval = cubicspline_1d(periodic[1] ? 3 : shapetypey, ry);
            for (int l = lmin; l < lmax; ++l)
            {
                int shapetypex;
                if (i + l == lo[0] || i + l == hi[0] + 1)
                    shapetypex = 1;
                else if (i + l == lo[0] + 1)
                    shapetypex = 2;
                else if (i + l == hi[0])
                    shapetypex = 4;
                else
                    shapetypex = 3;
                amrex::Real rx = (xp[0] - (plo[0] + (i + l) * dx[0])) * inv_dx;
                amrex::Real lval =
                    cubicspline_1d(periodic[0] ? 3 : shapetypex, rx);
                amrex::Real weight = lval * mval * nval;
                vals[0] += weight * nodaldata(i + l, j + m, k + n, comp0);
                vals[1] += weight * nodaldata(i + l, j + m, k + n, comp1);
            }
        }
    }
#endif
}

/**
 * @brief Evaluates a tensor-product basis function (linear, quadratic, or
 * cubic).
 *
 * For each dimension:
 *   - order 1 → linear basis
 *   - order 2 → quadratic B‑spline
 *   - order 3 → cubic B‑spline
 *
 * Boundary conditions and periodicity are handled per dimension.
 *
 * @param[in] stencil     Stencil offset per dimension
 * @param[in] cell        Base cell index
 * @param[in] xp          Particle position
 * @param[in] plo         Domain lower bounds
 * @param[in] dx          Grid spacing
 * @param[in] order_scheme_directional  Basis order per dimension
 * @param[in] periodic    Periodicity flags per dimension
 * @param[in] lo          Lower boundary indices
 * @param[in] hi          Upper boundary indices
 *
 * @return Basis function value at xp.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
basisval(const IntVect stencil,
         const IntVect cell,
         const amrex::Real xp[AMREX_SPACEDIM],
         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &plo,
         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &dx,
         const GpuArray<int, AMREX_SPACEDIM> &order_scheme_directional,
         const GpuArray<int, AMREX_SPACEDIM> &periodic,
         const int *lo,
         const int *hi)
{
    amrex::Real contrib[AMREX_SPACEDIM];

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        const int i = cell[d];
        const int l = stencil[d];
        const int od = order_scheme_directional[d];

        amrex::Real r = (xp[d] - (plo[d] + i * dx[d])) / dx[d];

        if (od == 1)
        {
            contrib[d] = (l == 0) ? (eka - r) : r;
        }
        else if (od == 2)
        {
            int shapefunctype = 3; // interior default
            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) // → 1
                    + (idx == lo[d] + 1) * (2 - 3)                   // → 2
                    + (idx == hi[d]) * (4 - 3);                      // → 4
            }

            r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];
            contrib[d] = quadraticspline_1d(shapefunctype, r);
        }
        else if (od == 3)
        {
            int shapefunctype = 3; // interior default
            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) // → 1
                    + (idx == lo[d] + 1) * (2 - 3)                   // → 2
                    + (idx == hi[d]) * (4 - 3);                      // → 4
            }

            r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];

            contrib[d] = cubicspline_1d(shapefunctype, r);
        }
    }

    amrex::Real val = eka;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        val *= contrib[d];
    }
    return val;
}

/**
 * @brief Evaluates the derivative of a tensor-product basis function.
 *
 * Computes ∂N/∂x_dir by differentiating the basis in the specified direction
 * and multiplying by the basis values in all other directions.
 *
 * @param[in] dir         Direction of derivative (0=x, 1=y, 2=z)
 * @param[in] stencil     Stencil offset per dimension
 * @param[in] cell        Base cell index
 * @param[in] xp          Particle position
 * @param[in] plo         Domain lower bounds
 * @param[in] dx          Grid spacing
 * @param[in] order_scheme_directional  Basis order per dimension
 * @param[in] periodic    Periodicity flags per dimension
 * @param[in] lo          Lower boundary indices
 * @param[in] hi          Upper boundary indices
 *
 * @return Derivative of basis function in direction dir.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
basisvalder(const int dir,
            const IntVect stencil,
            const IntVect cell,
            const amrex::Real xp[AMREX_SPACEDIM],
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &plo,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &dx,
            const GpuArray<int, AMREX_SPACEDIM> &order_scheme_directional,
            const GpuArray<int, AMREX_SPACEDIM> &periodic,
            const int *lo,
            const int *hi)
{
    amrex::Real prod = eka;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        const int i = cell[d];
        const int l = stencil[d];
        const int od = order_scheme_directional[d];
        const amrex::Real dxinv = eka / dx[d]; // computed once per dimension

        int shapefunctype = 3;

        // Local coordinate r: for linear, centered at cell i; for cubic, at
        // stencil (i+l)
        amrex::Real r;

        amrex::Real val_d = 0.0; // value contribution in dimension d
        amrex::Real der_d = 0.0; // derivative contribution in dimension d

        if (od == 1)
        {
            // Linear shape
            r = (xp[d] - (plo[d] + i * dx[d])) / dx[d];
            val_d = (l == 0) ? (eka - r) : r;
            der_d = (l == 0) ? -dxinv : dxinv;
        }
        else if (od == 2)
        {
            // quadratic spline shape (order 2): classify boundary/periodic type

            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) // → 1
                    + (idx == lo[d] + 1) * (2 - 3)                   // → 2
                    + (idx == hi[d]) * (4 - 3);                      // → 4
            }

            r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];
            val_d = quadraticspline_1d(shapefunctype, r);
            der_d = quadraticspline_1d_der(shapefunctype, r) * dxinv;
        }
        else
        {

            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) // → 1
                    + (idx == lo[d] + 1) * (2 - 3)                   // → 2
                    + (idx == hi[d]) * (4 - 3);                      // → 4
            }

            r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];
            val_d = cubicspline_1d(shapefunctype, r);
            der_d = cubicspline_1d_der(shapefunctype, r) * dxinv;
        }

        // Multiply derivative in the requested direction; value otherwise
        prod *= (d == dir) ? der_d : val_d;
    }

    return prod;
}

/**
 * @brief Evaluates a tensor-product basis function and all its spatial
 * derivatives in a single pass over dimensions.
 *
 * Computes the per-dimension 1D basis value N_d and derivative dN_d once,
 * then assembles:
 *   N        = N_0 * N_1 * ... * N_{SPACEDIM-1}
 *   dN/dx_d  = N_0 * ... * dN_d * ... * N_{SPACEDIM-1}
 *
 * This replaces one call to basisval() plus SPACEDIM calls to basisvalder(),
 * each of which would independently re-evaluate all 1D basis contributions.
 *
 * @param[in]  stencil                    Stencil offset per dimension
 * @param[in]  cell                       Base cell index
 * @param[in]  xp                         Particle position
 * @param[in]  plo                        Domain lower bounds
 * @param[in]  dx                         Grid spacing
 * @param[in]  order_scheme_directional   Basis order per dimension
 * @param[in]  periodic                   Periodicity flags per dimension
 * @param[in]  lo                         Lower boundary indices
 * @param[in]  hi                         Upper boundary indices
 * @param[out] grad                       dN/dx_d for each direction d
 *
 * @return Basis function value N at xp.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
basisval_and_grad(const IntVect stencil,
                  const IntVect cell,
                  const amrex::Real xp[AMREX_SPACEDIM],
                  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &plo,
                  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &dx,
                  const GpuArray<int, AMREX_SPACEDIM> &order_scheme_directional,
                  const GpuArray<int, AMREX_SPACEDIM> &periodic,
                  const int *lo,
                  const int *hi,
                  amrex::Real grad[AMREX_SPACEDIM])
{
    amrex::Real val1d[AMREX_SPACEDIM];
    amrex::Real der1d[AMREX_SPACEDIM];

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        const int i = cell[d];
        const int l = stencil[d];
        const int od = order_scheme_directional[d];
        const amrex::Real dxinv = eka / dx[d];

        if (od == 1)
        {
            amrex::Real r = (xp[d] - (plo[d] + i * dx[d])) / dx[d];
            val1d[d] = (l == 0) ? (eka - r) : r;
            der1d[d] = (l == 0) ? -dxinv : dxinv;
        }
        else if (od == 2)
        {
            int shapefunctype = 3;
            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) +
                    (idx == lo[d] + 1) * (2 - 3) + (idx == hi[d]) * (4 - 3);
            }
            amrex::Real r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];
            val1d[d] = quadraticspline_1d(shapefunctype, r);
            der1d[d] = quadraticspline_1d_der(shapefunctype, r) * dxinv;
        }
        else // od == 3
        {
            int shapefunctype = 3;
            int idx = i + l;
            if (periodic[d] == 0)
            {
                shapefunctype =
                    3 + (idx == lo[d] || idx == hi[d] + 1) * (1 - 3) +
                    (idx == lo[d] + 1) * (2 - 3) + (idx == hi[d]) * (4 - 3);
            }
            amrex::Real r = (xp[d] - (plo[d] + (i + l) * dx[d])) / dx[d];
            val1d[d] = cubicspline_1d(shapefunctype, r);
            der1d[d] = cubicspline_1d_der(shapefunctype, r) * dxinv;
        }
    }

    // Assemble N = product of all 1D basis values
    amrex::Real N = eka;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        N *= val1d[d];

    // Assemble dN/dx_dir: replace val1d[dir] with der1d[dir] in the product
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        amrex::Real g = eka;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            g *= (d == dir) ? der1d[d] : val1d[d];
        grad[dir] = g;
    }

    return N;
}

/**
 * @brief Extracts a symmetric tensor stored in particle rdata into a full
 * matrix.
 *
 * The particle stores only the upper-triangular components of a symmetric
 * tensor in packed order. This routine reconstructs the full tensor by
 * filling both (i,j) and (j,i) entries.
 *
 * @param[in]  p            Particle whose rdata contains the packed tensor
 * @param[in]  start_index  Starting index in p.rdata where the packed tensor
 * begins
 * @param[out] tens         Output full tensor of size (AMREX_SPACEDIM ×
 * AMREX_SPACEDIM)
 *
 * @return None.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
get_tensor(MPMParticleContainer::ParticleType &p,
           int start_index,
           amrex::Real tens[AMREX_SPACEDIM * AMREX_SPACEDIM])
{
    int ind = 0;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int j = i; j < AMREX_SPACEDIM; j++)
        {
            tens[i * AMREX_SPACEDIM + j] = p.rdata(start_index + ind);
            tens[j * AMREX_SPACEDIM + i] = p.rdata(start_index + ind);
            ind++;
        }
    }
}

/**
 * @brief Extracts a vector stored in particle rdata.
 *
 * Reads AMREX_SPACEDIM consecutive entries from p.rdata starting at
 * start_index and writes them into vect[].
 *
 * @param[in]  p            Particle containing vector data
 * @param[in]  start_index  Starting index in p.rdata
 * @param[out] vect         Output vector of length AMREX_SPACEDIM
 *
 * @return None.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
get_vector(MPMParticleContainer::ParticleType &p,
           int start_index,
           amrex::Real vect[AMREX_SPACEDIM])
{
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        vect[i] = p.rdata(start_index + i);
    }
}

/**
 * @brief Updates the deformation gradient tensor F using the velocity gradient.
 *
 * Computes:
 *      Lp = I + (grad v_p) * dt
 *      F_new = Lp * F_old
 *
 * The updated deformation gradient is written back into the particle's rdata.
 *
 * @param[in,out] p            Particle whose deformation gradient is updated
 * @param[in]     start_index  Starting index of F in p.rdata (row-major
 * storage)
 * @param[in]     gradvp       Velocity gradient tensor ∂v/∂x at the particle
 * @param[in]     dt           Time step
 *
 * @return None.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_deformation_gradient_tensor(
    MPMParticleContainer::ParticleType &p,
    int start_index,
    amrex::Real gradvp[AMREX_SPACEDIM][AMREX_SPACEDIM],
    amrex::Real dt)
{
    amrex::Real gradvp_norm = 0.0;
    for (int d1 = 0; d1 < AMREX_SPACEDIM; d1++)
        for (int d2 = 0; d2 < AMREX_SPACEDIM; d2++)
            gradvp_norm += gradvp[d1][d2] * gradvp[d1][d2];

    if (gradvp_norm >= 1.0e10)
    {
        return;
    }

    amrex::Real Lp[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Fp[AMREX_SPACEDIM][AMREX_SPACEDIM];

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            if (i == j)
            {
                Lp[i][j] = 1.0 + gradvp[i][j] * dt;
            }
            else
            {
                Lp[i][j] = gradvp[i][j] * dt;
            }
        }
    }

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            Fp[i][j] = p.rdata(start_index + AMREX_SPACEDIM * i + j);
        }
    }

    // Perform matrix multiplication: Result = Lp * Fp
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
        {
            p.rdata(start_index + i * AMREX_SPACEDIM + j) = 0.0;
            for (int k = 0; k < AMREX_SPACEDIM; ++k)
            {
                p.rdata(start_index + i * AMREX_SPACEDIM + j) +=
                    Lp[i][k] * Fp[k][j];
            }
        }
    }
}

/**
 * @brief Computes the product of a tensor and a vector.
 *
 * Computes:
 *      tensvect[i] = Σ_j tens[i,j] * vect[j]
 *
 * where tens is stored in row-major flattened form.
 *
 * @param[in]  tens      Input tensor of size (AMREX_SPACEDIM × AMREX_SPACEDIM)
 * @param[in]  vect      Input vector of length AMREX_SPACEDIM
 * @param[out] tensvect  Output vector = tens * vect
 *
 * @return None.
 */

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
tensor_vector_pdt(amrex::Real tens[AMREX_SPACEDIM * AMREX_SPACEDIM],
                  amrex::Real vect[AMREX_SPACEDIM],
                  amrex::Real tensvect[AMREX_SPACEDIM])
{
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        tensvect[i] = shunya;
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            tensvect[i] += tens[i * AMREX_SPACEDIM + j] * vect[j];
        }
    }
}

#endif
