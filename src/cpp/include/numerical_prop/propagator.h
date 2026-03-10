#ifndef NUMERICAL_PROPAGATOR_PLUGIN_H
#define NUMERICAL_PROPAGATOR_PLUGIN_H

/**
 * Numerical Propagator SDN Plugin — Powered by Tudat
 *
 * High-fidelity orbit propagation with configurable force models:
 *   - Point-mass gravity (two-body)
 *   - Spherical harmonics (J2, J4, full degree/order)
 *   - Atmospheric drag (exponential, NRLMSISE-00)
 *   - Solar radiation pressure (cannonball)
 *   - Third-body perturbations (Sun, Moon, planets)
 *   - Relativistic corrections (Schwarzschild, Lense-Thirring)
 *
 * Integrators:
 *   - Fixed-step: RK4, Euler
 *   - Variable-step: RKF45, RKF78, RKDP87
 *   - Extrapolation: Bulirsch-Stoer
 *   - Multi-step: Adams-Bashforth-Moulton
 *
 * Two build modes:
 *   1. Native C++ — links against tudat-wasm source directly
 *   2. WASM — compiles with Emscripten alongside tudat
 *
 * Units: SI (meters, seconds, kg) internally, km/s for SDN interface
 */

#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace numerical_prop {

// ─── State Vector ───

struct StateVector {
    double epoch_jd;       // Julian date (TDB)
    double x, y, z;        // km (J2000 ECI)
    double vx, vy, vz;     // km/s
};

// ─── Force Model Configuration ───

enum class GravityModel {
    POINT_MASS,            // Two-body only
    J2,                    // Zonal J2
    J4,                    // Zonal J2+J3+J4
    SPHERICAL_HARMONICS    // Full degree/order (configurable)
};

enum class AtmosphereModel {
    NONE,
    EXPONENTIAL,           // Simple exponential density
    NRLMSISE00             // NRLMSISE-00 empirical model
};

enum class SRPModel {
    NONE,
    CANNONBALL             // Simple cannonball SRP
};

struct ForceModelConfig {
    // Gravity
    GravityModel gravity = GravityModel::SPHERICAL_HARMONICS;
    int gravity_degree = 8;    // Spherical harmonics degree
    int gravity_order = 8;     // Spherical harmonics order

    // Atmosphere
    AtmosphereModel atmosphere = AtmosphereModel::EXPONENTIAL;
    double drag_coefficient = 2.2;
    double drag_area_m2 = 10.0;     // Cross-section area (m²)
    double mass_kg = 1000.0;         // Spacecraft mass (kg)

    // Solar radiation pressure
    SRPModel srp = SRPModel::CANNONBALL;
    double srp_coefficient = 1.2;    // Radiation pressure coefficient
    double srp_area_m2 = 10.0;      // SRP area (m²)

    // Third-body
    bool sun_gravity = true;
    bool moon_gravity = true;
    bool jupiter_gravity = false;

    // Relativistic
    bool relativistic_corrections = false;
};

// ─── Integrator Configuration ───

enum class IntegratorType {
    EULER,
    RK4,                   // Fixed-step Runge-Kutta 4
    RKF45,                 // Variable-step Runge-Kutta-Fehlberg 4(5)
    RKF78,                 // Variable-step Runge-Kutta-Fehlberg 7(8)
    RKDP87,                // Variable-step Dormand-Prince 8(7)
    BULIRSCH_STOER,        // Bulirsch-Stoer extrapolation
    ABM                    // Adams-Bashforth-Moulton multistep
};

struct IntegratorConfig {
    IntegratorType type = IntegratorType::RKF78;

    // Fixed-step
    double fixed_step_s = 30.0;      // Step size for fixed-step integrators

    // Variable-step
    double initial_step_s = 30.0;    // Initial step size
    double min_step_s = 0.001;       // Minimum step size
    double max_step_s = 3600.0;      // Maximum step size
    double rel_tolerance = 1e-12;    // Relative tolerance
    double abs_tolerance = 1e-12;    // Absolute tolerance
};

// ─── Propagation Configuration ───

struct PropagationConfig {
    double start_epoch_jd = 0.0;     // Start epoch (JD TDB)
    double end_epoch_jd = 0.0;       // End epoch (JD TDB)
    double duration_days = 7.0;      // Duration if end_epoch not set
    double output_step_s = 60.0;     // Output interpolation step

    ForceModelConfig forces;
    IntegratorConfig integrator;

    bool propagate_stm = false;      // Propagate State Transition Matrix
};

// ─── Propagation Result ───

struct PropagationResult {
    std::vector<StateVector> states;
    bool success = false;
    std::string error_message;
    double computation_time_ms = 0.0;
    int function_evaluations = 0;

    // Optional STM at each output time (6×6 row-major)
    std::vector<std::array<double, 36>> stm_history;
};

// ─── Core API ───

/**
 * Propagate from initial Cartesian state.
 * State in km + km/s (J2000 ECI). Internally converts to SI for Tudat.
 */
PropagationResult propagate(
    const StateVector& initial_state,
    const PropagationConfig& config);

/**
 * Propagate to a single epoch (point propagation).
 * Implements the Propagator interface for plugin composition.
 */
StateVector propagate_to_epoch(
    const StateVector& initial_state,
    double target_jd,
    const ForceModelConfig& forces = {},
    const IntegratorConfig& integrator = {});

/**
 * Propagate and output as OEM FlatBuffers binary.
 */
int32_t propagate_to_oem(
    const StateVector& initial_state,
    const PropagationConfig& config,
    uint8_t* output, uint32_t output_capacity);

// ─── Presets ───

namespace presets {

/// LEO satellite with standard force model
ForceModelConfig leo_standard();

/// GEO satellite (no drag, SRP important)
ForceModelConfig geo_standard();

/// Deep space (point mass + Sun/Moon/Jupiter)
ForceModelConfig deep_space();

/// High-fidelity LEO (full 70×70 gravity, NRLMSISE, SRP, relativistic)
ForceModelConfig leo_high_fidelity();

}  // namespace presets

// ─── Utility ───

/// Convert JD to seconds since J2000 (Tudat time system)
double jd_to_seconds_since_j2000(double jd);

/// Convert seconds since J2000 to JD
double seconds_since_j2000_to_jd(double seconds);

}  // namespace numerical_prop

#endif  // NUMERICAL_PROPAGATOR_PLUGIN_H
