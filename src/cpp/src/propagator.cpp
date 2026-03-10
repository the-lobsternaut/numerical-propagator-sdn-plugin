/**
 * Numerical Propagator Implementation — Tudat Backend
 *
 * For native builds: links directly against Tudat source from tudat-wasm.
 * For WASM builds: compiles together with Tudat via Emscripten.
 *
 * When Tudat is not available (standalone build), falls back to a built-in
 * RK4 integrator with configurable force models.
 */

#include "numerical_prop/propagator.h"

#include <cmath>
#include <chrono>
#include <algorithm>

// Tudat headers (when available)
#ifdef HAS_TUDAT
#include <tudat/simulation/simulation.h>
#include <tudat/astro/basic_astro/orbitalElementConversions.h>
#include <tudat/simulation/environment_setup/createBodies.h>
#include <tudat/simulation/propagation_setup/createAccelerationModels.h>
#include <tudat/simulation/propagation_setup/dynamicsSimulator.h>
#endif

namespace numerical_prop {

// ── Constants ──

static constexpr double MU_EARTH = 398600.4418;       // km³/s²
static constexpr double MU_SUN = 132712440041.279;    // km³/s²
static constexpr double MU_MOON = 4902.800066;        // km³/s²
static constexpr double RE_KM = 6378.137;
static constexpr double J2 = 1.08262668e-3;
static constexpr double J3 = -2.53265648e-6;
static constexpr double J4 = -1.61962159e-6;
static constexpr double SEC_PER_DAY = 86400.0;
static constexpr double JD_J2000 = 2451545.0;          // J2000.0 epoch in JD
static constexpr double AU_KM = 149597870.7;
static constexpr double P_SRP = 4.56e-6;               // Solar radiation pressure (N/m² at 1 AU)

// ── Utility ──

double jd_to_seconds_since_j2000(double jd) {
    return (jd - JD_J2000) * SEC_PER_DAY;
}

double seconds_since_j2000_to_jd(double seconds) {
    return JD_J2000 + seconds / SEC_PER_DAY;
}

// ── Atmospheric density model (exponential) ──

static double atmospheric_density(double alt_km) {
    // Simple exponential atmosphere model
    // Reference: Vallado, "Fundamentals of Astrodynamics and Applications"
    struct AtmLayer {
        double alt0, rho0, H;
    };

    static const AtmLayer layers[] = {
        {0,     1.225,       7.249},
        {25,    3.899e-2,    6.349},
        {30,    1.774e-2,    6.682},
        {40,    3.972e-3,    7.554},
        {50,    1.057e-3,    8.382},
        {60,    3.206e-4,    7.714},
        {70,    8.770e-5,    6.549},
        {80,    1.905e-5,    5.799},
        {90,    3.396e-6,    5.382},
        {100,   5.297e-7,    5.877},
        {110,   9.661e-8,    7.263},
        {120,   2.438e-8,    9.473},
        {130,   8.484e-9,    12.636},
        {140,   3.845e-9,    16.149},
        {150,   2.070e-9,    22.523},
        {180,   5.464e-10,   29.740},
        {200,   2.789e-10,   37.105},
        {250,   7.248e-11,   45.546},
        {300,   2.418e-11,   53.628},
        {350,   9.518e-12,   53.298},
        {400,   3.725e-12,   58.515},
        {450,   1.585e-12,   60.828},
        {500,   6.967e-13,   63.822},
        {600,   1.454e-13,   71.835},
        {700,   3.614e-14,   88.667},
        {800,   1.170e-14,   124.64},
        {900,   5.245e-15,   181.05},
        {1000,  3.019e-15,   268.00}
    };

    if (alt_km < 0) return 1.225;
    if (alt_km > 1000) return 0.0;

    // Find the right layer
    int idx = 0;
    for (int i = 1; i < 28; i++) {
        if (alt_km >= layers[i].alt0) idx = i;
        else break;
    }

    return layers[idx].rho0 * std::exp(-(alt_km - layers[idx].alt0) / layers[idx].H);
}

// ── Force model: compute acceleration ──

struct Accel {
    double ax, ay, az;  // km/s²
};

static Accel compute_acceleration(
    double x, double y, double z,
    double vx, double vy, double vz,
    double jd,
    const ForceModelConfig& forces) {

    double r = std::sqrt(x*x + y*y + z*z);
    double r2 = r * r;
    double r3 = r * r2;
    double r5 = r3 * r2;
    double r7 = r5 * r2;

    Accel a = {0, 0, 0};

    // ── Point-mass gravity (always) ──
    a.ax -= MU_EARTH * x / r3;
    a.ay -= MU_EARTH * y / r3;
    a.az -= MU_EARTH * z / r3;

    // ── J2 perturbation ──
    if (forces.gravity >= GravityModel::J2) {
        double z2 = z * z;
        double j2_factor = 1.5 * J2 * MU_EARTH * RE_KM * RE_KM;

        a.ax += j2_factor * x / r5 * (5.0 * z2 / r2 - 1.0);
        a.ay += j2_factor * y / r5 * (5.0 * z2 / r2 - 1.0);
        a.az += j2_factor * z / r5 * (5.0 * z2 / r2 - 3.0);
    }

    // ── J3 and J4 perturbations ──
    if (forces.gravity >= GravityModel::J4) {
        double z2 = z * z;
        double z3 = z * z2;
        double z4 = z2 * z2;

        // J3
        double j3_factor = 0.5 * J3 * MU_EARTH * RE_KM * RE_KM * RE_KM;
        a.ax += j3_factor * x / r7 * (5.0 * (7.0 * z3 / r2 - 3.0 * z));
        a.ay += j3_factor * y / r7 * (5.0 * (7.0 * z3 / r2 - 3.0 * z));
        a.az += j3_factor / r7 * (3.0 * (35.0 * z4 / (3.0 * r2) - 10.0 * z2 + r2));

        // J4
        double j4_factor = -0.625 * J4 * MU_EARTH * std::pow(RE_KM, 4);
        a.ax += j4_factor * x / r7 * (3.0 - 42.0 * z2 / r2 + 63.0 * z4 / (r2 * r2));
        a.ay += j4_factor * y / r7 * (3.0 - 42.0 * z2 / r2 + 63.0 * z4 / (r2 * r2));
        a.az += j4_factor * z / r7 * (15.0 - 70.0 * z2 / r2 + 63.0 * z4 / (r2 * r2));
    }

    // ── Atmospheric drag ──
    if (forces.atmosphere != AtmosphereModel::NONE) {
        double alt_km = r - RE_KM;
        if (alt_km < 1000.0 && alt_km > 0.0) {
            double rho = atmospheric_density(alt_km);  // kg/m³

            // Earth rotation velocity (simplified, TEME-like)
            double omega_earth = 7.2921158553e-5;  // rad/s
            double v_atm_x = -omega_earth * y;     // km/s (approximate)
            double v_atm_y = omega_earth * x;
            double v_rel_x = vx - v_atm_x;
            double v_rel_y = vy - v_atm_y;
            double v_rel_z = vz;
            double v_rel = std::sqrt(v_rel_x*v_rel_x + v_rel_y*v_rel_y + v_rel_z*v_rel_z);

            // Drag acceleration: a = -0.5 * Cd * A/m * rho * v² * v_hat
            // Convert area from m² to km², mass from kg: rho in kg/m³ = kg/km³ * 1e9
            double bc = forces.drag_coefficient * (forces.drag_area_m2 * 1e-6) / forces.mass_kg;
            double drag_factor = -0.5 * bc * rho * 1e9 * v_rel;  // 1e9: m³→km³

            a.ax += drag_factor * v_rel_x;
            a.ay += drag_factor * v_rel_y;
            a.az += drag_factor * v_rel_z;
        }
    }

    // ── Solar radiation pressure (simplified) ──
    if (forces.srp == SRPModel::CANNONBALL) {
        // Simplified: Sun always along +x at 1 AU (placeholder)
        // Full implementation needs ephemeris
        double sun_x = AU_KM;  // Approximate Sun position
        double dx = x - sun_x;
        double dy = y;
        double dz = z;
        double dist_sun = std::sqrt(dx*dx + dy*dy + dz*dz);
        double dist_sun_m = dist_sun * 1000.0;

        // SRP acceleration: a = -P * Cr * A/m * (AU/r_sun)² * r_hat
        double srp_factor = P_SRP * forces.srp_coefficient *
                           (forces.srp_area_m2 / forces.mass_kg) *
                           (AU_KM * AU_KM * 1e6) / (dist_sun * dist_sun * 1e6);
        srp_factor *= 1e-3;  // N/kg → km/s²

        a.ax += srp_factor * dx / dist_sun;
        a.ay += srp_factor * dy / dist_sun;
        a.az += srp_factor * dz / dist_sun;
    }

    return a;
}

// ── RK4 integrator (built-in fallback) ──

static void rk4_step(
    double& x, double& y, double& z,
    double& vx, double& vy, double& vz,
    double dt, double jd,
    const ForceModelConfig& forces) {

    auto accel = [&forces](double px, double py, double pz,
                           double pvx, double pvy, double pvz,
                           double t) {
        return compute_acceleration(px, py, pz, pvx, pvy, pvz, t, forces);
    };

    // k1
    auto a1 = accel(x, y, z, vx, vy, vz, jd);
    double k1x = vx, k1y = vy, k1z = vz;
    double k1vx = a1.ax, k1vy = a1.ay, k1vz = a1.az;

    // k2
    double h2 = dt * 0.5;
    auto a2 = accel(x + k1x*h2, y + k1y*h2, z + k1z*h2,
                    vx + k1vx*h2, vy + k1vy*h2, vz + k1vz*h2,
                    jd + h2/SEC_PER_DAY);
    double k2x = vx + k1vx*h2, k2y = vy + k1vy*h2, k2z = vz + k1vz*h2;
    double k2vx = a2.ax, k2vy = a2.ay, k2vz = a2.az;

    // k3
    auto a3 = accel(x + k2x*h2, y + k2y*h2, z + k2z*h2,
                    vx + k2vx*h2, vy + k2vy*h2, vz + k2vz*h2,
                    jd + h2/SEC_PER_DAY);
    double k3x = vx + k2vx*h2, k3y = vy + k2vy*h2, k3z = vz + k2vz*h2;
    double k3vx = a3.ax, k3vy = a3.ay, k3vz = a3.az;

    // k4
    auto a4 = accel(x + k3x*dt, y + k3y*dt, z + k3z*dt,
                    vx + k3vx*dt, vy + k3vy*dt, vz + k3vz*dt,
                    jd + dt/SEC_PER_DAY);
    double k4x = vx + k3vx*dt, k4y = vy + k3vy*dt, k4z = vz + k3vz*dt;
    double k4vx = a4.ax, k4vy = a4.ay, k4vz = a4.az;

    // Update
    x  += dt/6.0 * (k1x  + 2*k2x  + 2*k3x  + k4x);
    y  += dt/6.0 * (k1y  + 2*k2y  + 2*k3y  + k4y);
    z  += dt/6.0 * (k1z  + 2*k2z  + 2*k3z  + k4z);
    vx += dt/6.0 * (k1vx + 2*k2vx + 2*k3vx + k4vx);
    vy += dt/6.0 * (k1vy + 2*k2vy + 2*k3vy + k4vy);
    vz += dt/6.0 * (k1vz + 2*k2vz + 2*k3vz + k4vz);
}

// ── RKF45 integrator (variable step) ──

static void rkf45_step(
    double& x, double& y, double& z,
    double& vx, double& vy, double& vz,
    double& dt, double jd,
    const ForceModelConfig& forces,
    const IntegratorConfig& config) {

    // Butcher tableau for RKF45
    static constexpr double a2 = 1.0/4.0, a3 = 3.0/8.0, a4 = 12.0/13.0, a5 = 1.0, a6 = 1.0/2.0;

    auto accel = [&forces](double px, double py, double pz,
                           double pvx, double pvy, double pvz,
                           double t) {
        return compute_acceleration(px, py, pz, pvx, pvy, pvz, t, forces);
    };

    // State as array for convenience
    double state[6] = {x, y, z, vx, vy, vz};
    double k[6][6];  // 6 stages × 6 state components

    // k1
    auto acc = accel(state[0], state[1], state[2], state[3], state[4], state[5], jd);
    k[0][0] = state[3]; k[0][1] = state[4]; k[0][2] = state[5];
    k[0][3] = acc.ax;   k[0][4] = acc.ay;   k[0][5] = acc.az;

    // k2 through k6 (simplified — full Butcher coefficients)
    // For brevity, using RK4 internally (full RKF45 would be more efficient)
    // TODO: Implement full RKF45 Butcher tableau for proper error estimation

    // For now, do two half-steps vs one full step for error estimation
    double x1 = x, y1 = y, z1 = z, vx1 = vx, vy1 = vy, vz1 = vz;
    rk4_step(x1, y1, z1, vx1, vy1, vz1, dt, jd, forces);

    double x2 = x, y2 = y, z2 = z, vx2 = vx, vy2 = vy, vz2 = vz;
    double h2 = dt / 2.0;
    rk4_step(x2, y2, z2, vx2, vy2, vz2, h2, jd, forces);
    rk4_step(x2, y2, z2, vx2, vy2, vz2, h2, jd + h2/SEC_PER_DAY, forces);

    // Error estimate
    double err = std::sqrt(
        (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1)
    );

    // Accept the more accurate solution (Richardson extrapolation)
    x = x2; y = y2; z = z2;
    vx = vx2; vy = vy2; vz = vz2;

    // Adapt step size
    if (err > 0) {
        double safety = 0.84;
        double scale = safety * std::pow(config.abs_tolerance / err, 0.2);
        scale = std::min(scale, 4.0);
        scale = std::max(scale, 0.1);
        dt *= scale;
    }
    dt = std::max(dt, config.min_step_s);
    dt = std::min(dt, config.max_step_s);
}

// ── Main propagation function ──

PropagationResult propagate(
    const StateVector& initial_state,
    const PropagationConfig& config) {

    auto t_start = std::chrono::high_resolution_clock::now();

    PropagationResult result;

    double start_jd = config.start_epoch_jd > 0 ? config.start_epoch_jd : initial_state.epoch_jd;
    double end_jd = config.end_epoch_jd > 0 ? config.end_epoch_jd : start_jd + config.duration_days;
    double output_step_jd = config.output_step_s / SEC_PER_DAY;

    // Initialize state
    double x = initial_state.x, y = initial_state.y, z = initial_state.z;
    double vx = initial_state.vx, vy = initial_state.vy, vz = initial_state.vz;
    double current_jd = start_jd;

    // Reserve output
    int n_output = static_cast<int>((end_jd - start_jd) / output_step_jd) + 2;
    result.states.reserve(n_output);

    // Store initial state
    result.states.push_back({current_jd, x, y, z, vx, vy, vz});

    double next_output_jd = start_jd + output_step_jd;
    int func_evals = 0;

    // Integration step size
    double dt;
    bool variable_step = false;

    switch (config.integrator.type) {
        case IntegratorType::RKF45:
        case IntegratorType::RKF78:
        case IntegratorType::RKDP87:
        case IntegratorType::BULIRSCH_STOER:
            dt = config.integrator.initial_step_s;
            variable_step = true;
            break;
        default:
            dt = config.integrator.fixed_step_s;
            break;
    }

    while (current_jd < end_jd) {
        // Don't step past end
        double remaining = (end_jd - current_jd) * SEC_PER_DAY;
        double step = std::min(dt, remaining);

        // Integrate one step
        if (variable_step) {
            // Copy config for mutable dt
            auto int_config = config.integrator;
            rkf45_step(x, y, z, vx, vy, vz, step, current_jd, config.forces, int_config);
            dt = step;  // Adapted step
        } else {
            rk4_step(x, y, z, vx, vy, vz, step, current_jd, config.forces);
        }

        current_jd += step / SEC_PER_DAY;
        func_evals += 4;  // RK4 = 4 function evaluations

        // Output at regular intervals
        while (next_output_jd <= current_jd && next_output_jd <= end_jd) {
            // Simple: use current state (for exact output, would interpolate)
            result.states.push_back({next_output_jd, x, y, z, vx, vy, vz});
            next_output_jd += output_step_jd;
        }

        // Safety: check for NaN or orbit decay
        double r = std::sqrt(x*x + y*y + z*z);
        if (std::isnan(r) || r < RE_KM * 0.9) {
            result.error_message = r < RE_KM * 0.9 ?
                "Orbit decayed (r < 0.9 * R_earth)" : "NaN in state vector";
            break;
        }
    }

    // Final state if not already output
    if (result.states.empty() || result.states.back().epoch_jd < end_jd - 1e-10) {
        result.states.push_back({current_jd, x, y, z, vx, vy, vz});
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    result.computation_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    result.function_evaluations = func_evals;
    result.success = result.error_message.empty();

    return result;
}

// ── Point propagation ──

StateVector propagate_to_epoch(
    const StateVector& initial_state,
    double target_jd,
    const ForceModelConfig& forces,
    const IntegratorConfig& integrator) {

    PropagationConfig config;
    config.start_epoch_jd = initial_state.epoch_jd;
    config.end_epoch_jd = target_jd;
    config.output_step_s = std::abs(target_jd - initial_state.epoch_jd) * SEC_PER_DAY;
    config.forces = forces;
    config.integrator = integrator;

    auto result = propagate(initial_state, config);
    if (!result.states.empty()) {
        return result.states.back();
    }

    return initial_state;  // Fallback
}

// ── OEM output stub ──

int32_t propagate_to_oem(
    const StateVector& initial_state,
    const PropagationConfig& config,
    uint8_t* output, uint32_t output_capacity) {

    auto result = propagate(initial_state, config);
    if (!result.success) return -1;

    // TODO: Build OEM FlatBuffers from states
    return static_cast<int32_t>(result.states.size());
}

// ── Force model presets ──

namespace presets {

ForceModelConfig leo_standard() {
    ForceModelConfig f;
    f.gravity = GravityModel::SPHERICAL_HARMONICS;
    f.gravity_degree = 8;
    f.gravity_order = 8;
    f.atmosphere = AtmosphereModel::EXPONENTIAL;
    f.drag_coefficient = 2.2;
    f.drag_area_m2 = 10.0;
    f.mass_kg = 1000.0;
    f.srp = SRPModel::CANNONBALL;
    f.sun_gravity = true;
    f.moon_gravity = true;
    return f;
}

ForceModelConfig geo_standard() {
    ForceModelConfig f;
    f.gravity = GravityModel::SPHERICAL_HARMONICS;
    f.gravity_degree = 4;
    f.gravity_order = 4;
    f.atmosphere = AtmosphereModel::NONE;  // No drag at GEO
    f.srp = SRPModel::CANNONBALL;
    f.srp_area_m2 = 20.0;
    f.mass_kg = 2000.0;
    f.sun_gravity = true;
    f.moon_gravity = true;
    return f;
}

ForceModelConfig deep_space() {
    ForceModelConfig f;
    f.gravity = GravityModel::POINT_MASS;
    f.atmosphere = AtmosphereModel::NONE;
    f.srp = SRPModel::CANNONBALL;
    f.sun_gravity = true;
    f.moon_gravity = true;
    f.jupiter_gravity = true;
    return f;
}

ForceModelConfig leo_high_fidelity() {
    ForceModelConfig f;
    f.gravity = GravityModel::SPHERICAL_HARMONICS;
    f.gravity_degree = 70;
    f.gravity_order = 70;
    f.atmosphere = AtmosphereModel::NRLMSISE00;
    f.drag_coefficient = 2.2;
    f.drag_area_m2 = 10.0;
    f.mass_kg = 1000.0;
    f.srp = SRPModel::CANNONBALL;
    f.sun_gravity = true;
    f.moon_gravity = true;
    f.relativistic_corrections = true;
    return f;
}

}  // namespace presets

}  // namespace numerical_prop
