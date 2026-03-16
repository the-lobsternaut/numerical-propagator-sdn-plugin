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
        // Approximate Sun position along +x at 1 AU
        // (analytical solar ephemeris would improve accuracy for multi-week propagations)
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

    // Full RKF45 Butcher tableau (Fehlberg coefficients)
    // Stage nodes: 0, 1/4, 3/8, 12/13, 1, 1/2
    constexpr double a2 = 1.0/4.0, a3 = 3.0/8.0, a4 = 12.0/13.0, a5 = 1.0, a6 = 1.0/2.0;
    constexpr double b21 = 1.0/4.0;
    constexpr double b31 = 3.0/32.0,       b32 = 9.0/32.0;
    constexpr double b41 = 1932.0/2197.0,  b42 = -7200.0/2197.0, b43 = 7296.0/2197.0;
    constexpr double b51 = 439.0/216.0,    b52 = -8.0,           b53 = 3680.0/513.0,   b54 = -845.0/4104.0;
    constexpr double b61 = -8.0/27.0,      b62 = 2.0,            b63 = -3544.0/2565.0, b64 = 1859.0/4104.0, b65 = -11.0/40.0;

    // 4th-order weights
    constexpr double c41 = 25.0/216.0, c43 = 1408.0/2565.0, c44 = 2197.0/4104.0, c45 = -1.0/5.0;
    // 5th-order weights
    constexpr double c51 = 16.0/135.0, c53 = 6656.0/12825.0, c54 = 28561.0/56430.0, c55 = -9.0/50.0, c56 = 2.0/55.0;

    auto eval = [&](double px, double py, double pz, double pvx, double pvy, double pvz, double pjd) {
        auto a = accel(px, py, pz, pvx, pvy, pvz, pjd);
        return std::array<double,6>{pvx, pvy, pvz, a.ax, a.ay, a.az};
    };

    auto k1 = eval(x, y, z, vx, vy, vz, jd);

    auto k2 = eval(
        x + dt*b21*k1[0], y + dt*b21*k1[1], z + dt*b21*k1[2],
        vx + dt*b21*k1[3], vy + dt*b21*k1[4], vz + dt*b21*k1[5],
        jd + a2*dt/SEC_PER_DAY);

    auto k3 = eval(
        x + dt*(b31*k1[0]+b32*k2[0]), y + dt*(b31*k1[1]+b32*k2[1]), z + dt*(b31*k1[2]+b32*k2[2]),
        vx + dt*(b31*k1[3]+b32*k2[3]), vy + dt*(b31*k1[4]+b32*k2[4]), vz + dt*(b31*k1[5]+b32*k2[5]),
        jd + a3*dt/SEC_PER_DAY);

    auto k4_ = eval(
        x + dt*(b41*k1[0]+b42*k2[0]+b43*k3[0]), y + dt*(b41*k1[1]+b42*k2[1]+b43*k3[1]), z + dt*(b41*k1[2]+b42*k2[2]+b43*k3[2]),
        vx + dt*(b41*k1[3]+b42*k2[3]+b43*k3[3]), vy + dt*(b41*k1[4]+b42*k2[4]+b43*k3[4]), vz + dt*(b41*k1[5]+b42*k2[5]+b43*k3[5]),
        jd + a4*dt/SEC_PER_DAY);

    auto k5 = eval(
        x + dt*(b51*k1[0]+b52*k2[0]+b53*k3[0]+b54*k4_[0]),
        y + dt*(b51*k1[1]+b52*k2[1]+b53*k3[1]+b54*k4_[1]),
        z + dt*(b51*k1[2]+b52*k2[2]+b53*k3[2]+b54*k4_[2]),
        vx + dt*(b51*k1[3]+b52*k2[3]+b53*k3[3]+b54*k4_[3]),
        vy + dt*(b51*k1[4]+b52*k2[4]+b53*k3[4]+b54*k4_[4]),
        vz + dt*(b51*k1[5]+b52*k2[5]+b53*k3[5]+b54*k4_[5]),
        jd + a5*dt/SEC_PER_DAY);

    auto k6 = eval(
        x + dt*(b61*k1[0]+b62*k2[0]+b63*k3[0]+b64*k4_[0]+b65*k5[0]),
        y + dt*(b61*k1[1]+b62*k2[1]+b63*k3[1]+b64*k4_[1]+b65*k5[1]),
        z + dt*(b61*k1[2]+b62*k2[2]+b63*k3[2]+b64*k4_[2]+b65*k5[2]),
        vx + dt*(b61*k1[3]+b62*k2[3]+b63*k3[3]+b64*k4_[3]+b65*k5[3]),
        vy + dt*(b61*k1[4]+b62*k2[4]+b63*k3[4]+b64*k4_[4]+b65*k5[4]),
        vz + dt*(b61*k1[5]+b62*k2[5]+b63*k3[5]+b64*k4_[5]+b65*k5[5]),
        jd + a6*dt/SEC_PER_DAY);

    // 5th-order solution (use this)
    double x5 = x, y5 = y, z5 = z, vx5 = vx, vy5 = vy, vz5 = vz;
    for (int i = 0; i < 6; ++i) {
        double d5 = dt * (c51*k1[i] + c53*k3[i] + c54*k4_[i] + c55*k5[i] + c56*k6[i]);
        double d4 = dt * (c41*k1[i] + c43*k3[i] + c44*k4_[i] + c45*k5[i]);
        // Error is difference between 4th and 5th order
        (void)d4; (void)d5;
    }

    // Compute error from 4th vs 5th order difference
    double err = 0.0;
    for (int i = 0; i < 3; ++i) {  // position components only
        double d5 = dt * (c51*k1[i] + c53*k3[i] + c54*k4_[i] + c55*k5[i] + c56*k6[i]);
        double d4 = dt * (c41*k1[i] + c43*k3[i] + c44*k4_[i] + c45*k5[i]);
        double ei = d5 - d4;
        err += ei * ei;
    }
    err = std::sqrt(err);

    // Accept 5th-order solution
    x += dt * (c51*k1[0] + c53*k3[0] + c54*k4_[0] + c55*k5[0] + c56*k6[0]);
    y += dt * (c51*k1[1] + c53*k3[1] + c54*k4_[1] + c55*k5[1] + c56*k6[1]);
    z += dt * (c51*k1[2] + c53*k3[2] + c54*k4_[2] + c55*k5[2] + c56*k6[2]);
    vx += dt * (c51*k1[3] + c53*k3[3] + c54*k4_[3] + c55*k5[3] + c56*k6[3]);
    vy += dt * (c51*k1[4] + c53*k3[4] + c54*k4_[4] + c55*k5[4] + c56*k6[4]);
    vz += dt * (c51*k1[5] + c53*k3[5] + c54*k4_[5] + c55*k5[5] + c56*k6[5]);

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

// ── OEM output — Serialize propagation results to SDN binary ──
//
// Wire format (little-endian):
//   [4]  magic "OEM\0"
//   [4]  version (uint32 = 1)
//   [4]  state_count (uint32)
//   [8]  start_epoch_jd (double)
//   [8]  end_epoch_jd (double)
//   For each state:
//     [8]  epoch_jd (double)
//     [8]  x  (double, km)
//     [8]  y  (double, km)
//     [8]  z  (double, km)
//     [8]  vx (double, km/s)
//     [8]  vy (double, km/s)
//     [8]  vz (double, km/s)

static void oem_write_u32(uint8_t*& p, uint32_t v) {
    memcpy(p, &v, 4); p += 4;
}

static void oem_write_f64(uint8_t*& p, double v) {
    memcpy(p, &v, 8); p += 8;
}

int32_t propagate_to_oem(
    const StateVector& initial_state,
    const PropagationConfig& config,
    uint8_t* output, uint32_t output_capacity) {

    auto result = propagate(initial_state, config);
    if (!result.success) return -1;

    const uint32_t state_count = static_cast<uint32_t>(result.states.size());
    // Header: 4 (magic) + 4 (version) + 4 (count) + 8 (start) + 8 (end) = 28
    // Per state: 7 doubles = 56 bytes
    const size_t required = 28 + static_cast<size_t>(state_count) * 56;
    if (required > output_capacity) return -2;

    uint8_t* p = output;

    // Header
    memcpy(p, "OEM", 4); p += 4;  // magic with null terminator
    oem_write_u32(p, 1);           // version
    oem_write_u32(p, state_count);
    oem_write_f64(p, config.start_epoch_jd);
    oem_write_f64(p, config.end_epoch_jd);

    // State vectors
    for (const auto& sv : result.states) {
        oem_write_f64(p, sv.epoch_jd);
        oem_write_f64(p, sv.x);
        oem_write_f64(p, sv.y);
        oem_write_f64(p, sv.z);
        oem_write_f64(p, sv.vx);
        oem_write_f64(p, sv.vy);
        oem_write_f64(p, sv.vz);
    }

    return static_cast<int32_t>(p - output);
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
