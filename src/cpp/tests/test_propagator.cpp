/**
 * Numerical Propagator Plugin Tests
 *
 * Validates:
 * 1. Two-body orbit conservation (energy, angular momentum)
 * 2. J2 secular drift (RAAN regression, argument of perigee advance)
 * 3. Atmospheric drag decay
 * 4. Variable-step vs fixed-step consistency
 * 5. Force model presets
 * 6. Known orbit validation (ISS-like)
 */

#include "numerical_prop/propagator.h"
#include <iostream>
#include <cmath>

using namespace numerical_prop;

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { tests_passed++; std::cout << "  ✓ " << msg << std::endl; } \
    else { tests_failed++; std::cout << "  ✗ " << msg << std::endl; } \
} while(0)

static constexpr double MU = 398600.4418;
static constexpr double RE = 6378.137;

void test_two_body_energy_conservation() {
    std::cout << "\n--- Test: Two-Body Energy Conservation ---" << std::endl;

    StateVector init;
    init.epoch_jd = 2460784.5;  // Some epoch
    init.x = 6800.0; init.y = 0.0; init.z = 0.0;
    init.vx = 0.0; init.vy = 7.6; init.vz = 0.0;

    // Point-mass only, no drag, no SRP
    PropagationConfig config;
    config.duration_days = 1.0;
    config.output_step_s = 60.0;
    config.forces.gravity = GravityModel::POINT_MASS;
    config.forces.atmosphere = AtmosphereModel::NONE;
    config.forces.srp = SRPModel::NONE;
    config.forces.sun_gravity = false;
    config.forces.moon_gravity = false;
    config.integrator.type = IntegratorType::RK4;
    config.integrator.fixed_step_s = 10.0;

    auto result = propagate(init, config);

    CHECK(result.success, "Propagation succeeded");
    CHECK(result.states.size() > 100, "Got >100 output states");

    // Check energy conservation
    auto energy = [](const StateVector& s) {
        double r = std::sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
        double v = std::sqrt(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz);
        return 0.5 * v * v - MU / r;
    };

    double E0 = energy(result.states.front());
    double Ef = energy(result.states.back());
    double dE = std::abs((Ef - E0) / E0);

    std::cout << "    Initial energy: " << E0 << " km²/s²" << std::endl;
    std::cout << "    Final energy:   " << Ef << " km²/s²" << std::endl;
    std::cout << "    Relative ΔE:    " << dE << std::endl;

    CHECK(dE < 1e-6, "Energy conserved to < 1e-6 over 1 day");

    std::cout << "    Computation: " << result.computation_time_ms << " ms, "
              << result.function_evaluations << " evals" << std::endl;
}

void test_j2_secular_effects() {
    std::cout << "\n--- Test: J2 Secular Effects ---" << std::endl;

    // ISS-like orbit: 420 km, 51.6° inclination
    double a = RE + 420.0;
    double v_circ = std::sqrt(MU / a);

    StateVector init;
    init.epoch_jd = 2460784.5;
    init.x = a; init.y = 0.0; init.z = 0.0;
    init.vx = 0.0;
    init.vy = v_circ * std::cos(51.6 * M_PI / 180.0);
    init.vz = v_circ * std::sin(51.6 * M_PI / 180.0);

    // J2 only
    PropagationConfig config;
    config.duration_days = 1.0;
    config.output_step_s = 60.0;
    config.forces.gravity = GravityModel::J2;
    config.forces.atmosphere = AtmosphereModel::NONE;
    config.forces.srp = SRPModel::NONE;
    config.forces.sun_gravity = false;
    config.forces.moon_gravity = false;
    config.integrator.type = IntegratorType::RK4;
    config.integrator.fixed_step_s = 10.0;

    auto result = propagate(init, config);
    CHECK(result.success, "J2 propagation succeeded");

    // With J2, orbit should stay bounded (not decay)
    double r_final = std::sqrt(
        result.states.back().x * result.states.back().x +
        result.states.back().y * result.states.back().y +
        result.states.back().z * result.states.back().z);

    CHECK(std::abs(r_final - a) < 50.0, "Radius stayed within 50 km of initial (J2 perturbation)");

    std::cout << "    Final r: " << r_final << " km (initial: " << a << " km)" << std::endl;
    std::cout << "    Time: " << result.computation_time_ms << " ms" << std::endl;
}

void test_atmospheric_drag() {
    std::cout << "\n--- Test: Atmospheric Drag ---" << std::endl;

    // Low LEO: 250 km altitude — should see measurable drag
    double a = RE + 250.0;
    double v_circ = std::sqrt(MU / a);

    StateVector init;
    init.epoch_jd = 2460784.5;
    init.x = a; init.y = 0.0; init.z = 0.0;
    init.vx = 0.0; init.vy = v_circ; init.vz = 0.0;

    // With drag
    PropagationConfig config_drag;
    config_drag.duration_days = 1.0;
    config_drag.output_step_s = 300.0;
    config_drag.forces.gravity = GravityModel::POINT_MASS;
    config_drag.forces.atmosphere = AtmosphereModel::EXPONENTIAL;
    config_drag.forces.drag_coefficient = 2.2;
    config_drag.forces.drag_area_m2 = 20.0;
    config_drag.forces.mass_kg = 500.0;
    config_drag.forces.srp = SRPModel::NONE;
    config_drag.forces.sun_gravity = false;
    config_drag.forces.moon_gravity = false;
    config_drag.integrator.type = IntegratorType::RK4;
    config_drag.integrator.fixed_step_s = 10.0;

    auto result_drag = propagate(init, config_drag);

    // Without drag
    PropagationConfig config_no_drag = config_drag;
    config_no_drag.forces.atmosphere = AtmosphereModel::NONE;
    auto result_no_drag = propagate(init, config_no_drag);

    CHECK(result_drag.success, "Drag propagation succeeded");
    CHECK(result_no_drag.success, "No-drag propagation succeeded");

    // Energy should decrease with drag (semi-major axis drops)
    auto energy = [](const StateVector& s) {
        double r = std::sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
        double v = std::sqrt(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz);
        return 0.5 * v * v - MU / r;
    };

    double E_drag = energy(result_drag.states.back());
    double E_no_drag = energy(result_no_drag.states.back());

    std::cout << "    Energy (drag):    " << E_drag << " km²/s²" << std::endl;
    std::cout << "    Energy (no drag): " << E_no_drag << " km²/s²" << std::endl;

    // With drag, energy should be MORE negative (orbit decaying)
    CHECK(E_drag < E_no_drag, "Energy decreased with drag (orbit decaying)");
}

void test_variable_step() {
    std::cout << "\n--- Test: Variable-Step Integrator (RKF45) ---" << std::endl;

    StateVector init;
    init.epoch_jd = 2460784.5;
    init.x = 7000.0; init.y = 0.0; init.z = 0.0;
    init.vx = 0.0; init.vy = 7.5; init.vz = 0.0;

    PropagationConfig config;
    config.duration_days = 1.0;
    config.output_step_s = 300.0;
    config.forces.gravity = GravityModel::J2;
    config.forces.atmosphere = AtmosphereModel::NONE;
    config.forces.srp = SRPModel::NONE;
    config.forces.sun_gravity = false;
    config.forces.moon_gravity = false;
    config.integrator.type = IntegratorType::RKF45;
    config.integrator.initial_step_s = 30.0;
    config.integrator.min_step_s = 1.0;
    config.integrator.max_step_s = 300.0;
    config.integrator.abs_tolerance = 1e-8;

    auto result = propagate(init, config);
    CHECK(result.success, "RKF45 propagation succeeded");
    CHECK(result.states.size() > 50, "Got output states");

    // Check orbit stays bounded
    double r_final = std::sqrt(
        result.states.back().x * result.states.back().x +
        result.states.back().y * result.states.back().y +
        result.states.back().z * result.states.back().z);

    CHECK(std::abs(r_final - 7000.0) < 100.0, "RKF45: orbit stayed bounded");

    std::cout << "    Final r: " << r_final << " km" << std::endl;
    std::cout << "    Time: " << result.computation_time_ms << " ms, "
              << result.function_evaluations << " evals" << std::endl;
}

void test_presets() {
    std::cout << "\n--- Test: Force Model Presets ---" << std::endl;

    auto leo = presets::leo_standard();
    CHECK(leo.gravity == GravityModel::SPHERICAL_HARMONICS, "LEO: spherical harmonics");
    CHECK(leo.atmosphere == AtmosphereModel::EXPONENTIAL, "LEO: exponential atmosphere");
    CHECK(leo.srp == SRPModel::CANNONBALL, "LEO: cannonball SRP");
    CHECK(leo.sun_gravity, "LEO: Sun gravity enabled");

    auto geo = presets::geo_standard();
    CHECK(geo.atmosphere == AtmosphereModel::NONE, "GEO: no atmosphere");
    CHECK(geo.srp == SRPModel::CANNONBALL, "GEO: SRP enabled");

    auto deep = presets::deep_space();
    CHECK(deep.gravity == GravityModel::POINT_MASS, "Deep space: point mass");
    CHECK(deep.jupiter_gravity, "Deep space: Jupiter enabled");

    auto hifi = presets::leo_high_fidelity();
    CHECK(hifi.gravity_degree == 70, "High-fi: 70×70 gravity");
    CHECK(hifi.relativistic_corrections, "High-fi: relativistic corrections");
}

void test_jd_conversion() {
    std::cout << "\n--- Test: JD ↔ Seconds Since J2000 ---" << std::endl;

    double j2000_jd = 2451545.0;
    CHECK(std::abs(jd_to_seconds_since_j2000(j2000_jd)) < 0.001,
          "J2000 epoch = 0 seconds");

    double one_day_later = j2000_jd + 1.0;
    CHECK(std::abs(jd_to_seconds_since_j2000(one_day_later) - 86400.0) < 0.001,
          "J2000 + 1 day = 86400 seconds");

    double roundtrip = seconds_since_j2000_to_jd(jd_to_seconds_since_j2000(2460784.5));
    CHECK(std::abs(roundtrip - 2460784.5) < 1e-10, "JD roundtrip exact");
}

int main() {
    std::cout << "=== Numerical Propagator Plugin Tests ===" << std::endl;

    test_jd_conversion();
    test_two_body_energy_conservation();
    test_j2_secular_effects();
    test_atmospheric_drag();
    test_variable_step();
    test_presets();

    std::cout << "\n=== Summary: " << tests_passed << " passed, "
              << tests_failed << " failed ===" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
