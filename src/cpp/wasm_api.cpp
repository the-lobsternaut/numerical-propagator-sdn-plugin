/**
 * numerical-propagator-sdn-plugin WASM API
 *
 * High-fidelity orbit propagation with configurable force models.
 * JSON-in/JSON-out for propagation requests and results.
 */

#include "numerical_prop/propagator.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;
#endif

#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

// ============================================================================
// Version
// ============================================================================

static std::string version() {
    return "0.1.0";
}

// ============================================================================
// JSON-in/JSON-out wrappers
// ============================================================================

/// Propagate from JSON config → JSON result
/// Input JSON: {
///   "state": {"epoch_jd":..., "x":..., "y":..., "z":..., "vx":..., "vy":..., "vz":...},
///   "config": {"duration_days":7, "output_step_s":60, ...}
/// }
/// Output JSON: {"success":true, "num_states":N, "computation_time_ms":..., "states":[...]}
static std::string propagate_json(const std::string& input_json) {
    // In WASM context, the JSON parsing is handled by a lightweight parser.
    // This is a structural placeholder — the actual JSON parse uses the
    // project's JSON utility (nlohmann or flatbuffers JSON).
    // For the SDN plugin ABI, the primary path is the C-API below.
    return "{\"error\":\"use C-API propagate() for binary I/O\"}";
}

/// Propagate to single epoch — lightweight point query
static std::string propagate_to_epoch_json(double epoch_jd,
                                            double x, double y, double z,
                                            double vx, double vy, double vz,
                                            double target_jd) {
    numerical_prop::StateVector sv;
    sv.epoch_jd = epoch_jd;
    sv.x = x; sv.y = y; sv.z = z;
    sv.vx = vx; sv.vy = vy; sv.vz = vz;

    auto result = numerical_prop::propagate_to_epoch(sv, target_jd);

    std::ostringstream oss;
    oss.precision(15);
    oss << "{\"epoch_jd\":" << result.epoch_jd
        << ",\"x\":" << result.x << ",\"y\":" << result.y << ",\"z\":" << result.z
        << ",\"vx\":" << result.vx << ",\"vy\":" << result.vy << ",\"vz\":" << result.vz
        << "}";
    return oss.str();
}

/// Get preset force model config as JSON
static std::string get_preset_json(const std::string& preset_name) {
    numerical_prop::ForceModelConfig fm;
    if (preset_name == "leo_standard") fm = numerical_prop::presets::leo_standard();
    else if (preset_name == "geo_standard") fm = numerical_prop::presets::geo_standard();
    else if (preset_name == "deep_space") fm = numerical_prop::presets::deep_space();
    else if (preset_name == "leo_high_fidelity") fm = numerical_prop::presets::leo_high_fidelity();
    else return "{\"error\":\"unknown preset\"}";

    std::ostringstream oss;
    oss << "{\"gravity_degree\":" << fm.gravity_degree
        << ",\"gravity_order\":" << fm.gravity_order
        << ",\"drag_coefficient\":" << fm.drag_coefficient
        << ",\"drag_area_m2\":" << fm.drag_area_m2
        << ",\"mass_kg\":" << fm.mass_kg
        << ",\"srp_coefficient\":" << fm.srp_coefficient
        << ",\"srp_area_m2\":" << fm.srp_area_m2
        << ",\"sun_gravity\":" << (fm.sun_gravity ? "true" : "false")
        << ",\"moon_gravity\":" << (fm.moon_gravity ? "true" : "false")
        << ",\"jupiter_gravity\":" << (fm.jupiter_gravity ? "true" : "false")
        << ",\"relativistic\":" << (fm.relativistic_corrections ? "true" : "false")
        << "}";
    return oss.str();
}

/// Utility: JD ↔ seconds since J2000
static double jd_to_j2000(double jd) {
    return numerical_prop::jd_to_seconds_since_j2000(jd);
}

static double j2000_to_jd(double seconds) {
    return numerical_prop::seconds_since_j2000_to_jd(seconds);
}

// ============================================================================
// SDN Plugin ABI (C exports)
// ============================================================================

#ifdef __EMSCRIPTEN__

extern "C" {

EMSCRIPTEN_KEEPALIVE
void* sdn_malloc(size_t size) {
    return malloc(size);
}

EMSCRIPTEN_KEEPALIVE
void sdn_free(void* ptr) {
    free(ptr);
}

/// Propagate and output OEM FlatBuffers binary
EMSCRIPTEN_KEEPALIVE
int32_t propagate_to_oem(const uint8_t* state_buf, size_t state_len,
                         double duration_days, double step_s,
                         uint8_t* output, size_t output_len) {
    try {
        if (state_len < 56) return -3; // need 7 doubles
        const double* sv = reinterpret_cast<const double*>(state_buf);
        numerical_prop::StateVector initial;
        initial.epoch_jd = sv[0];
        initial.x = sv[1]; initial.y = sv[2]; initial.z = sv[3];
        initial.vx = sv[4]; initial.vy = sv[5]; initial.vz = sv[6];

        numerical_prop::PropagationConfig config;
        config.duration_days = duration_days;
        config.output_step_s = step_s;

        return numerical_prop::propagate_to_oem(
            initial, config, output, static_cast<uint32_t>(output_len));
    } catch (...) {
        return -1;
    }
}

} // extern "C"

// ============================================================================
// Embind Exports
// ============================================================================

EMSCRIPTEN_BINDINGS(sdn_numerical_propagator) {
    function("version", &version);
    function("propagateToEpoch", &propagate_to_epoch_json);
    function("getPreset", &get_preset_json);
    function("jdToJ2000", &jd_to_j2000);
    function("j2000ToJd", &j2000_to_jd);
}

#endif
