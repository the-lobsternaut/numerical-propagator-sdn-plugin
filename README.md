# Numerical Propagator SDN Plugin

High-fidelity numerical orbit propagator for the [Space Data Network](https://github.com/the-lobsternaut/space-data-network), powered by [Tudat](https://github.com/tudat-team/tudat). Configurable force models, multiple integrators, and optional state transition matrix propagation — compiled to WebAssembly.

## Force Models

| Model | Options | Default |
|---|---|---|
| **Gravity** | Point-mass, J2, J4, full spherical harmonics (degree/order) | SH 8×8 |
| **Atmosphere** | None, Exponential, NRLMSISE-00 | Exponential |
| **SRP** | None, Cannonball | Cannonball |
| **Third-body** | Sun, Moon, Jupiter (toggleable) | Sun + Moon |
| **Relativistic** | Schwarzschild, Lense-Thirring | Off |

## Integrators

| Type | Step | Use Case |
|---|---|---|
| Euler | Fixed | Testing only |
| RK4 | Fixed | Fast, moderate accuracy |
| RKF45 | Variable | General purpose |
| **RKF78** | Variable | **Default — high accuracy** |
| RKDP87 | Variable | Near-machine precision |
| Bulirsch-Stoer | Variable | Stiff problems, long arcs |
| ABM | Multi-step | Long-duration propagation |

## Usage

```cpp
#include "numerical_prop/propagator.h"

using namespace numerical_prop;

// Initial state (km, km/s, J2000 ECI)
StateVector initial = { epoch_jd, 6778.0, 0, 0, 0, 7.669, 0 };

// Configure
PropagationConfig config;
config.duration_days = 7.0;
config.output_step_s = 60.0;
config.forces = presets::leo_standard();
config.integrator.type = IntegratorType::RKF78;
config.integrator.rel_tolerance = 1e-12;

// Propagate
auto result = propagate(initial, config);
// result.states: vector of StateVector
// result.computation_time_ms, result.function_evaluations

// Single-point propagation
auto state = propagate_to_epoch(initial, target_jd);

// OEM FlatBuffer output
uint8_t buffer[1024*1024];
int32_t size = propagate_to_oem(initial, config, buffer, sizeof(buffer));
```

## Presets

```cpp
// Ready-made force model configurations
auto leo   = presets::leo_standard();      // SH 8×8, exp atm, SRP, Sun+Moon
auto geo   = presets::geo_standard();      // SH 8×8, no drag, SRP, Sun+Moon
auto deep  = presets::deep_space();        // Point mass, Sun+Moon+Jupiter
auto hifi  = presets::leo_high_fidelity(); // SH 70×70, NRLMSISE, relativistic
```

## STM Propagation

State Transition Matrix propagation for covariance mapping:

```cpp
config.propagate_stm = true;
auto result = propagate(initial, config);
// result.stm_history: 6×6 matrix at each output time
```

## Data Flow

```
Initial state (km, km/s)
        │
        ▼
┌─────────────────┐
│ ForceModelConfig │ ← presets or custom
│ IntegratorConfig │
└────────┬────────┘
         │
         ▼
  Tudat propagation engine
  (SI internally, km/s at I/O)
         │
         ▼
  ┌──────────────┐
  │ OEM $OEM     │ ← FlatBuffer binary
  │ StateVector  │ ← C++ structs
  │ STM history  │ ← optional 6×6 matrices
  └──────────────┘
```

## Building

```bash
cd src/cpp && mkdir -p build && cd build

# Native (requires Tudat source)
cmake .. -DTUDAT_DIR=/path/to/tudat
make -j4
./test_propagator

# WASM (via emsdk)
source /path/to/emsdk/emsdk_env.sh
emcmake cmake .. -DTUDAT_DIR=/path/to/tudat-wasm
emmake make -j4
```

## Dependencies

- [Tudat](https://github.com/tudat-team/tudat) — TU Delft Astrodynamics Toolbox (BSD-3-Clause)
- [spacedatastandards.org](https://spacedatastandards.org) — OEM FlatBuffers schema

## WASM

- Export: `NUMERICAL_PROP_WASM`
- Initial memory: 64 MB
- Max memory: 512 MB
- Supports pthreads for integrator parallelism

## Units

- **Interface**: km, km/s (J2000 ECI)
- **Internal**: SI (m, m/s, kg) — matches Tudat convention
- **Time**: Julian Date (TDB)

## License

BSD-3-Clause (matches Tudat)
