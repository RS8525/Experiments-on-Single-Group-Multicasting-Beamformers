# Multicast Beamforming Algorithms

Implementation and comparison of multicast beamforming algorithms from the paper:
**"Design of Single-Group Multicasting-Beamformers"**

## Project Structure

```
.
|-- README.md                  # This file
|-- MASTER_HINTS.md           # Development guidelines and design principles
|-- paper/                    # Original paper and reading notes
|-- doc/                      # Theory overview, notation mapping, algorithm glossary
|-- config/                   # Experiment configuration files
|   |-- experiments/          # Specific experiment configs
|   |-- defaults/             # Default parameter configs
|-- src/                      # Source code
|   |-- algorithms/           # Algorithm implementations (no plotting)
|   |   |-- baselines/        # Baseline algorithms
|   |   |-- utils/            # Algortithm subfunctions
|   |-- channels/             # Channel generation functions
|   |-- utils/                # Utility functions
|   |   |-- config/           # Config loading
|   |   |-- experiments/      # Experiment execution helpers
|   |   |-- logging/          # Logging utilities
|   |   |-- metrics/          # SNR and rate metrics
|   |   |-- random/           # Random seed management
|   |   |-- registry/         # Algorithm and channel registries
|   |   |-- validation/       # Constraint checks
|   |-- solvers/              # CVX-based and external solver wrappers
|-- experiments/              # Top-level experiment scripts
|-- plotting/                 # Plotting scripts (load results, generate figures)
|-- results/                  # Experimental results
|   |-- raw/                  # Raw .mat files from experiments
|   |-- processed/            # Post-processed data
|   |-- figures/              # Generated figures
|-- tests/                    # Sanity checks and regression tests
```

## Getting Started

### Prerequisites

- MATLAB R2020b or later
- CVX (for convex optimization algorithms)
- Optimization Toolbox (for SQP-based algorithms)

### Running Experiments

1. **Configure your experiment**: Edit or create a configuration file in `config/experiments/`
   
2. **Run an experiment**:
   ```matlab
   run_experiment('config/experiments/fig4.json')
   ```
   
   Results are saved to `results/raw/<experiment_name>/raw_results.mat`

3. **Generate plots**:
   ```matlab
   cd plotting
   plot_results_example
   ```
   
   Figures are saved to `results/figures/`

## Implemented Algorithms

### Main Algorithms
- **FF-C2**: Fixed-point algorithm (Criterion 2)
- **RC-C2**: Randomized coordinate descent (Criterion 2)
- **SBFC**: Sequential beamformer computation
- **SDR**: Semi-definite relaxation
- **SQP**: Sequential quadratic programming

### Baselines
- **Zero-Forcing**: Zero-forcing beamformer
- **MRT**: Maximum ratio transmission

## Design Principles

- **Pure functions**: All algorithms are implemented as pure functions with clear signatures
- **Config-driven**: Experiment parameters are loaded from JSON/TOML configuration files
- **Reproducible**: Explicit random seed management in all experiments
- **Separation of concerns**: Algorithms in `src/`, experiments in `experiments/`, plotting in `plotting/`
- **lower_snake_case**: For all functions and variables

## Example Usage

```matlab
% Set random seed for reproducibility
set_global_seed(42);

% Generate channel
[H, config] = generate_channel_iid_rayleigh(N, K, M);

% Run algorithm
[W, metrics] = ff_c2(H, config);

% Compute performance
snr = compute_snr(W, H, config);
rate = compute_sum_rate(snr);
```

## Adding New Algorithms

1. Implement the algorithm as a pure function in `src/algorithms/your_algorithm.m`
2. Add function header documentation (inputs, outputs, assumptions)
3. Add tests in `tests/test_your_algorithm.m`
4. Update this README with algorithm description

## Contributing

- Follow naming conventions (lower_snake_case)
- Add header comments to all functions
- Write tests for new functionality
- Keep algorithms separate from plotting code

## References

- Paper: "Design of Single-Group Multicasting-Beamformers"
- See `paper/` directory for the original paper and notes


