# Multicast Beamforming Algorithms

Implementation and comparison of multicast beamforming algorithms from the paper:
**"Design of Single-Group Multicasting-Beamformers"**

## Project Structure

```
.
|-- README.md                  # This file
|-- MASTER_HINTS.md           # Development guidelines and design principles
|-- paper/                    # Original paper and reading notes
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
|   |   |-- parameters/       # Parameter generation
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

3. **Process the results**: Compute summary statistics from raw experiment data
   ```matlab
   process_raw_results('results/raw/fig4/raw_results.mat', ...
                       'results/processed/fig4/processed_results.mat')
   ```
   
   This computes mean, median, percentiles (p25, p75), min, and max for all metrics.

4. **Generate plots**:
   
   **Option A - Quick plotting (single method/scenario)**:
   ```matlab
   cd plotting
   plot_example_results  % Interactive - automatically finds your results
   ```
   
   **Option B - Multi-method comparison plots**:
   ```matlab
   % Plot mean power vs K for all methods
   plotting_processed('results/processed/fig4/processed_results.mat', ...
                      'mean', 'power', ...
                      'SavePath', 'results/figures/power_comparison.png')
   
   % Plot mean minimum SNR vs K
   plotting_processed('results/processed/fig4/processed_results.mat', ...
                      'mean', 'min_snr', ...
                      'SavePath', 'results/figures/min_snr_comparison.png')
   
   % Plot feasibility rate
   plotting_processed('results/processed/fig4/processed_results.mat', ...
                      'mean', 'feasible', ...
                      'SavePath', 'results/figures/feasibility_comparison.png')
   ```
   
   Figures are saved to `results/figures/`

### Available Metrics for Plotting

- **`power`** - Total transmit power
- **`min_snr`** - Minimum SNR across all users
- **`snr_mean`** - Average SNR across all users
- **`feasible`** - Feasibility rate (0 = infeasible, 1 = feasible)
- **`solve_time`** - Algorithm computation time

### Available Statistics

- **`mean`** - Mean across trials
- **`median`** - Median across trials
- **`p25`** - 25th percentile
- **`p75`** - 75th percentile
- **`min`** - Minimum value
- **`max`** - Maximum value

## Implemented Algorithms

### Main Algorithms
- **FF-C2** (`ff_c2`): Full Featured Combine-2 algorithm
- **RC-C2** (`rc_c2`): Reduced Complexity Combine-2 with heuristic
- **RC-C2-IT** (`rc_c2_it_update`): RC-C2 with iterative SNR increasing updates (Algorithm 2)
- **SBFC** (`sbfc`): Sequential beamformer computation
- **SDR** (`sdr_beamformer`): Semi-definite relaxation with randomization

### Baselines
- **Zero-Forcing**: Zero-forcing beamformer (to be implemented)
- **MRT**: Maximum ratio transmission (to be implemented)

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


