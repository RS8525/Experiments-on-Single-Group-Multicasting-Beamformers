# Project TODOs

## Priority: Critical (Must-Have for Basic Functionality)

### 1. Algorithm Registry Updates
- [ ] **Update `resolve_algorithm.m`** to include all implemented algorithms:
  - ✅ `sdr_beamformer` - Already registered
  - ✅ `ff_c2` - Already registered
  - ✅ `rc_c2` - Already registered
  - ✅ `sbfc` - Already registered
  - ❌ `rc_c2_it_update` - **NOT registered** (registry has 'rc_c2_update' but file is 'rc_c2_it_update')
  - **Action**: Fix name mismatch in registry or rename algorithm
Done

### 2. Example Configuration File
- [ ] **Fix `example_experiment.json`** - Current structure doesn't match expected format
  - Current: Has nested "configs" object with single algorithm
  - Expected: Should follow fig4.json structure with `scenarios`, `algorithms` array, `k_list`
  - **Action**: Rewrite to match the run_raw_experiment.m expectations

### 3. Baseline Algorithms
- [ ] **Implement baseline algorithms** in `src/algorithms/baselines/`
  - README mentions: Zero-Forcing (ZF), Maximum Ratio Transmission (MRT)
  - Currently: Folder is empty
  - **Action**: 
    - Implement `mrt_beamformer.m`
    - Implement `zf_beamformer.m`
    - Register in `resolve_algorithm.m`

### 4. Results Saving Enhancement
- [ ] **Add CSV export functionality** for Python plotting compatibility
  - **Current**: Results saved only as `.mat` files
  - **Recommended**: Add optional CSV export in `process_raw_results.m`
  - **Benefits**:
    - ✅ Easy integration with Python (pandas, matplotlib, seaborn)
    - ✅ Human-readable format for quick inspection
    - ✅ Version control friendly (text-based)
    - ✅ Cross-platform compatibility
  - **Proposed structure**:
    ```
    results/processed/<experiment_name>/
      ├── summary_power.csv        # K vs power for all methods
      ├── summary_min_snr.csv      # K vs min SNR for all methods
      ├── summary_feasibility.csv  # K vs feasibility rate for all methods
      └── metadata.json            # Experiment metadata
    ```

---

## Priority: High (Important for Usability)

### 5. Documentation Updates
- [ ] **Update README.md**:
  - Remove SQP from "Implemented Algorithms" (not found in codebase)
  - Update actual implemented algorithms list
  - Add section on rc_c2_it_update algorithm
  - Fix paths and examples

### 6. Configuration Validation
- [ ] **Add config schema validation** in `load_config.m`
  - Check required fields: `experiment_name`, `seed`, `num_trials`, `k_list`, `scenarios`, `algorithms`
  - Validate `gamma_db` or `gamma_linear` is present
  - Provide helpful error messages for missing/invalid fields

### 7. Plotting Scripts Enhancement
- [ ] **Improve `plot_example_results.m`**:
  - Add automatic scenario selection menu if multiple scenarios exist
  - Add method comparison plots (all methods on same figure)
  - Add error bars/confidence intervals
  - Export plots in multiple formats (PNG, PDF, SVG)

- [ ] **Enhance `plotting_processed.m`**:
  - Add more plot types (CDF, boxplots, heatmaps)
  - Add automatic color scheme selection
  - Improve legend formatting

### 8. Process Results Automation
- [ ] **Create wrapper script** `process_and_plot.m`:
  - Automatically process raw results after experiments
  - Generate standard plots
  - Create summary report
  - Example usage: `process_and_plot('fig4')`

---

## Priority: Medium (Nice-to-Have)

### 9. Testing Infrastructure
- [ ] **Add test scripts** in `tests/`:
  - Unit tests for parameter generation utilities
  - Integration tests for algorithm pipeline
  - Regression tests comparing against known results
  - Example: `test_algorithms.m`, `test_parameter_generation.m`

### 10. Performance Metrics
- [ ] **Add computational complexity tracking**:
  - Track FLOPs per algorithm
  - Measure actual runtime
  - Memory usage profiling
  - Add to metrics output structure

### 11. Algorithm Comparison
- [ ] **Create comparison config** `config/experiments/algorithm_comparison.json`:
  - Run all algorithms on same channel realizations
  - Fair comparison with same parameters
  - Generate comparison plots automatically

### 12. Channel Models
- [ ] **Implement additional channel models**:
  - Rician fading
  - Correlated channels
  - Line-of-sight (LoS) component
  - Register in `resolve_channel.m`

---

## Priority: Low (Future Enhancements)

### 13. Python Integration (Optional)
- [ ] **Create Python plotting package** `plotting_py/`:
  - `load_results.py` - Load MATLAB .mat or CSV files
  - `plot_results.py` - Generate publication-quality plots
  - `process_results.py` - Post-processing in Python
  - Requirements: `pandas`, `matplotlib`, `seaborn`, `scipy`

### 14. Batch Experiments
- [ ] **Create batch runner** `run_batch_experiments.m`:
  - Run multiple configs in sequence
  - Parallel processing support
  - Progress tracking and ETA
  - Email notification on completion

### 15. Interactive Visualization
- [ ] **Create MATLAB App** for interactive results exploration:
  - Load and browse multiple experiments
  - Interactive parameter exploration
  - Real-time plot generation
  - Export selected results

### 16. Documentation Website
- [ ] **Generate documentation** using MATLAB's publish feature:
  - Algorithm descriptions with equations
  - Usage examples
  - API reference
  - Tutorial notebooks

---

## Clarifications Needed

### Q1: CSV Export Format
**Question**: What level of detail should CSV exports include?

**Options**:
- **A) Summary only** (mean, median, std per K value) - Lightweight, easier to plot
- **B) Full trial data** (all trials for all K) - Complete data, large files
- **C) Both** (separate files for summary and raw data)

**Recommendation**: **Option C** - Provide flexibility for both quick plotting and detailed analysis

### Q2: Baseline Algorithm Priority
**Question**: Which baseline algorithms are most important?

**Options**:
- MRT (Maximum Ratio Transmission)
- ZF (Zero Forcing) 
- RZF (Regularized Zero Forcing)
- Equal power allocation
- Random beamforming

**Recommendation**: Start with **MRT** and **ZF** as they're mentioned in README

### Q3: Python vs MATLAB Plotting
**Question**: Should we prioritize Python or MATLAB plotting?

**Considerations**:
- MATLAB: Already integrated, familiar, good for rapid prototyping
- Python: Better for publication-quality plots, more flexible, larger ecosystem

**Recommendation**: 
- Keep MATLAB plotting for quick checks during experiments
- Add Python plotting for publication-quality figures
- Support both through CSV export

### Q4: Algorithm Naming Convention
**Question**: Should we standardize algorithm names?

**Current inconsistency**:
- File: `rc_c2_it_update.m`
- Function: `rc_c2_it_update`
- Registry: `rc_c2_update` ❌

**Recommendation**: Use consistent naming across all:
- File name = Function name = Registry name
- Decide on: `rc_c2_it_update` or `rc_c2_update`

---

## Quick Start Checklist

To get the project running end-to-end:

1. ✅ Fix algorithm registry (Task #1)
2. ✅ Fix example_experiment.json (Task #2)
3. ✅ Run a simple experiment: `run_experiment('config/experiments/fig4.json')`
4. ✅ Process results: `process_raw_results('results/raw/fig4/raw_results.mat', 'results/processed/fig4/processed_results.mat')`
5. ✅ Generate plots: `cd plotting; plot_example_results`

---

## Dependency Check

**MATLAB Toolboxes Required**:
- ✅ Optimization Toolbox (confirmed needed for `eigs`, optimization functions)
- ✅ Statistics and Machine Learning Toolbox (for percentile calculations)
- ⚠️ CVX (external) - Required for `sdr_beamformer.m`

**Action**: Add toolbox dependency check in `run_experiment.m` startup

---

## Notes

- **Parameter Generation**: ✅ Already refactored and working (Option 2 implementation complete)
- **Algorithm Interface**: ✅ All algorithms follow consistent `[W, metrics] = algo(H, config)` pattern
- **Config System**: ✅ JSON-based configuration system is functional
- **Results Pipeline**: ⚠️ Partially functional (raw → processed works, but needs CSV export)
- **Plotting**: ⚠️ Functional but needs enhancement for multiple methods comparison

---

## Current Project Status

### ✅ Working Components
- Core algorithm implementations (ff_c2, rc_c2, sbfc, sdr_beamformer, rc_c2_it_update)
- Parameter generation utilities (gamma, noise_power)
- Channel generation (iid_rayleigh)
- Raw experiment execution
- Results saving (.mat format)
- Basic plotting

### ⚠️ Needs Attention
- Algorithm registry name mismatch
- Example config file format
- CSV export for Python compatibility
- Baseline algorithms missing
- Documentation updates

### ❌ Not Yet Implemented
- Baseline algorithms (MRT, ZF)
- Python plotting package
- Automated testing
- Batch experiment runner
- Additional channel models

---

**Last Updated**: January 24, 2026
**Version**: 1.0
