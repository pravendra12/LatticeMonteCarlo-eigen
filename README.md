# LatticeMC

LatticeMC is a C++ package for Monte Carlo simulations and analysis of chemical
configurations on crystal lattices. It combines cluster-expansion energy models,
vacancy-mediated kinetic Monte Carlo, equilibrium sampling, simulated annealing,
and tools for studying chemical order and local environments in alloys.

The implementation uses Eigen for numerical operations, MPI and OpenMP for
parallel execution, and an embedded Python interface for parts of the symmetry
and cluster-expansion setup. BCC and FCC supercell generation is available through
the configuration API.

## Capabilities

| Area | Functionality |
| --- | --- |
| Equilibrium sampling | Canonical Monte Carlo using exchanges between lattice sites at fixed composition |
| Kinetic simulation | Vacancy exchanges with predicted migration barriers, event rates, and elapsed physical time |
| Annealing | Temperature-dependent Monte Carlo sampling with the implemented cooling schedule |
| Energy models | Cluster-expansion bases, correlations, symmetry-equivalent clusters, and symmetric cluster-expansion predictions |
| Migration models | Kinetically resolved activation (KRA) and energy-change predictors; local vacancy formation-energy tools |
| Structural analysis | Short-range order, B2 order parameters, B2 cluster analysis, and cluster-to-atom mapping |
| Insertion analysis | Site-resolved energy changes from substituting a specified ghost species |
| Transport data | Vacancy trajectories, collective species displacements, online tracer MSDs, and restart checkpoints |

The executable selects a workflow through a text parameter file. The active
command-line methods are listed below; additional classes and utilities are
available in the source tree.

## Build and dependencies

Use a C++20 compiler with OpenMP support, CMake, Git, an MPI implementation, and
Python development headers/libraries. MPI is a build dependency even when running
a simulation with one process.

| Dependency | Use |
| --- | --- |
| Boost filesystem, system, and iostreams | File handling and compressed input/output |
| Eigen | Vectors, matrices, and numerical operations |
| MPI and OpenMP | Distributed and shared-memory execution |
| pybind11 and Python | Embedded Python interface |
| nlohmann/json | Model parameters and temperature profiles |
| spglib | Crystal symmetry operations |
| MKL or BLAS/LAPACK, when available | Numerical acceleration selected by CMake |

[CMakeLists.txt](CMakeLists.txt) defines dependency discovery and pinned downloads.
It fetches Eigen, nlohmann/json, and spglib, and fetches pybind11 if it cannot find
an installed copy. The first configuration therefore needs network access unless
the dependency sources are already supplied through CMake's FetchContent options.

The Python symmetry helper imports NumPy, ASE, and the Python spglib package.
Install them in the Python environment used by the embedded interpreter:

```bash
python3 -m pip install numpy ase spglib
```

Pandas is optional for postprocessing. The C++ spglib dependency and its Python
package serve different parts of the program; both are used by these workflows.
Keep `lmc/symce/py/` available with the source tree because the embedded interface
locates its Python helper relative to the compiled source path.

From the repository root:

```bash
cmake -S . -B cmake-build -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build --target lmc -j 4
```

The executable is written to `bin/lmc.exe`. For development, use a separate build
directory with `-DCMAKE_BUILD_TYPE=Debug`.

The repository also provides a compiler-selection wrapper:

```bash
bash build.sh release gcc
```

The wrapper accepts `release`/`debug` and `default`/`gcc`/`clang`/`intel`. It deletes
and recreates `cmake-build`, so use the direct CMake commands for incremental builds.

## Prepare and run a simulation

A typical workflow is:

1. Prepare a starting configuration with the desired lattice, composition, and
   vacancy content.
2. Supply a fitted coefficient JSON compatible with that structure and chemistry.
3. Copy the relevant parameter template from [script/](script/) and edit its paths,
   temperatures, neighbor cutoffs, and sampling lengths.
4. Run from a dedicated working directory. Output files and relative input paths
   are resolved from the process's working directory, not the parameter-file directory.
5. Inspect the logs and saved configurations before launching longer calculations.

For example, from the repository root:

```bash
mkdir -p run
cp script/cmc_param.txt run/params.txt
# Edit run/params.txt and provide the configuration and coefficient files.
cd run
../bin/lmc.exe -p params.txt
```

Both `-p` and `--p` select the parameter file. The templates contain example paths
and settings; they are not complete material-specific datasets.

### Available methods

| `simulation_method` | Purpose | Execution |
| --- | --- | --- |
| `CanonicalMcSerial` | Fixed-temperature canonical sampling | One MPI process |
| `SimulatedAnnealing` | Cooling-based configurational sampling | Use one MPI process |
| `KineticMcFirstOmp` | First-order residence-time KMC | One MPI process, OpenMP threads |
| `KineticMcFirstMpi` | First-order residence-time KMC | MPI ranks matched to the vacancy event count |
| `KineticMcChainOmpi` | Second-order chain KMC | MPI ranks matched to the vacancy event count, with OpenMP |
| `Ansys` | Analyze saved configurations and simulation logs | OpenMP frame processing |
| `WidomInsertion` | Ghost-species substitution energies for saved configurations | Snapshot traversal |

These names are the branches currently dispatched by
[lmc/api/src/Home.cpp](lmc/api/src/Home.cpp). A class or deprecated implementation
elsewhere in the repository is not necessarily an available command-line method.

For OpenMP KMC, select `KineticMcFirstOmp` and set the thread count:

```bash
OMP_NUM_THREADS=4 ../bin/lmc.exe -p kmc_params.txt
```

For MPI KMC, the number of ranks must equal the number of first-neighbor vacancy
jump events in the configured neighbor list. For a BCC first-neighbor list with
eight events, for example:

```bash
OMP_NUM_THREADS=1 mpirun -np 8 ../bin/lmc.exe -p kmc_params.txt
```

Select an MPI method in that parameter file. For `KineticMcChainOmpi`, increase
`OMP_NUM_THREADS` according to the threads allocated per rank. Match process and
thread counts to the resources allocated by your workstation or scheduler.

## Input files and parameters

Parameter files contain one key followed by its values on each line. Use spaces
between tokens, lowercase `true`/`false` for booleans, and separate comment lines
beginning with `#`. Write vector values as space-separated numbers, for example
`cutoffs 3 4 5`.

| Parameter | Meaning |
| --- | --- |
| `simulation_method` | Workflow name from the method table |
| `config_filename` | Starting configuration for MC, KMC, or annealing |
| `json_coefficients_filename` | Fitted energy/migration model parameters |
| `structure_type`, `lattice_param` | Reference-lattice settings used during model setup |
| `cutoffs` | Neighbor-shell cutoff distances |
| `temperature` | Fixed temperature for canonical MC or KMC |
| `initial_temperature` | Starting temperature for simulated annealing |
| `maximum_steps` | Terminal simulation step |
| `log_dump_steps` | Logging interval; early adaptive sampling may be denser |
| `config_dump_steps` | Configuration snapshot interval |
| `thermodynamic_averaging_steps` | Energy averaging-window size for MC/KMC |
| `restart_steps`, `restart_energy`, `restart_time` | State counters for restart; physical time applies to KMC |

Use positive dump intervals and, where used, a positive averaging-window size.
Neighbor cutoffs and model settings must agree with the fitted model and lattice.
Temperatures are in kelvin, the energy convention is eV, and CFG lengths are in
Angstrom. KMC times are in seconds with the attempt-frequency convention in
[lmc/constant/include/Constants.hpp](lmc/constant/include/Constants.hpp).

The complete parsed parameter list is in
[Parameter.h](lmc/api/include/Parameter.h) and
[Parameter.cpp](lmc/api/src/Parameter.cpp). Some parameters are specific to one
workflow or retained for older implementations; being parsed does not imply that
every method uses them.

### Configurations and fitted models

The configuration reader dispatches by filename extension:

- `.cfg`, `.cfg.gz`, `.cfg.bz2`
- `.POSCAR`, `.POSCAR.gz`, `.POSCAR.bz2`
- `.xyz`, `.xyz.gz`, `.xyz.bz2`

Use the recognized suffixes, including the uppercase `.POSCAR`. `X` represents a
vacancy. KMC follows a single vacancy exchanging with neighboring atoms. An atom
ID follows the atom as its lattice site changes; atom IDs and lattice IDs are
therefore distinct.

The coefficient JSON must come from a compatible fitted model. Relevant sections
include `symCE` for symmetric cluster-expansion settings and `kra` for KMC
migration settings. Other predictors use their corresponding CE/LVFE parameters.
Coefficient vectors depend on the basis, species, cluster definitions, and ordering;
an arbitrary array of energies cannot substitute for a fitted model. See
[ClusterExpansionParameters.cpp](lmc/ce/src/ClusterExpansionParameters.cpp) for the
fields read by each predictor.

For current simulation entry points, use `config_filename` and leave the legacy
`map_filename` unset.

### Workflow templates

| Template | Settings to adapt |
| --- | --- |
| [cmc_param.txt](script/cmc_param.txt) | Starting configuration, fitted energy model, fixed temperature, and sampling |
| [sa_param.txt](script/sa_param.txt) | Starting configuration, energy model, initial temperature, and run length |
| [kmc_param.txt](script/kmc_param.txt) | Vacancy configuration, energy/migration model, KMC method, temperature, and output intervals |
| [ansys_param.txt](script/ansys_param.txt) | Log/configuration types, frame range, neighbor cutoffs, and analysis flags |
| [widom_insertion_param.txt](script/widom_insertion_param.txt) | Model, starting frame, frame increment, and ghost species |

KMC accepts `time_temperature_filename` for a JSON array of `[time, temperature]`
pairs. Omit it for constant temperature. Simulated annealing uses its own cooling
schedule implemented in [SimulatedAnnealing.cpp](lmc/mc/src/SimulatedAnnealing.cpp).

## Outputs and postprocessing

| Output | Contents |
| --- | --- |
| `cmc_log.txt` | Canonical MC sampling statistics |
| `sa_log.txt` | Annealing step, temperature, and energy statistics |
| `kmc_log.txt` | KMC time, energies, selected-event information, vacancy trajectory, species displacement vectors, and tracer MSDs |
| `<step>.cfg.gz` | Periodically saved configurations |
| `end.cfg.gz` / `end.cfg` | Final configuration for canonical MC / KMC and annealing, respectively |
| `<step>.displacements.gz` | Per-atom KMC displacement checkpoint for tracer restart |
| `ansys_frame_log.txt` | Frame-level structural analysis output |
| `processedConfig/`, `b2ClusterAtomMap/` | Outputs enabled by B2 cluster-analysis options |
| `widomInsertion/<step>.txt.gz` | Site-resolved ghost-species insertion/substitution energy changes |

### Structural and insertion analysis

`Ansys` reads saved frames using `initial_steps` and `increment_steps`. Select
`log_type` from `canonical_mc`, `kinetic_mc`, or `simulated_annealing`, and
`config_type` from `config` or `xyz`. Available switches are `enable_sro`,
`enable_b2_order_param`, `enable_b2_cluster_ansys`, and `enable_b2_cluster_atom_map`.
The cluster-atom mapping option is used with B2 cluster analysis.

Run analysis in the directory containing the matching log and numbered frames.
The current analysis log reader assumes scalar whitespace-separated data and a
single column-header row. KMC vector fields need preprocessing into a compatible
scalar log before using this reader.

`WidomInsertion` currently traverses numbered `.cfg.gz` snapshots. Set
`initial_steps`, a positive `increment_steps`, and `ghost_element` (for example
`X`). The helper follows consecutive available frames at that interval; it does
not launch a Monte Carlo trajectory.

### KMC transport data

Set `log_dump_mode linear` to sample every `log_dump_steps` events. The default
`adaptive` mode retains denser early sampling. Initial, restart, configuration-dump, and final states
are also recorded even when off the regular sampling grid. Equal event intervals
do not imply equal elapsed times; retain the recorded physical time in analysis.

`kmc_log.txt` includes one `dR_<element>` vector per non-vacancy species. This is
the sum of the cumulative unwrapped displacements of that species' atoms.
Each vector occupies one tab-delimited field with space-separated components.
These collective vectors are used for Onsager transport correlations.

The following `MSD_<element>` scalar columns contain `sum_i |dR_i|^2 / N_element`,
where the sum includes every atom of that species, including unmoved atoms, and
excludes the vacancy. Displacements are unwrapped Cartesian vectors relative to
the tracer measurement origin. MSD has squared-distance units, not diffusivity
units. In three dimensions, obtain tracer diffusivity from one sixth of the
long-time linear MSD slope against physical time.

Every exchange replaces the moving atom's old squared-displacement contribution
with its new one. This preserves backtracking: an atom returning to its origin
contributes zero. It is neither a sum of squared jump lengths nor the square of
the collective species vector. Individual vectors remain in memory, but ordinary
log rows only write the species sums divided by fixed species populations.

For example, using optional pandas:

```python
import pandas as pd

transport = pd.read_csv("kmc_log.txt", sep="\t")
time = transport["time"]  # Subtract the tracer origin time for elapsed time.
msd_fe = transport["MSD_Fe"]
```

Full `<step>.displacements.gz` files are written only alongside numbered
configuration checkpoints and the terminal `end.cfg`, using its final step
number. Initial/restart log rows alone do not trigger per-atom output.

Each `<step>.displacements.gz` is independently readable compressed text:

```text
# step 100000
# time 0.025
# displacement_origin_step 0
# displacement_origin_time 0
# segment_start_step 0
atom_id element dx dy dz
0 Fe -2.1 0 0
1 Ni 0 1.4 0
3 Fe 0 0 0
```

Rows are sorted by persistent atom ID, omit the vacancy without renumbering, and
include atoms with zero displacement. Vectors are unwrapped Cartesian
displacements from the recorded measurement origin. They are updated every jump
and streamed directly into gzip with sufficient precision for numeric round trips.
The log and displacement snapshot describe the same state before the next jump.

For example, using optional pandas:

```python
import pandas as pd

df = pd.read_csv("100000.displacements.gz", sep=r"\s+", comment="#", compression="gzip")
df["r_squared"] = df["dx"]**2 + df["dy"]**2 + df["dz"]**2
msd_by_element = df.groupby("element")["r_squared"].mean()
```

This calculates species MSDs relative to the measurement origin. Differences
between matched atom vectors at two snapshots support other time origins.
Collective transport instead uses correlations between species-total vectors.
Read header metadata separately with `gzip.open(path, "rt")`, sort filenames by
numeric step, and process one snapshot at a time to limit memory use.

## Restarting and organizing runs

Keep the configuration, logged state counters, coefficient model, and simulation
settings consistent when restarting. Set `restart_steps` and `restart_energy` to
match the selected configuration, and `restart_time` for KMC. Use the restart
energy in the same reference convention as the original log.

For KMC, supply the vacancy and collective species vectors in the parameter file,
using the values from the matching log row. For example (illustrative vectors):

```text
config_filename 100000.cfg.gz
restart_steps 100000
# Set restart_energy and restart_time to the matching recorded values.
vacancy_trajectory 1.0 2.0 3.0
species_displacement Fe -0.5 -1.0 -1.5
species_displacement Ni -0.5 -1.0 -1.5
# Optional: continue the individual-atom tracer measurement as well.
displacement_restart_filename 100000.displacements.gz
```

Use one `species_displacement <element> dx dy dz` line per non-vacancy species.
Omitted species vectors start at zero. Collective transport only needs these
species totals; no per-atom positions or separate JSON restart files are written.
The vacancy vector and species totals are always taken from the parameters,
independently of the optional atom snapshot.

For continuous tracer measurements, supply the gzip snapshot from the same step
as the configuration. The reader checks step/time, atom IDs/types, finite vectors,
and the tracer measurement origin. It restores the individual displacement vectors
and their original origin, then reconstructs the species squared-displacement sums
once from the vectors. Preserve CFG atom ordering and use the corresponding
configuration: positions are not stored in the displacement snapshot.

Without an atom snapshot, a restart begins a new tracer measurement with zero
individual vectors and zero tracer MSD. This does not reset the collective
species totals supplied in the parameters. Each gzip snapshot records the tracer origin and
`segment_start_step`, the start of the current invocation.

Configuration-dump steps also write a matching KMC log row and atom snapshot,
even when they fall outside the regular logging schedule. The terminal
`end.cfg` matches the atom snapshot named by the final step number.

A same-step displacement snapshot is replaced only after its new gzip stream
finishes. Other existing snapshots are not deleted. When branching from an earlier
configuration, use separate run directories or select files by segment metadata.
KMC logs contain one column header and append only data rows on restart. A log
with an older or different header is rejected before appending; archive the old
log or restart in a new directory when upgrading to the MSD columns. Tracer
origin and segment metadata are stored only in the gzip atom snapshots. The
random-number generator state is not restored, so restarting does not promise an
identical future stochastic trajectory.

[script/restart.py](script/restart.py) uses an older map-file/log convention. Use
the parameter-and-snapshot procedure above for the current KMC output.

## Source layout and development

| Path | Responsibility |
| --- | --- |
| [lmc/api/](lmc/api/) | Parameter parsing, reporting, and workflow construction |
| [lmc/cfg/](lmc/cfg/) | Configurations, lattice/atom mappings, geometry, neighbors, and structure I/O |
| [lmc/ce/](lmc/ce/) | Cluster-expansion bases, cluster enumeration, correlations, and model parameters |
| [lmc/symce/](lmc/symce/) | Symmetric cluster expansion and supporting Python interface |
| [lmc/pred/](lmc/pred/) | Energy, migration-barrier, vacancy-energy, and temperature predictors |
| [lmc/mc/](lmc/mc/) | Monte Carlo algorithms, event selection, sampling, and trajectory output |
| [lmc/ansys/](lmc/ansys/) | Order parameters and cluster analysis |
| [lmc/utility/](lmc/utility/) | Geometry, encoding, training support, and insertion utilities |
| [lmc/constant/](lmc/constant/) | Elements and physical constants |
| [lmc/testing/](lmc/testing/) | Standalone development and diagnostic programs |
| [script/](script/) | Example parameter files and helper scripts |
| [deprecated/](deprecated/) | Older implementations retained for reference |

Public declarations and function documentation are generally in each module's
`include/` directory, with implementations under `src/`. Build the `lmc` target
after changes and validate behavior with an appropriate small configuration.
The current top-level CMake configuration does not enable an automated CTest suite;
programs under `lmc/testing/` may require their own inputs and build setup.

To generate API documentation with Doxygen installed:

```bash
doxygen Doxyfile
```

The supplied configuration writes HTML documentation under `html/`.

## Acknowledgments and references

This package builds on the original lattice Monte Carlo implementation developed by Zhucong Xi and uses ideas and methodology from the `icet` cluster-expansion framework.

* **Original LatticeMonteCarlo code — Zhucong Xi:**
  [GitHub repository](https://github.com/zhucongx/LatticeMonteCarlo)

* **icet — cluster expansion construction and sampling:**
  [Source code](https://gitlab.com/materials-modeling/icet) and [documentation](https://icet.materialsmodeling.org/)

* **icet publication:**
  M. Ångqvist et al., *ICET – A Python Library for Constructing and Sampling Alloy Cluster Expansions* (2019).
  [DOI: 10.1002/adts.201900015](https://doi.org/10.1002/adts.201900015)
  [Open manuscript](https://arxiv.org/abs/1901.08790)


## Troubleshooting

- **Dependency download fails:** check access to the repositories specified in
  CMake, or provide local dependency sources using `FETCHCONTENT_SOURCE_DIR_<NAME>`.
- **Embedded Python import fails:** check the Python interpreter/library selected
  by CMake, availability of NumPy/ASE/spglib in that environment, and the
  `lmc/symce/py/` helper path.
- **MPI method exits with a process-count message:** match the rank count to the
  configured vacancy neighbor-event count, or use `KineticMcFirstOmp` with one rank.
- **Configuration cannot be opened:** check the working directory, filename suffix,
  and paths in the parameter file.
- **Model dimensions or species do not match:** check fitted coefficients, basis,
  chemical species, reference lattice, and neighbor/cluster definitions together.
- **Restart snapshot is rejected:** use the configuration and gzip atom snapshot
  from the same step, with matching counters and preserved atom ordering.

### Tracer MSD regression checks

After configuring CMake with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`, run
`python3 script/test_tracer_msd.py build`. The harness compiles the production
transport code with prescribed events and a predictor stand-in. It checks
backtracking, multiple atoms/species, brute-force MSD agreement, checkpoint
continuation, output columns/cadence, and two MPI ranks. The continuation check
uses identical prescribed jumps, since production RNG state is not checkpointed.
