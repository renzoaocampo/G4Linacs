# G4Linacs

Simulations of linear accelerators (Linacs) and clinical heads in Geant4.

## Description

G4Linacs is a collection of scenarios and configurations for simulating radiation sources and medical linear-accelerator heads using Geant4. The goal is to provide reusable models (e.g., `Halcyon6MV`, `Infinity6MV`, `Infinity6MV_PSF`) that allow dose-distribution studies, detector-response profiling, and the generation of phase-space files (PSF) and other parameters relevant to medical physics.

Each subfolder contains a variant of the model (source code, CMake build files, macros, and output data) designed for different linacs or simulation configurations.

## Repository Structure

* `Halcyon6MV/` — Halcyon 6 MV model
* `Infinity6MV/` — Infinity 6 MV model
* `Infinity6MV_PSF/` — Variant for PSF generation
* Each folder contains:

  * `CMakeLists.txt` — build file
  * `src/` — C++ source (DetectorConstruction, PrimaryGenerator, PhysicsList, etc.)
  * `include/` — headers
  * `macros/` — example macros to run simulations
  * `build/` — project and compilation outputs (should not be versioned if using `.gitignore`)
  * `data.ipynb`, `merged_output.csv` — example notebooks and data

## Requirements

* Geant4 (compatible version; 10.x or newer recommended)
* CMake
* Modern C++ compiler (MSVC, GCC, or Clang)
* Python (optional, for notebooks and analysis)

## Build (Windows)

Open the **Developer Command Prompt for Visual Studio**, create a build directory inside the project folder, and run:

```powershell
cd Halcyon6MV    # or the desired subproject
mkdir build; cd build
cmake ..
cmake --build . --config Release
```

To select a specific generator or Visual Studio version, use
`cmake .. -G "Visual Studio 16 2019" -A x64` instead of plain `cmake ..`.

On Linux or macOS, use the system generator and a standard terminal (e.g., `cmake .. && make`).

## Run

After building, run the simulation executable with:

```
sim --t <num_threads> --s <seed> [macro_filename]
```

* `--t <num_threads>`: number of threads. If omitted, all CPU cores are used.
* `--s <seed>`: random seed. If omitted, a random seed is selected.
* `[macro_filename]`: path to a Geant4 macro file (e.g., `macros/default.mac`). If omitted, the default macro is used.

Examples:

```powershell
.\sim.exe --t 8 --s 12345 macros\default.mac
.\sim.exe                 # runs with default macro, random seed, all threads
```

Notebooks (`data.ipynb`) are provided for post-processing and output analysis.

## Best Practices

* Do not version the `build/` folder (add an entry in `.gitignore`).
* Keep detector code (`DetectorConstruction`) and physics list (`PhysicsList`) separated to simplify testing.
* Specify units explicitly (`cm`, `mm`, `m`, etc.) when defining dimensions.

## Contributions

To contribute:

* Open an issue describing the improvement or bug.
* Create a branch for your work.
* Submit a pull request with tests or usage examples.

## License

Add your preferred license here (MIT, Apache-2.0, etc.) or leave it under your group/institution policy.

---

If you want, I can:

* Add a recommended `.gitignore` for Geant4/CMake projects on Windows.
* Generate a basic `CONTRIBUTING.md`.
* Format the `README.md` files inside each subfolder with model-specific instructions.
