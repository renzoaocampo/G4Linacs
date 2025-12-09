# G4Linacs

Simulations of medical linear accelerators (Linacs) and clinical heads using Geant4.

## Description

G4Linacs is a collection of simulation scenarios and configurations for modeling radiation sources and medical linear accelerator heads with Geant4. The project aims to provide reusable models (for example: `Halcyon6MV`, `Infinity6MV`, `Infinity6MV_PSF`) to study dose distributions, characterize detector responses, compute point spread functions (PSFs), and extract other parameters relevant to medical physics.

Each subfolder contains a variant of the model (source code, CMake build files, example macros, and sample outputs) tailored for different linacs or simulation setups.

## Repository structure

- `Halcyon6MV/` — Halcyon 6 MV model implementation
- `Infinity6MV/` — Infinity 6 MV model implementation
- `Infinity6MV_PSF/` — Variant for PSF computations

Each subfolder typically contains:

- `CMakeLists.txt` — build configuration
- `src/` — C++ source files (DetectorConstruction, PrimaryGenerator, PhysicsList, etc.)
- `include/` — header files
- `macros/` — example Geant4 macro files to run simulations
- `build/` — build artifacts and project files (not recommended to commit)
- `data.ipynb`, `merged_output.csv` — example notebooks and output data

## Requirements

- Geant4 (compatible version; Geant4 10.x or newer recommended)
- CMake
- A modern C++ compiler (MSVC, GCC or Clang)
- Python (optional, for notebooks and analysis)

## Build (Windows)

Open the **Developer Command Prompt for Visual Studio** and create a build directory inside the chosen subproject. Then run:

```powershell
cd Halcyon6MV    # or whichever subproject you want to build
mkdir build; cd build
cmake ..
cmake --build . --config Release
```

If you need to select a specific generator / Visual Studio version, run `cmake .. -G "Visual Studio 16 2019" -A x64` instead of `cmake ..`.

Note: on Linux or macOS use the system generator and standard terminal (for example `cmake .. && make`).

## Run

After building, run the simulation executable with the following command-line interface:

```
sim --t <num_threads> --s <seed> [macro_filename]
```

- `--t <num_threads>` : number of threads to use. If omitted, the program will use all available CPU cores.
- `--s <seed>`        : random seed. If omitted, a random seed will be chosen.
- `[macro_filename]`  : path to a Geant4 macro file to execute (for example `macros/default.mac`). If omitted, the default macro will be used.

Examples:

```powershell
.\sim.exe --t 8 --s 12345 macros\default.mac
.\sim.exe                 # runs with default macro, random seed, all threads
```

Notebooks (`data.ipynb`) are provided for post-processing and analysis of simulation outputs.

## Best practices

- Do not commit `build/` to the repository (add it to `.gitignore`).
- Keep detector geometry code (`DetectorConstruction`) separate from physics lists (`PhysicsList`) to facilitate testing.
- Always annotate units explicitly (`cm`, `mm`, `m`, etc.) when defining geometrical dimensions.

## Contributing

If you want to contribute:

- Open an issue describing the improvement or bug.
- Create a feature branch for your work.
- Send a pull request including tests or usage examples when appropriate.

## License

Add your preferred license here (MIT, Apache-2.0, etc.) or state the institutional policy.

---

I can also help with:

- Adding a recommended `.gitignore` for Geant4/CMake projects on Windows.
- Generating a basic `CONTRIBUTING.md`.
- Creating per-subfolder `README.md` files with tailored build/run instructions.

Tell me which of these you'd like and I'll add them.
