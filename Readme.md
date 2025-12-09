# G4Linacs

Simulaciones de aceleradores lineales (Linacs) y cabezales clínicos en Geant4.

## Descripción

G4Linacs es una colección de escenarios y configuraciones para simular fuentes de radiación y cabezales de aceleradores lineales médicos usando Geant4. El objetivo del proyecto es proporcionar modelos reutilizables (ej. `Halcyon6MV`, `Infinity6MV`, `Infinity6MV_PSF`) que permitan estudiar la distribución de dosis, perfilar la respuesta de detectores y generar funciones de dispersión puntual (PSF) y otros parámetros relevantes para física médica.

Cada subcarpeta contiene una variante del modelo (código fuente, archivos de construcción CMake, macros y datos de salida) pensada para diferentes linacs o configuraciones de simulación.

## Estructura del repositorio

- `Halcyon6MV/` — Implementación del modelo Halcyon a 6 MV
- `Infinity6MV/` — Implementación del modelo Infinity a 6 MV
- `Infinity6MV_PSF/` — Variante para cálculo de PSF
- Cada carpeta contiene:
  - `CMakeLists.txt` — archivo de construcción
  - `src/` — código fuente C++ (DetectorConstruction, PrimaryGenerator, PhysicsList, etc.)
  - `include/` — cabeceras
  - `macros/` — macros de ejemplo para correr simulaciones
  - `build/` — archivos de proyecto y salidas de compilación (no versionar si se usa .gitignore)
  - `data.ipynb`, `merged_output.csv` — notebooks y datos de ejemplo

## Requisitos

- Geant4 (versión compatible; recomendado 10.x o superior)
- CMake
- Compilador C++ moderno (MSVC, GCC o Clang)
- Python (opcional, para notebooks y análisis)

## Cómo compilar

Ejemplo rápido (desde Windows PowerShell):

```powershell
mkdir build; cd build
cmake .. -G "Visual Studio 16 2019" -A x64
cmake --build . --config Release
```

Ajústate al generador y la configuración de Visual Studio que tengas instalada. En Linux/macOS usa el generador y compilador correspondiente.

## Cómo ejecutar

Una vez compilado, puedes ejecutar el binario principal (por ejemplo `sim.exe` o `sim`) y pasarle una macro de Geant4:

```powershell
.\sim.exe macros\default.mac
```

También se incluyen notebooks (`data.ipynb`) para análisis de resultados.

## Buenas prácticas

- No versionar la carpeta `build/` (añade una entrada en `.gitignore`).
- Mantener separados el código del detector (`DetectorConstruction`) y la lista de física (`PhysicsList`) para facilitar pruebas.
- Anotar unidades explícitamente (`cm`, `mm`, `m`, etc.) al definir dimensiones.

## Contribuciones

Si deseas contribuir:

- Abre un issue describiendo la mejora o bug.
- Crea una rama para tu trabajo.
- Envía un pull request con pruebas o ejemplos de uso.

## Licencia

Añade aquí la licencia de tu preferencia (MIT, Apache-2.0, etc.) o déjalo bajo la política del grupo/institución.

---

Si quieres, puedo:

- Añadir un `.gitignore` recomendado para proyectos Geant4/CMake en Windows.
- Generar un `CONTRIBUTING.md` básico.
- Formatear los `README.md` en cada subcarpeta con instrucciones específicas.

Dime qué prefieres y lo hago.
