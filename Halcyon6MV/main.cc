#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <thread>

#include <random>
#include <filesystem>
#include "G4Types.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4AnalysisManager.hh"
#include "G4MTRunManager.hh"
#include "G4ScoringManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "GlobalSeed.hh"

int main(int argc, char** argv) {
    // ------------------ SEED y THREADS ------------------
        
        time_t current_time = time(0);
        
        // Usar el tiempo actual como semilla para el generador Mersenne Twister
        std::mt19937 mersenne_generator(current_time);
        
        // Definir un rango para los números aleatorios que queremos generar
        std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());  // Rango completo de un int
        
        // Generar el número aleatorio utilizando el generador Mersenne Twister
        seed = dist(mersenne_generator);
        
        // Mostrar la semilla generada
        std::cout << "Semilla generada: " << seed << std::endl;

        int numThreads = std::thread::hardware_concurrency();

        std::string macroFile;
        
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--s" && i + 1 < argc) {
                seed = std::stoi(argv[++i]);
            } else if (arg == "--t" && i + 1 < argc) {
                numThreads = std::stoi(argv[++i]);
            } else if (arg.rfind(".mac") != std::string::npos) {
                macroFile = arg;
            }
        }

    // ------------------ RANDOM SETUP ------------------
        G4Random::setTheEngine(new CLHEP::RanecuEngine);
        G4Random::setTheSeed(seed);

    // ------------------ LOG DE LA EJECUCIÓN ------------------
        std::filesystem::create_directory("run_logs");
        auto now = std::chrono::system_clock::now();
        std::time_t timestamp = std::chrono::system_clock::to_time_t(now);
        std::stringstream filename;
        filename << "run_logs/run_" << std::put_time(std::localtime(&timestamp), "%Y%m%d_%H%M%S")
                << "_seed_" << seed << "_threads_" << numThreads << ".txt";

        std::ofstream logfile(filename.str());
        logfile << "Seed: " << seed << std::endl;
        logfile << "Threads: " << numThreads << std::endl;
        logfile.close();

    // ------------------ RUN MANAGER ------------------
        G4MTRunManager *runManager = new G4MTRunManager();
        runManager->SetNumberOfThreads(numThreads);
        runManager->SetUserInitialization(new MyDetectorConstruction());
        runManager->SetUserInitialization(new MyPhysicsList());
        runManager->SetUserInitialization(new MyActionInitialization());
        runManager->Initialize();

    // ------------------ UI Y VISUALIZACIÓN ------------------
        G4UIExecutive *ui= new G4UIExecutive(argc,argv);
        G4VisManager *visManager = new G4VisExecutive();
        visManager->Initialize();
        G4UImanager *UImanager = G4UImanager::GetUIpointer();

    // ------------------ MACRO ------------------
        if (!macroFile.empty()) {
            UImanager->ApplyCommand("/control/execute " + macroFile);
        } else {
            std::string defaultMacro = "../../macros/default.mac";
            std::ifstream checkFile(defaultMacro);
            if (checkFile.good()) {
                G4cout << "Cargando macro por defecto: " << defaultMacro << G4endl;
                UImanager->ApplyCommand("/control/execute " + defaultMacro);
                

            } else {
                G4cerr << "No se encontró el macro. Ejecutando en modo interactivo." << G4endl;
                UImanager->ApplyCommand("/control/execute vis.mac");
                if (ui) {
                    ui->SessionStart();
                    delete ui;
                }
            }
        } 

        ui->SessionStart();
    return 0;
}
