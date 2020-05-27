#include "header.h"
using namespace std;

int main ( int argc, char *argv[] ) {
    // Initialize MPI
    int ierr, id, np;
    ierr = MPI_Init ( &argc, &argv );
    // Initialize Kokkos
    Kokkos::initialize(argc, argv); {
        
        // Get number of processes
        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &np );
        // Get individual process ID
        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
        
        if (id == 0) Kokkos::DefaultExecutionSpace::print_configuration(std::cout);
        if (id == 0) cout << "Number of MPI ranks = " << np << endl;
        
        char* InputFileEnv = getenv ("CAINPUT");
        if (InputFileEnv == NULL) {
            if (id == 0) cout << "Error: Input file not found, program will not run" << endl;
        }
        else {
            // Run CA code using reduced temperature data format
            string InputFile(InputFileEnv);
            RunProgram_Reduced(id, np, ierr, InputFile);
        }
        
    }
    // Finalize Kokkos
    Kokkos::finalize();
    // Finalize MPI
    ierr = MPI_Finalize ( );
    return 0;
}
