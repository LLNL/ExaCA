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
        // Variables with values taken from input file "MasterInputs.txt"
        int DecompositionStrategy;
        double deltax, AConst, BConst, CConst, NMax, dTN, dTsigma;
        string BaseFileName, GrainOrientationFile, TemperatureDataType;
        
        MasterInputRead(DecompositionStrategy, deltax, AConst, BConst, CConst, NMax, dTN, dTsigma, BaseFileName, GrainOrientationFile, TemperatureDataType);
        
        // Run CA code based on bulky of reduced version of data
        if (TemperatureDataType == "B") {
            RunProgram_Bulky(id,np,ierr,DecompositionStrategy, deltax, AConst, BConst, CConst, NMax, dTN, dTsigma, BaseFileName, GrainOrientationFile);
        }
        else {
            RunProgram_Reduced(id,np,ierr,DecompositionStrategy, deltax, AConst, BConst, CConst, NMax, dTN, dTsigma, BaseFileName, GrainOrientationFile, TemperatureDataType);
        }
        
    }
    // Finalize Kokkos
    Kokkos::finalize();
    // Finalize MPI
    ierr = MPI_Finalize ( );
    return 0;
}
