#ifndef FORCE_WRITER_H_
#define FORCE_WRITER_H_

#include <fstream>
#include <iostream>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

class ForceWriter {
public:
    // Constructor: filename to write to
    ForceWriter(const std::string& filename)
        : filename_(filename) {}

    // Write a single double value to file (rank 0 only if MPI)
    void write(double value) const {
        int rank = 0;

        #ifdef USE_MPI
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #endif
        
        if (rank == 0) {
            std::ofstream outfile(filename_, std::ios::app);
            if (!outfile) {
                std::cerr << "Error opening file: " << filename_ << "\n";
                return;
            }
            outfile << value << "\n";  // append value to file
        }
    }

private:
    std::string filename_;
};

#endif // FORCE_WRITER_H_