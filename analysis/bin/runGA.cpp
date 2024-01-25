// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAorientation.hpp"
#include "CAparsefiles.hpp"
#include "GArepresentativeregion.hpp"
#include "GAutils.hpp"

#include "ExaCA.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// The path to/name of the ExaCA analysis input file and the path to/ name of the base filename for the ExaCA
// microstructure data (without extension) are given on the command line
int main(int argc, char *argv[]) {

    // Initialize Kokkos
    Kokkos::initialize();
    {
        using memory_space = typename Kokkos::DefaultExecutionSpace::memory_space;
        // Read command line input to obtain name of analysis file
        std::string analysis_file;
        if (argc < 2) {
            throw std::runtime_error("Error: Full path to and name of analysis file must be given on the command line");
        }
        analysis_file = argv[1];
        std::string base_filename = argv[2];
        std::string log_file = base_filename + ".json";
        std::string microstructure_file = base_filename + ".vtk";
        std::cout << "Performing analysis of " << microstructure_file << " , using the log file " << log_file
                  << " and the options specified in " << analysis_file << std::endl;

        std::string grain_unit_vector_file;
        double deltax;
        int nx, ny, nz, number_of_layers;
        std::vector<double> xyz_bounds(6);
        parseLogFile(log_file, nx, ny, nz, deltax, number_of_layers, xyz_bounds, grain_unit_vector_file, false);

        // Allocate memory blocks for grain_id and layer_id data
        Kokkos::View<int ***, Kokkos::HostSpace> grain_id(Kokkos::ViewAllocateWithoutInitializing("grain_id"), nz, nx,
                                                          ny);
        Kokkos::View<short ***, Kokkos::HostSpace> layer_id(Kokkos::ViewAllocateWithoutInitializing("layer_id"), nz, nx,
                                                            ny);

        // Fill arrays with data from paraview file
        initializeData(microstructure_file, nx, ny, nz, grain_id, layer_id);

        // Grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        Orientation<memory_space> orientation(grain_unit_vector_file, true);

        // Representative region creation
        std::ifstream analysis_data_stream(analysis_file);
        nlohmann::json analysis_data = nlohmann::json::parse(analysis_data_stream);
        nlohmann::json regions_data = analysis_data["Regions"];
        int number_of_regions = analysis_data["Regions"].size();
        std::cout << "There are " << number_of_regions << " regions to analyze" << std::endl;
        for (auto it = regions_data.begin(); it != regions_data.end(); it++) {
            // Create region
            std::string region_name = it.key();
            RepresentativeRegion representativeregion(analysis_data, region_name, nx, ny, nz, deltax, xyz_bounds,
                                                      grain_id);
            std::string base_filename_this_region = base_filename + "_" + region_name;

            // Output file stream for quantities of interest
            std::ofstream qois;
            std::string qois_fname = base_filename_this_region + "_QoIs.txt";
            qois.open(qois_fname);
            // Header data for qois file
            representativeregion.printAnalysisHeader(qois);

            // Fraction of region consisting of nucleated grains, unmelted material
            if (representativeregion.analysis_options_stats_yn[0])
                representativeregion.printGrainTypeFractions(qois, grain_id, layer_id);

            // Calculate and if specified, print misorientation data
            std::vector<float> grain_misorientation_x_vector =
                representativeregion.getGrainMisorientation("X", orientation);
            std::vector<float> grain_misorientation_y_vector =
                representativeregion.getGrainMisorientation("Y", orientation);
            std::vector<float> grain_misorientation_z_vector =
                representativeregion.getGrainMisorientation("Z", orientation);
            if (representativeregion.analysis_options_stats_yn[1])
                representativeregion.printMeanMisorientations(
                    qois, grain_misorientation_x_vector, grain_misorientation_y_vector, grain_misorientation_z_vector);

            // Print mean size data if specified
            if (representativeregion.analysis_options_stats_yn[2])
                representativeregion.printMeanSize(qois);

            // If XExtent, YExtent, ZExtent, or build_trans_aspect_ratio/Extent are toggled for general stats printing
            // or per grain printing, calculate grain extents for the necessary direction(s) (otherwise don't, since
            // it can be slow for large volumes) Extents are calculated in microns
            // If options StatsYN[3] or PerGrainStatsYN[5] is toggled, grain extents are needed
            // If StatsYN[4] or PerGrainStatsYN[2] is toggled, X extents are needed
            // If StatsYN[5] or PerGrainStatsYN[3] is toggled, Y extents are needed
            // If StatsYN[6] or PerGrainStatsYN[4] is toggled, Z extents are needed
            representativeregion.calcNecessaryGrainExtents(grain_id, deltax);
            std::vector<float> build_trans_aspect_ratio(representativeregion.number_of_grains);
            if ((representativeregion.analysis_options_stats_yn[3]) ||
                (representativeregion.analysis_options_per_grain_stats_yn[5]))
                representativeregion.calcBuildTransAspectRatio(build_trans_aspect_ratio);
            if (representativeregion.analysis_options_stats_yn[3])
                representativeregion.printMeanBuildTransAspectRatio(qois);

            if (representativeregion.analysis_options_stats_yn[4])
                representativeregion.printMeanExtent(qois, "X");
            if (representativeregion.analysis_options_stats_yn[5])
                representativeregion.printMeanExtent(qois, "Y");
            if (representativeregion.analysis_options_stats_yn[6])
                representativeregion.printMeanExtent(qois, "Z");

            // Determine IPF-Z color of each grain relative to each direction: 0 (red), 1 (green), 2 (blue)
            std::vector<float> grain_red = representativeregion.getIPFZColor(0, orientation);
            std::vector<float> grain_green = representativeregion.getIPFZColor(1, orientation);
            std::vector<float> grain_blue = representativeregion.getIPFZColor(2, orientation);

            // Write grain area data as a function of Z location in the representative volume if the options are
            // toggled, writing to files
            // "[base_filename_this_region]_WeightedGrainAreas.csv" and "[base_filename_this_region]_GrainAreas.csv",
            // respectively
            if ((representativeregion.analysis_options_layerwise_stats_yn[0]) ||
                (representativeregion.analysis_options_layerwise_stats_yn[1]))
                representativeregion.writeAreaSeries(base_filename_this_region, deltax, grain_id);
            qois.close();

            // Write per-grain stats for the analysis types specified to the file
            // "[base_filename_this_region]_grains.csv"
            if (representativeregion.print_per_grain_stats_yn)
                representativeregion.writePerGrainStats(base_filename_this_region, grain_misorientation_x_vector,
                                                        grain_misorientation_y_vector, grain_misorientation_z_vector,
                                                        build_trans_aspect_ratio, grain_red, grain_green, grain_blue);

            // ExaConstit print a file named "[base_filename_this_region]_ExaConstit.csv"
            if (representativeregion.print_exaconstit_yn) {
                representativeregion.writeExaConstitRVE(base_filename_this_region, deltax, grain_id);
            }

            // Pole figure print a file named "[base_filename_this_region]_PoleFigureData.txt"
            if (representativeregion.print_pole_figure_yn) {
                auto go_histogram =
                    representativeregion.getOrientationHistogram(orientation.n_grain_orientations, grain_id, layer_id);
                representativeregion.writePoleFigure(base_filename_this_region, orientation, go_histogram);
            }

            // IPF map for area print a file named "[base_filename_this_region]_IPFCrossSectionData.txt"
            if (representativeregion.print_inverse_pole_figure_map_yn) {
                representativeregion.writeIPFColoredCrossSection(base_filename_this_region, grain_id, orientation,
                                                                 deltax);
            }
            std::cout << "Finished analysis for region " << region_name << std::endl;
        } // end loop over all representative regions in analysis file
    }     // end scope for kokkos
    // Finalize kokkos and end program
    Kokkos::finalize();
    return 0;
}
