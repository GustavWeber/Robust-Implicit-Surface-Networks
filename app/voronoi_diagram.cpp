//
// Created by Charles Du on 1/15/22.
//
#include <simplicial_arrangement/lookup_table.h>

#include <iostream>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

//#define WRITEACTIVEFUNCS

#include "voronoi_diagram.h"

#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;

using namespace simplicial_arrangement;


int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
        bool robust_test = false;
    } args;
    CLI::App app{"Material Interface Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_flag("-T,--timing-only", args.timing_only, "Record timing without saving results");
    app.add_flag("-R,--robust-test",args.robust_test, "Perform robustness test");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    Config config = parse_config_file(args.config_file);
    if (config.use_lookup) {
        // load lookup table
        std::cout << "load table ..." << std::endl;
        bool loaded = load_lookup_table(simplicial_arrangement::MATERIAL_INTERFACE);
        if (loaded) {
            std::cout << "loading finished." << std::endl;
        } else {
            std::cout << "loading failed." << std::endl;
            return -1;
        }
    } else {
        disable_lookup_table();
        config.use_secondary_lookup = false;
    }

    output_dir = config.output_dir;

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    if (config.tet_mesh_file != "") {
        std::cout << "load mesh file " << config.tet_mesh_file << std::endl;
        load_tet_mesh(config.tet_mesh_file, pts, tets);
    } else {
        std::cout << "generating mesh with resolution "
            << config.tet_mesh_resolution << std::endl;
        generate_tet_mesh(config.tet_mesh_resolution, config.tet_mesh_bbox_min,
                config.tet_mesh_bbox_max, pts, tets);
    }

    //Load Points
    std::vector<Point_3> points;
    if (load_points(config.func_file, points)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }

    // compute implicit arrangement
    std::vector<std::array<double, 3>> MI_pts;
    std::vector<PolygonFace> MI_faces;
    std::vector<std::vector<size_t>> patches;
    std::vector<std::pair<size_t, size_t>> patch_function_label;
    std::vector<Edge> MI_edges;
    std::vector<std::vector<size_t>> chains;
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> shells;
    std::vector<std::vector<size_t>> material_cells;
    std::vector<size_t> cell_function_label;
    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;
    // record stats
    std::vector<std::string> stats_labels;
    std::vector<size_t> stats;

    if (!voronoi_diagram<Kernel>(
            args.robust_test,
            config.use_lookup,
            config.use_secondary_lookup,
            config.use_topo_ray_shooting,
            config.use_original_filter,
            //
            pts, tets, points,
            //
            MI_pts,MI_faces,patches,
            patch_function_label,
            MI_edges,chains,
            non_manifold_edges_of_vert,
            shells,material_cells,cell_function_label,
            timing_labels,timings,
            stats_labels,stats)) {
        return -1;
    }
    if (args.robust_test) return 0;

    // test: export MI_mesh, patches, chains
    if (!args.timing_only && material_cells.size() > 0) {
        save_result_MI(config.output_dir + "/mesh.json",
                                         MI_pts,
                                         MI_faces,
                                         patches,
                                         patch_function_label,
                                         MI_edges,
                                         chains,
                                         non_manifold_edges_of_vert,
                                         shells,
                                         material_cells,
                                         cell_function_label);
        //
        if (material_cells.front().front() != Mesh_None){
            save_result_msh(config.output_dir + "/mesh",
                            MI_pts,
                            MI_faces,
                            patches,
                            MI_edges,
                            chains,
                            non_manifold_edges_of_vert,
                            shells,
                            material_cells);
        }
    }
    // save timing records
    save_timings(config.output_dir + "/timings.json", timing_labels, timings);
    // save statistics
    save_statistics(config.output_dir + "/stats.json", stats_labels, stats);

    return 0;
}
