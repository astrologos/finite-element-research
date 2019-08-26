/* ---------------------------------------------------------------------
 *
 * This is a test file adapted from the deal.ii library examples
 * by Jack Alvarez.
 * 
 * This file is used for performance testing
 * involving parallelism and mesh refinement.
 *
 * ---------------------------------------------------------------------
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Jack Alvarez
 * Institution: Virginia Polytechnic Institute and State University
 * Date: 3/19/2-019
 *
 */

/*
 * NOTE: Using trilinos SSOR preconditioner
 */
#include <string>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_solver.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <array>


// The last step is as in all previous programs:
namespace current
{
  using namespace dealii;

  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem (const unsigned int poly_degree);
    void run();

  private:
    void setup_system();
    void assemble_and_solve();
    void solve();

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    MappingQ<dim>        mapping;

    SparsityPattern      sparsity;
    SparseMatrix<double> A;
    ConstraintMatrix     mvcs;

    Vector<double>       x;
    Vector<double>       b;

    TableHandler         table_out;
    TrilinosWrappers::SparseMatrix trilinosA;
  };



  // Initialize LaplaceProblem
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int poly_degree):
    fe(1),
    dof_handler(triangulation),
    mapping(poly_degree)
  {
    std::cout << "Using mapping with degree " << poly_degree << ":"
              << std::endl
              << "============================"
              << std::endl;
  }


  // Set up linear system
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    // Distribute DOFs then adjust system accordingly
    dof_handler.distribute_dofs(fe);
    x.reinit(dof_handler.n_dofs());
    b.reinit(dof_handler.n_dofs());

    // List nodes at the boundary
    std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
    DoFTools::extract_boundary_dofs(dof_handler,
                                     ComponentMask(),
                                     boundary_dofs);

    // Select the first boundary node
    const unsigned int first_boundary_dof
      = std::distance(boundary_dofs.begin(),
                       std::find(boundary_dofs.begin(),
                                  boundary_dofs.end(),
                                  true));

    // Clear previous constraints
    // Then add line constraining the first boundary node DOF to the sum
    // of other boundary node DOFs with weight -1.
    mvcs.clear();
    mvcs.add_line(first_boundary_dof);
    for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true)
        mvcs.add_entry(first_boundary_dof,
                                          i, -1);
    mvcs.close();

    // Use dynamic sparsity patterns
    // Use Cuthill-McKee method to reduce bandwidth
        DoFRenumbering::Cuthill_McKee (dof_handler);
    DynamicSparsityPattern dsp(dof_handler.n_dofs(),
                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    mvcs.condense(dsp);

    // Copy new sparsity pattern from dsp
    // Then reinitialize matrix with boundary value constraints
    sparsity.copy_from(dsp);
    A.reinit(sparsity);
    std::ofstream out ("sparsity_pattern1.svg");
    sparsity.print_svg(out);
    
  }



  // Solve the Laplace problem
  template <int dim>
  void LaplaceProblem<dim>::assemble_and_solve()
  {
    // out << Number of active cells
    table_out.add_value("cells,", std::to_string(triangulation.n_active_cells())+",");
    // Get gaussian quadrature degree
    const unsigned int gauss_degree
      = std::max(static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)),
                  2U);
    // Use MatrixTools to create matrix and RHS
    MatrixTools::create_laplace_matrix(mapping, dof_handler,
                                        QGauss<dim>(gauss_degree),
                                        A);
    // f=(-2)
    VectorTools::create_right_hand_side(mapping, dof_handler,
                                         QGauss<dim>(gauss_degree),
                                         ConstantFunction<dim>(-2),
                                         b);
    // g=(1)
    Vector<double> tmp(b.size());
    VectorTools::create_boundary_right_hand_side(mapping, dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(1),
                                                  tmp);
    // Now add boundary contributions
    b += tmp;
    
    // Eliminate constrained DOFs
    mvcs.condense(A);
    mvcs.condense(b);
    // Solve
    solve();
    
    // Then redistribute constraints to the rest of the DOFs
    mvcs.distribute(x);

    // Find norm of solution
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping, dof_handler,
                                       x,
                                       ZeroFunction<dim>(),
                                       norm_per_cell,
                                       QGauss<dim>(gauss_degree+1),
                                       VectorTools::H1_seminorm);

    const double norm = VectorTools::compute_global_error(triangulation,
                                                          norm_per_cell,
                                                          VectorTools::H1_seminorm);

    // ----------------- Generate Output -----------------
    table_out.add_value("|u|_1,", std::to_string(norm)+",");
    table_out.add_value("error;", std::to_string(std::fabs(norm-std::sqrt(3.14159265358/2)))+";");
  }



  // Solve the system
  template <int dim>
  void LaplaceProblem<dim>::solve()
  {
    SolverControl                       solver_control(100000, 1e-10);
    SolverCG<>			        cg(solver_control);
    trilinosA.reinit(A);

    std:: cout << ".." << std::flush;
    TrilinosWrappers::PreconditionSSOR preconditioner;
    Timer timer0;
    timer0.start();
    preconditioner.initialize(trilinosA, 1.2);
    timer0.stop();
    std:: cout << ".." << std::flush;
    table_out.add_value("elapsed CPU time in Preconditioning (sec),",std::to_string(timer0())+",");      
    table_out.add_value("elapsed Wall time in Preconditioning (sec),",std::to_string(timer0.wall_time())+","); 

    // ----------------- Time -----------------
    Timer timer1;
    timer1.start();
    cg.solve(trilinosA, x, b, preconditioner);
    timer1.stop();
    table_out.add_value("elapsed CPU time in CG (sec),",std::to_string(timer1())+",");      
    table_out.add_value("elapsed Wall time in CG (sec),",std::to_string(timer1.wall_time())+","); 
    std:: cout << ".." << std::flush;
  }



  // Run the Laplace problem for increasingly finer mesh
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    GridGenerator::hyper_ball(triangulation);
    static const SphericalManifold<dim> boundary;
    triangulation.set_all_manifold_ids_on_boundary(0);
    triangulation.set_manifold(0, boundary);
    triangulation.refine_global(4);
    std:: cout << "<------------------------------>\n";
    std::cout << "<";
    for (unsigned int cycle=0; cycle<5; ++cycle, triangulation.refine_global(1))
      {
        setup_system();
        assemble_and_solve();
      };
      std::cout << ">\n\n";

    // Write table out to console:
    table_out.set_precision("|u|_1,", 6);
    table_out.set_precision("error;", 6);
    table_out.set_precision("elapsed CPU time in Preconditioning (sec),",3);
    table_out.set_precision("elapsed Wall time in Preconditioning (sec),",3);
    table_out.set_precision("elapsed CPU time in CG (sec),",3);
    table_out.set_precision("elapsed Wall time in CG (sec),",3);
    table_out.write_text(std::cout);
    std::cout << std::endl;
  }
}



// Execute
int main(int argc, char **argv)
{
    using namespace dealii;
    using namespace current;
    
    Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv,
                                                           numbers::invalid_unsigned_int);
    AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)==1,
        ExcMessage("This program can only be run in serial"));
    // RUN ON MACHINE WITH 2 8-CORE CPUS WITH SHARED MEMORY, 2 THREADS PER CORE,
    // I.E. 32 THREAD CAPACITY
    for (int threads = 1; threads<33; threads*=2) {

        dealii::MultithreadInfo::set_thread_limit(threads);
        std:: cout << "\nStarting run with " << threads <<  " threads::\n";    

        try
        {
          std::cout.precision(5);

          // Run problem for mapping up to degree 3
          for (unsigned int poly_degree=1; poly_degree<=3; ++poly_degree) {
          
            current::LaplaceProblem<2>(poly_degree).run();
            }
        }
        catch (std::exception &exc)
        {
          std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "Exception on processing: " << std::endl
                    << exc.what() << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          return 1;
        }
        catch (...)
        {
          std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "Unknown exception!" << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          return 1;
            };
    }
    return 0;
}

