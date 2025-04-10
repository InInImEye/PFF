/* ---------------------------------------------------------------------
 *
 * 
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Sparsh Sinha
 * Program working discontinuous staggered step */

// Base Headers
#include <deal.II/base/quadrature_lib.h> // Quadrature
#include <deal.II/base/function.h> // Funtion class
#include <deal.II/base/logstream.h> // Logging dealii
#include <deal.II/base/tensor.h> // Tensor
#include <deal.II/base/symmetric_tensor.h> // Implements symmetric tensors of rank 2 and 4

// Vector and Matrix libraries
#include <deal.II/lac/vector.h> // Vector
#include <deal.II/lac/full_matrix.h> // Full Dense matrix, whether 0 or non 0
#include <deal.II/lac/sparse_matrix.h> // Sparse Matrix, only stored non 0
#include <deal.II/lac/dynamic_sparsity_pattern.h> // Dynamic matrix, grows, at end needs to be converted to Sparse Matrix
#include <deal.II/lac/solver_cg.h> // Conjugate Gradient solver
#include <deal.II/lac/precondition.h> // Preconditions, simplifies linear systems for solving
#include <deal.II/lac/affine_constraints.h> // Applies constraints on dofs boundary or refinement

// Parallel headers
#include <deal.II/base/multithread_info.h> // Useful in multithread programming
#include <deal.II/lac/petsc_vector.h> // Parallelize vector and operations
#include <deal.II/lac/petsc_sparse_matrix.h> // Partitions and parallelize matrices
#include <deal.II/lac/petsc_solver.h> // Parallelize solver
#include <deal.II/lac/petsc_precondition.h> // Again parallel precondition
#include <deal.II/lac/sparsity_tools.h> // Performs function on sparsity patterns, renumbering, etc

// Triangulation objects, generators, tools, iterators
#include <deal.II/grid/tria.h> // Triangulation object
#include <deal.II/grid/grid_generator.h> // Grid Generator
#include <deal.II/grid/grid_refinement.h> // Grid refinement and functions
#include <deal.II/grid/manifold_lib.h> // Manifold lib, Polar, Spherical, Flat, etc.
#include <deal.II/grid/grid_tools.h> // Grid manipulation tools
#include <deal.II/grid/tria_accessor.h> // Triagulation object accessor, e.g. vertices, child, etc
#include <deal.II/grid/tria_iterator.h> // Iterator to face, cells, vertices, etc.
// And a header that implements filters for iterators looping over all
// cells. We will use this when selecting only those cells for output that are
// owned by the present process in a parallel program:
#include <deal.II/grid/filtered_iterator.h> // Filters interator using predicate, etc.

// DOF and tolls
#include <deal.II/dofs/dof_handler.h> // For a triangulation enumerates, and some functions
#include <deal.II/dofs/dof_accessor.h> // Accesses data for edges, faces and cells
#include <deal.II/dofs/dof_tools.h> // DOF manipulations tools, e.g. renumbering
#include <deal.II/dofs/dof_renumbering.h> // Renumbering algorithm

// FE and Q values and systems
#include <deal.II/fe/fe_values.h> // To assemble matrices, vectors. Link FE, Q, maps
#include <deal.II/fe/fe_system.h> // Vector values elements
#include <deal.II/fe/fe_q.h> // Implements scalar lagrange FE Q points
#include <deal.II/fe/fe_dgq.h> // Same thing but for discontinuous FE
#include <deal.II/fe/mapping_q.h> // Polynomial mappings Qp to get higher order mapping

// Tools, error and output
#include <deal.II/numerics/vector_tools.h> // Operations on vectors
#include <deal.II/numerics/matrix_tools.h> // Tools to apply on matrices, including BC
#include <deal.II/numerics/data_out.h> // Data out
#include <deal.II/numerics/error_estimator.h> // Error estimation per cell

// Solvers and others
#include <deal.II/lac/identity_matrix.h> // Identity matrix
#include <deal.II/lac/solver_minres.h> // Minimal residual method for symmetric matrices
#include <deal.II/lac/solver_gmres.h> // Generalized Minimal Residual Method, norm stops
#include <deal.II/lac/sparse_direct.h> // Sparse direct solver

// Parallel stuff
#include <deal.II/distributed/shared_tria.h> // Parallel triangulation every processor knows about every cell of the global mesh
#include <deal.II/base/conditional_ostream.h> // Conditional ostream pcout on mpi id
#include <deal.II/base/utilities.h> // Mainly MPI abstraction

// To compute rotation matrices of the local coordinate systems at specific
// points in the domain.
#include <deal.II/physics/transformations.h> // Assist in transformation of tensor quantities

// C++ headers
#include <fstream> // read/ write to file
#include <cmath> // cmath pow etc
#include <iostream> // cin cout cerr
#include <typeinfo> // Find type
#include <stdio.h> // printf scanf
#include <sstream> // string stream class
#include <unistd.h> // standard symbolic constants and types
#include <fcntl.h> // file control options
#include <cstdio> // C stdio
#include <iomanip> // manipulation i/o
#include "stdlib.h" // malloc free
#include <chrono> // time
#include <random> // random no.

// Removable headers
#include "utils.cc"
// Debug and other options
#include "inputFiles/options.h"

// Namespace Problem
namespace Problem
{
  using namespace dealii;
  
  // Parses and removes spaces from string
  std::string removeSpaces(std::string input)
  {
    input.erase(std::remove(input.begin(),input.end(),' '),input.end());
    return input;
  }
  
  // Point history class
  template <int spacedim>
  class PointHistory
  {   
    public:
    PointHistory();
    ~PointHistory();
    
    ////// Main State variables
    Tensor<2, spacedim> eps_p, eps_p_NR;//eps_e_NR for next NR iter, used for storing update quadrature history till energy convergence
    Tensor<2, spacedim> strss, strss_NR;
    
    double alpha, alpha_NR;
    
    // $$ p := \frac{\varepsilon^{p}_{eq}}{\varepsilon^{p}_{eq,crit}} $$
    double acc_eps_p, acc_eps_p_NR; //Accumulated plastic strain variable
    double Psi_plus_max, // $$ \mathcal{H}_{e} $$
           Psi_plus_max_NR; // Psi_plus_max_NR is changed in each staggered step
    
    ///// Mainly used for storage and output of other variables
    Tensor<2, spacedim> eps_e, eps_e_NR; //when solving for next timestep, eps_e means at current timestep
    
    double del_lambda;
    double g_hat; // $$ g(s, p) := s^{2p} + \eta $$
    
    //unsigned int Qindex_u; For check the index order
    //unsigned int Qindex_s;
  };

  template <int spacedim> PointHistory<spacedim>::PointHistory()
  {// A little shortcut but may need to be initialized in *_quadrature_point_history() fns
    
    for(unsigned int i=0;i<spacedim;i++)
      for(unsigned int j=0;j<spacedim;j++)
      {
        eps_e_NR[i][j]= 0.0; 
        eps_p_NR[i][j]= 0.0;
        strss_NR[i][j]= 0.0;
      }
    eps_e = eps_e_NR;
    eps_p = eps_p_NR;
    strss = strss_NR;
    
    alpha_NR = 0.0;
    alpha = alpha_NR; 
    
    acc_eps_p_NR = 0.0;
    acc_eps_p = acc_eps_p_NR;
    
    Psi_plus_max_NR = 0.;
    Psi_plus_max = Psi_plus_max_NR;
    
    del_lambda = 0;
    g_hat = 0;
    
    //Qindex_u = 0;
    //Qindex_s = 0;
    
  }

  template <int spacedim>
  PointHistory<spacedim>::~PointHistory()
  {}
  

  template <int dim>
  inline Tensor<2, dim> get_strain(const FEValues<dim> &fe_values,
                                            const unsigned int   shape_func,
                                            const unsigned int   q_point)
  {
    Tensor<2, dim> tmp;

    //for (unsigned int i = 0; i < dim; ++i)
    //  tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        tmp[i][j] =
          (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
           fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
          2;

    return tmp;
  }


  template <int dim>
  inline Tensor<2, dim>
  get_strain(const std::vector<Tensor<1, dim>> &grad)
  {
    Assert(grad.size() == dim, ExcInternalError());
    //std::cout<<typeid(grad).name()<<std::endl;
    Tensor<2, dim> strain;
    //for (unsigned int i = 0; i < dim; ++i)
    //  strain[i][i] = grad[i][i];

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
    //std::cout<<"grad"<< grad <<std::endl;
    //std::cout<<"strain"<< strain <<std::endl;
    return strain;
  }
  
  template <int dim, int spacedim>
  inline Tensor<2, spacedim> // has some redundancy to prevent rewriting the code
  tensordimtospacedim(const Tensor<2, dim> tendim)
  {
    Tensor<2, spacedim> tenspacedim;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        tenspacedim[i][j] = tendim[i][j];
        
    return tenspacedim;
  }
  
  
  // Double contract right C:e
  template <int spacedim>
  Tensor<2,spacedim> double_contract(const Tensor<4, spacedim> &Cijkl,
                                     const Tensor<2, spacedim> &eps)
  {
  	Tensor<2,spacedim> double_contract_tensor;
  	  	  	
  	for(unsigned int i=0; i<spacedim; i++)
      for(unsigned int j=0; j<spacedim; j++)
        for(unsigned int k=0; k<spacedim; k++)
          for(unsigned int l=0; l<spacedim; l++)
            double_contract_tensor[i][j] += Cijkl[i][j][k][l] * eps[k][l];
  	
    return double_contract_tensor;
  }
  
  // C1:C2
  template <int spacedim>
  Tensor<4,spacedim> double_contract(const Tensor<4, spacedim> &Cijkl1,
                                     const Tensor<4, spacedim> &Cijkl2)
  {
  	Tensor<4,spacedim> double_contract_tensor;
  	  	  	
  	for(unsigned int i=0; i<spacedim; i++)
      for(unsigned int j=0; j<spacedim; j++)
        for(unsigned int k=0; k<spacedim; k++)
          for(unsigned int l=0; l<spacedim; l++)
            for(unsigned int m=0; m<spacedim; m++)
              for(unsigned int n=0; n<spacedim; n++)
                double_contract_tensor[i][j][k][l] += Cijkl1[i][j][m][n] * Cijkl2[m][n][k][l];
  	
    return double_contract_tensor;
  }
  
  // Dyadic product: C_ijkl = sig_ij * sig_kl
  template <int spacedim>
  Tensor<4,spacedim> dyadic_product(const Tensor<2, spacedim> &eps1,
                                    const Tensor<2, spacedim> &eps2)
  {
  	Tensor<4,spacedim> dyadic_product_tensor;
  	  	  	
  	for(unsigned int i=0; i<spacedim; i++)
      for(unsigned int j=0; j<spacedim; j++)
        for(unsigned int k=0; k<spacedim; k++)
          for(unsigned int l=0; l<spacedim; l++)
            dyadic_product_tensor[i][j][k][l] = eps1[i][j] * eps2[k][l];
  	
    return dyadic_product_tensor;
  }
  
// Class Phase Field Crack
  template <int dim, int spacedim = dim>
  class PFC
  {
  public:
    PFC();
    ~PFC();
    void run();

  private:
    void calc_bounding_box();
    void make_grid();

    void setup_system(const int); //To setup the entire system or to reset only the system matrices, system rhs and correction vectors

    void assemble_system_disp( const int ); //Assemble system for displacement with Dirichlet BC
    void assemble_system_dam();//To assemble system for damage with no Dirichlet BC
    
    void solve_NR_disp(); 
    void solve_NR_dam();

    unsigned int solve_linear_problem_disp();
    unsigned int solve_linear_problem_dam();
    
    void output_results(const unsigned int); // Output results and averages of various variables
    void output_Qp_results(const unsigned int); // Output fields of state variables
    void output_stress_strain_plot(); // Output stress strain plot for the geometery
/*
    
    void refine_initial_grid();
*/
    void setup_Qp_history(); // Qp: = Quadrature point //Initiallizes based on setup_system

    void update_Qp_history_disp_NR(); //only disp_NR state variables updated because after NR convergence only disp convergence ensured but not energy convergence
    void update_Qp_history_dam_NR(); // similar for damage
    
    // Accept
    void accept_NR_solution(); //Accept NR solutions (disp_NR, dam_NR) to final solutions (disp, dam) after Energy convergence
    void accept_NR_state(); // Accept state variables: stress, plastic strain ... after convergence
    void reset_vark();
    
    // Functions for computation purpose
    Tensor<4, spacedim> calc_strss_strain_tensor_constant(const Tensor<2, spacedim> &,
                                                    const double &);
    Tensor<2, spacedim> calc_strss_from_strain(const Tensor <4, spacedim> &, const Tensor<2, spacedim> &);
    
    //$$ \Psi^{+}_{e}(\pmb{\varepsilon}_{e}) $$
    double Psi_plus(const Tensor<2, spacedim> &);
    
    //$$ \Psi^{-}_{e}(\pmb{\varepsilon}_{e}) $$
    double Psi_minus(const Tensor<2, spacedim> &);
    
    //$$ \Psi_{tot} = \hat{g}* \Psi^{+}_{e} +  \Psi^{-}_{e} + \Psi_{p} + \Psi_{s} $$
    double Psi_total(const Tensor<2, spacedim> &,
                      const Tensor<1,dim> &,
                      const double &,
                      const double &,
                      const double &);
    
    
    double yield_function(const Tensor<2, spacedim>&, const double &);
    double delta_lambda (const double &, const double &);
    
    // $$ E_{\ell}(\pmb{\varepsilon}^{e}, \pmb{\varepsilon}^{p}, \alpha, s)  := 
     //   \int_{\Omega} \Psi_{tot} d \pmb{x} $$
    double compute_Energy();
    void calc_timestep(); // Calculates dt from default_dt and eps_p
    void predictor_corrector (const Tensor<2, spacedim> &,
                              const Tensor<2, spacedim> &,
                              const Tensor<2, spacedim> &,
                              const double &,
                              const double &,
                              double &,
                              double &,
                              double &,
                              Tensor<2, spacedim> &,
                              Tensor<2, spacedim> &,
                              Tensor<2, spacedim> &,
                              Tensor<4, spacedim> &);
    
    
    // Class Vars
    parallel::shared::Triangulation<dim> triangulation;
    std::string mesh_type;
    std::string mesh_filename;
    DoFHandler<dim> dof_handler_dam, dof_handler_disp; //Separate dofhandler for damage and displacement solve // Enums all faces, edges, cells, vertices, provides bases for Vh
    FE_Q<dim> fe_dam;
    FESystem<dim> fe_disp;

    //const FEValuesExtractors::Vector displacement(0);
    //const FEValuesExtractors::Scalar damage(dim);

    AffineConstraints<double> hanging_node_constraints_dam, hanging_node_constraints_disp;

    const QGauss<dim> quadrature_formula;//quadrature_formula_u; //quadrature_formula_s
    
    // SparsityPattern sparsity_pattern;
    PETScWrappers::MPI::SparseMatrix system_matrix_dam;    //System matrix for damage 
    PETScWrappers::MPI::SparseMatrix system_matrix_disp;  //System matrix for displacement 
    
    PETScWrappers::MPI::Vector system_rhs_dam; // system rhs/residue for damage
    PETScWrappers::MPI::Vector system_rhs_disp;// system rhs/residue for displacement
    PETScWrappers::MPI::Vector system_rxn;// system vector for calculating the reaction force
    
    Vector<double> dam; // Final damage solution  // $$ s $$
    Vector<double> dam_NR; //NR step solutions
    Vector<double> d_dam;// damage correction at each NR step
    
    Vector<double> disp; //Final displacement solution // $$ \pmb{u} $$
    Vector<double> disp_NR;  //NR step solutions       //u_NR
    Vector<double> d_disp; // displacement correction at each NR step
    
    // Displacement conditions
    unsigned int load_type;
    double       present_strain;// current displacemnt
    double       strain_rate; // velocity = strain_rate * Lenght y
    double       strain_min;
    double       strain_max;
    bool         stored_time;
    bool         enable_linear;
    double       velocity; // load step = velocity * dt
    double       time;
    double       time_sin_start;
    double       end_time;
    double       dt; // adjusted timestep scaled with plastic strain, calculated from default_dt
    double       default_dt; // Default timestep provided by the user
    unsigned int timestep_no;
    
    
    ////////// Constants
    double lambda; //Lame's constant // $$ \lambda $$
    double MU; //Shear modulus // $$ \mu $$
    double Kn; //Bulk modulus // $$ K_{n} $$
    double Eta; // Constant to prevent numerical error due to damage =0 // $$ \eta $$
    double Ystrss;// yield stress // $$ \sigma_{y} $$
    double H_pl; //hardening const // $$ h $$
    // $$ \varepsilon^{p}_{eq,crit} $$
    double Ep_crit; // Critical von-Mises equivalent plastic strain 
    double Gc; //Material Fracture Toughness // $$ G_{c} $$
    double L_dam; //Transition length for damage // $$ \ell $$
    // Checking var
    
    ////////// Control parameters
    bool enable_bounding_box; //Only usefull if model oriented in x y z directions
    unsigned int refinement_level;
    unsigned int outputFrequency;
    
    unsigned int N_stagger; //Max Staggered step
    unsigned int N_NR_disp; //Max NR step for displacement solve
    unsigned int N_NR_dam; //Max NR step for damage solve
    
    double Energy_tol; //Energy error tolerance
    double rhs_disp_tol; //displacemnt system rhs/residue relative tolerance
    double rhs_dam_tol; //damage system rhs/residue relative tolerance
    
    double lin_solve_disp_tol;
    double lin_solve_dam_tol;
    
    double d_eps_p_limit;
    double d_dam_max_limit;
    int d_eps_p_count = 0;
    
    MPI_Comm mpi_communicator;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;
    ConditionalOStream pcout;
    
    // std::vector<types::global_dof_index> local_dofs_per_process;
    
    IndexSet locally_owned_dofs_dam, locally_owned_dofs_disp;
    IndexSet locally_relevant_dofs_dam, locally_relevant_dofs_disp;
    
    // unsigned int n_local_cells;
        
    std::ofstream ssFile;

    std::ofstream outfile;
    std::ofstream matrixfile, residual_vectorfile, data_vectorfile;
    
    //#ifdef PFC_DEBUG_STAGGER
    //  std::ofstream debug_stagger;
    //  std::ofstream debug_stagger_list;
    //#endif
    //#ifdef PFC_DEBUG_NR_U
    //  std::ofstream debug_nr_disp;
    //  std::ofstream debug_nr_disp_list;
    //#endif
    //#ifdef PFC_DEBUG_NR_S
    //  std::ofstream debug_nr_dam;
    //  std::ofstream debug_nr_dam_list;
    //#endif
    //std::ofstream outfile, plastic_strain; 
    

    std::vector<PointHistory<spacedim>> quadrature_point_history;

    double size_x, size_y, size_z; 
    double n_subdivisions_x, n_subdivisions_y, n_subdivisions_z; 

    bool is_reduced_dt,
         enable_dam_debug,
         enable_brittle_LEFM,
         enable_elastoplastic_brittle_fracture; 
    
    static_assert(dim <= spacedim,
                  "The dimension <dim> of a PFC must be less than or "
                  "equal to the space dimension <spacedim> in which it lives.");
    
  };


/*
  template <int dim>
  class BodyForce : public Function<dim>
  {
  public:
    BodyForce();

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &  values) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override;
  };


  template <int dim>
  BodyForce<dim>::BodyForce()
    : Function<dim>(dim)
  {}


  template <int dim>
  inline void BodyForce<dim>::vector_value(const Point<dim> &,
                                           Vector<double> &values) const
  {
    AssertDimension(values.size(), dim);

    const double g   = 9.81;
    const double rho = 7700;

    values          = 0;
    values(dim - 1) = -rho * g;
  }



  template <int dim>
  void BodyForce<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>> &  value_list) const
  {
    const unsigned int n_points = points.size();

    AssertDimension(value_list.size(), n_points);

    for (unsigned int p = 0; p < n_points; ++p)
      BodyForce<dim>::vector_value(points[p], value_list[p]);
  }
*/
  
  template <int dim>
  class BodyPotential : public Function<dim>
  {
  public:
    BodyPotential();
    
    virtual void scalar_value(const Point<dim> &p,
                              double & values) const override;

    virtual void
    scalar_value_list(const std::vector<Point<dim>> &points,
                      std::vector<double> &  value_list) const override;
    
    const double amp;
    const double width;
    const Point<dim> Point0;
  };


  template <int dim>
  BodyPotential<dim>::BodyPotential()
    : Function<dim>(dim),
    amp(10.),
    width (5.),
    Point0(Point<dim> (0,0))
  {}


  template <int dim>
  inline void BodyPotential<dim>::scalar_value(const Point<dim> &p,
                                              double &values) const
  {
    //AssertDimension(values.size(), 1);

    
    double dist = p.distance(Point0);

    values = amp * std::exp(-(2.*dist*dist)/(width*width));
  }



  template <int dim>
  void BodyPotential<dim>::scalar_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<double> &  value_list) const
  {
    const unsigned int n_points = points.size();

    AssertDimension(value_list.size(), n_points);

    for (unsigned int p = 0; p < n_points; ++p)
      BodyPotential<dim>::scalar_value(points[p], value_list[p]);
  }

// Check this changed to load step 
// Only insert space function no time function: 
//  u = f(x,t)
//  Δu |x = (∂f/∂t) Δt
//  Here velocity is (∂f/∂t)
  template <int dim>
  class IncrementalBoundaryValues : public Function<dim>
  {
  public:
    IncrementalBoundaryValues(const double timestep,
                              const Tensor<1,dim> velocity);

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &  values) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override;

  private:
    const double timestep;
    const Tensor<1, dim> velocity;
  };


  template <int dim>
  IncrementalBoundaryValues<dim>::IncrementalBoundaryValues(
    const double timestep,
    const Tensor<1, dim> velocity)
    : Function<dim>(dim)
    , timestep(timestep)
    , velocity(velocity)
  {}



  template <int dim>
  void
  IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/,
                                               Vector<double> &values) const
  {
    AssertDimension(values.size(), dim);

    values    = 0;
    
    for(unsigned int i = 0; i < dim; i++)
      values(i) = timestep * velocity[i];
  }



  template <int dim>
  void IncrementalBoundaryValues<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>> &  value_list) const
  {
    const unsigned int n_points = points.size();

    AssertDimension(value_list.size(), n_points);

    for (unsigned int p = 0; p < n_points; ++p)
      IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]);
  }



/// PFC functions
// Constructor
  template <int dim, int spacedim>
  PFC<dim, spacedim>::PFC()
    : triangulation(MPI_COMM_WORLD)
    , dof_handler_dam(triangulation)
    , dof_handler_disp(triangulation)
    , fe_dam(1)
    , fe_disp(FE_Q<dim>(1), dim)
    //, quadrature_formula_s(fe_s.degree + 1)
    , quadrature_formula(fe_disp.degree + 1)
    , time(0.)
    , time_sin_start(0.)
    , timestep_no(0)
    , mpi_communicator(MPI_COMM_WORLD)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, this_mpi_process == 0)
    , enable_linear(true)
    , stored_time(false)
    , is_reduced_dt(false)
    , enable_dam_debug(false)
  {
  	
  	//std::cout << quadrature_formula.size();
  	
  	std::ifstream file("inputFiles/params.in");
  	std::string line;
  	while( std::getline( file, line ) )   
    {
      std::istringstream iss( line );

      std::string result;

      if(line.find_first_of("#") != std::string::npos)
        continue;

      if( std::getline( iss, result , '=') )
      {
        result = removeSpaces(result);
        if( result == "mesh_type" )
        {
          std::string token;
          std::getline( iss, token );            
          mesh_type = token;
        }
        
        if( result == "mesh_filename" )
        {
          std::string token;
          std::getline( iss, token );            
          mesh_filename = token;
        }
        
        if( result == "load_type" )
        {
          std::string token;
          std::getline( iss, token );            
          load_type = stoi(token);
        }
        
        if( result == "present_strain" )
        {
          std::string token;
          std::getline( iss, token );            
          present_strain = stod(token);
        }

        if( result == "strain_rate" )
        {
          std::string token;
          std::getline( iss, token );            
          strain_rate = stod(token);
        }

        if( result == "strain_min" )
        {
          std::string token;
          std::getline( iss, token );            
          strain_min = stod(token);
        }
        
        if( result == "strain_max" )
        {
          std::string token;
          std::getline( iss, token );            
          strain_max = stod(token);
        }
        
        if( result == "end_time" )
        {
          std::string token;
          std::getline( iss, token );            
          end_time = stod(token);
        }
        
        if( result == "default_dt" )
        {
          std::string token;
          std::getline( iss, token );            
          default_dt = stod(token);
        }
        
        dt = default_dt;

        if( result == "MU" )
        {
          std::string token;
          std::getline( iss, token );            
          MU = stod(token);
        }

        if( result == "Kn" )
        {
          std::string token;
          std::getline( iss, token );            
          Kn = stod(token);
        }
		
		    //Lambda
		    lambda = Kn - (2. * MU / double(spacedim));//Kn - 2.0 * MU / double(spacedim);
		
        if( result == "Eta" )
        {
          std::string token;
          std::getline( iss, token );            
          Eta = stod(token);
        }
        
        #ifndef PFC_DAMAGE
          Eta = 0;
        #endif
        
        if( result == "Ystrss" )
        {
          std::string token;
          std::getline( iss, token );            
          Ystrss = stod(token);
        }

        if( result == "H_pl" )
        {
          std::string token;
          std::getline( iss, token );            
          H_pl = stod(token);
        }

        if( result == "Ep_crit" )
        {
          std::string token;
          std::getline( iss, token );            
          Ep_crit = stod(token);
        }

        if( result == "Gc" )
        {
          std::string token;
          std::getline( iss, token );            
          Gc = stod(token);
        }

        if( result == "L_dam" )
        {
          std::string token;
          std::getline( iss, token );            
          L_dam = stod(token);
        }
        
        if( result == "enable_bounding_box" )
        {
          std::string token;
          std::getline( iss, token );            
          enable_bounding_box = stoi(token);
        }
        
        if( result == "refinement_level" )
        {
          std::string token;
          std::getline( iss, token );            
          refinement_level = stoi(token);
        }
        
        if( result == "enable_brittle_LEFM" )
        {
          std::string token;
          std::getline( iss, token );            
          enable_brittle_LEFM = stoi(token);
        }
        
        if( result == "enable_elastoplastic_brittle_fracture" )
        {
          std::string token;
          std::getline( iss, token );            
          enable_elastoplastic_brittle_fracture = stoi(token);
          if (enable_elastoplastic_brittle_fracture)
            enable_brittle_LEFM = false;
        }
        
        if( result == "outputFrequency" )
        {
          std::string token;
          std::getline( iss, token );            
          outputFrequency = stoi(token);
        }
        
        if( result == "N_stagger" )
        {
          std::string token;
          std::getline( iss, token );            
          N_stagger = stoi(token);
        }
        
        if( result == "N_NR_disp" )
        {
          std::string token;
          std::getline( iss, token );            
          N_NR_disp = stoi(token);
        }
        
        if( result == "N_NR_dam" )
        {
          std::string token;
          std::getline( iss, token );            
          N_NR_dam = stoi(token);
        }
        
        if( result == "Energy_tol" )
        {
          std::string token;
          std::getline( iss, token );            
          Energy_tol = stod(token);
        }
        
        if( result == "rhs_disp_tol" )
        {
          std::string token;
          std::getline( iss, token );            
          rhs_disp_tol = stod(token);
        }
        
        if( result == "rhs_dam_tol" )
        {
          std::string token;
          std::getline( iss, token );            
          rhs_dam_tol = stod(token);
        }
    		
        if( result == "lin_solve_disp_tol" )
        {
          std::string token;
          std::getline( iss, token );            
          lin_solve_disp_tol= stod(token);
        }
        
        if( result == "lin_solve_dam_tol" )
        {
          std::string token;
          std::getline( iss, token );            
          lin_solve_dam_tol= stod(token);
        }
        
        if( result == "d_eps_p_limit" )
        {
          std::string token;
          std::getline( iss, token );            
          d_eps_p_limit= stod(token);
        } 

        if( result == "d_dam_max_limit" )
        {
          std::string token;
          std::getline( iss, token );            
          d_dam_max_limit= stod(token);
        }
        
        if( result == "size_x" )
        {
          std::string token;
          std::getline( iss, token );            
          size_x = stod(token);
        }

        if( result == "size_y" )
        {
          std::string token;
          std::getline( iss, token );            
          size_y = stod(token);
        }

        if( result == "size_z" )
        {
          std::string token;
          std::getline( iss, token );            
          size_z = stod(token);
        }
      
        if( result == "n_subdivisions_x" )
        {
          std::string token;
          std::getline( iss, token );            
          n_subdivisions_x= stoi(token);
        }

        if( result == "n_subdivisions_y" )
        {
          std::string token;
          std::getline( iss, token );            
          n_subdivisions_y= stoi(token);
        }  

        if( result == "n_subdivisions_z" )
        {
          std::string token;
          std::getline( iss, token );            
          n_subdivisions_z= stoi(token);
        }
      

      }
    }
  	
  	
  	pcout<<"-----------------------------------------------------"<<std::endl; 
  	pcout<<"Starting Strain \t\t\t = "<<present_strain<< std::endl;
  	pcout<<"Ending Time \t\t\t = "<<end_time<<std::endl;
  	pcout<<"Strain rate  \t\t(s^-1)\t\t = "<<strain_rate<<std::endl;
  	pcout<<"Default timestep  \t(s)\t\t = "<<default_dt<<std::endl;
  	pcout<<"--------------------- Constants ---------------------"<<std::endl; 
  	pcout<<"Lame's constant: \n"
  		<<"\tLambda \t\t\t(MPa)\t = "<<lambda
  		<<"\n\tMu \t\t\t(MPa)\t = "<<MU<<std::endl; 
    pcout<<"Bulk Modulus \t\t\t(MPa)\t = "<<Kn<<std::endl; 
    pcout<<"Yield stress \t\t\t(MPa)\t = "<<Ystrss<<std::endl; 
    pcout<<"Hardening modulus \t\t(MPa)\t = "<<H_pl<<std::endl; 
    pcout<<"Material Fracture Toughness \t(MPa mm) = "<<Gc<<std::endl; 
    pcout<<"------------------ Input variables --------------------"<<std::endl; 
    pcout<<"Artificial residual stiffness\t\t\t = "<<Eta<<std::endl; 
    pcout<<"von Mises equivalent threshold plastic strain\t = "<<Ep_crit<<std::endl; 
    pcout<<"Transition zone parameter (mm)\t\t\t = "<<L_dam<<std::endl; 
    pcout<<"------------------- Mesh parameters -------------------"<<std::endl;
    pcout<<"Refinement level for given mesh\t\t = "<<refinement_level<<std::endl;
    pcout<<"Output Frequency\t\t\t = "<<outputFrequency<<std::endl;
    pcout<<"----------------- Control parameters ------------------"<<std::endl;
    pcout<<"Max. Staggered steps\t\t\t = "<<N_stagger<<std::endl;
    pcout<<"Max. NR steps for displacement\t\t = "<<N_NR_disp<<std::endl;
    pcout<<"Max. NR steps for damage\t\t = "<<N_NR_dam<<std::endl; 
    pcout<<"Energy error tolerance\t\t\t = "<<Energy_tol<<std::endl; 
    pcout<<"RHS displacement eq. tolerance\t\t = "<<rhs_disp_tol<<std::endl; 
    pcout<<"RHS damage eq. tolerance\t\t = "<<rhs_dam_tol<<std::endl; 
    pcout<<"-------------------------------------------------------"<<std::endl;
    
    
    const std::string st = "results/stress_strain.txt";
    ssFile.open(st);

    //#ifdef PFC_DEBUG
      const std::string outfilename = "debug/debug_list.txt";
                                       //Utilities::int_to_string(this_mpi_process, 3) +
                                       //".txt");
      outfile.open(outfilename);
      /*
      const std::string matrixfilename = ("debug/debug_matrix" + 
                                       Utilities::int_to_string(this_mpi_process, 3) +
                                       ".txt");
      matrixfile.open(matrixfilename);
      
      const std::string cellmatrixfilename = ("debug/debug_cellmatrix" + 
                                       Utilities::int_to_string(this_mpi_process, 3) +
                                       ".txt");
      cellmatrixfile.open(cellmatrixfilename);*/
    //   if(this_mpi_process == 0)
    //     debug_stagger_list << stag <<std::endl;
    //#endif
        
  }

// Destructor
  template <int dim, int spacedim>
  PFC<dim, spacedim>::~PFC()
  {
    dof_handler_disp.clear();
    dof_handler_dam.clear();
    
    ssFile.close();
    outfile.close();
    //matrixfile.close();
    //cellmatrixfile.close();
    
  }

// Find bounding box, useful when model oriented in x, y, z, do not need to enter size x,y,z
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::calc_bounding_box()
  {
    Point<dim> minbound, maxbound;
    if (dim == 2)
    {
      minbound = Point<dim> (std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max());
      maxbound = Point<dim> (std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::lowest());
    }
    else if (dim ==3)
    {
      minbound = Point<dim> (std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max());
      maxbound = Point<dim> (std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::lowest());
    }
    else
    {
      pcout << "\n Bounding box not defined for this dimension"<< std::endl;
    }
    
    for (auto cell : triangulation.active_cell_iterators())
    for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
    {
      const Point<dim> &p = cell->vertex(v);
      
      for (unsigned int d = 0; d < dim; d++)
      {
        minbound[d] = std::min(minbound[d], p[d]);
        maxbound[d] = std::max(maxbound[d], p[d]);
      }
    }
    
    if (dim == 2)
    {
      size_x = maxbound[0] - minbound[0];
      size_y = maxbound[1] - minbound[1];
      size_z = 1.;
    }
    else
    {
      size_x = maxbound[0] - minbound[0];
      size_y = maxbound[1] - minbound[1];
      size_z = maxbound[2] - minbound[2];
    }
  }


// Create coarse grid
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::make_grid()
  {
  	
    
    double bottom_y(0.), top_y(0.);
    
    if (mesh_type == "rectangular_plate")
    {
      Point<dim> corner1, corner2;
      std::vector<unsigned int> repetitions;
      
      if (dim == 3)
      {
        
        corner1 = Point<dim>(0, 0, 0);
        corner2 = Point<dim>(size_x, size_y , size_z);
        
        repetitions.push_back(n_subdivisions_x); 
        repetitions.push_back(n_subdivisions_y); 
        repetitions.push_back(n_subdivisions_z);
        
      }
      else // For 2D
      {
        corner1 = Point<dim>(0, 0);
        corner2 = Point<dim>(size_x, size_y);
        
        repetitions.push_back(n_subdivisions_x); 
        repetitions.push_back(n_subdivisions_y);
        
      }
      
      bottom_y = corner1(1);
      top_y = bottom_y + size_y;
      
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                repetitions  ,
                                                corner1      ,
                                                corner2      ,
                                                false        );
    }
    
    if (mesh_type == "I_bar")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 20. / std::pow(2., refinement_level);
	    
	    size_x = 60.;
	    size_y = 120.;
	    
	    bottom_y = 0.;
	    top_y = bottom_y + size_y;
	    
	    grid_Ibar<dim>(triangulation, (n_subdivisions_z + 1), size_z);
    }
  	
  	if (mesh_type == "asymetrical_notch")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 4. / std::pow(2., refinement_level);
	    
	    size_x = 18.;
	    size_y = 50.;
	    
	    bottom_y = -25.;
	    top_y = bottom_y + size_y;
	    
	    grid_asymm_notches<dim>(triangulation, (n_subdivisions_z + 1), size_z);
    }
    
    if (mesh_type == "symmetrical_notch_large")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 9. / std::pow(2., refinement_level);
	    
	    size_x = 18.;
	    size_y = 50.;
	    
	    bottom_y = -25.;
	    top_y = bottom_y + size_y;
	    
	    grid_symm_notches_r<dim>(triangulation,
	                             5.,
	                             4.,
	                             (n_subdivisions_z + 1),
	                             size_z);
    }
    
    if (mesh_type == "symmetrical_notch_small")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 9. / std::pow(2., refinement_level);
      
	    size_x = 18.;
	    size_y = 50.;
	    
	    bottom_y = -25.;
	    top_y = bottom_y + size_y;
	    
	    grid_symm_notches_r<dim>(triangulation,
	                             2.5,
	                             6.5,
	                             (n_subdivisions_z + 1),
	                             size_z);
    }
    
    if (mesh_type == "single_notch")
    {
      
      Triangulation<2> triatmp, triatmp_notch;
      
      if (dim == 3)
	      size_z = n_subdivisions_z * 5. / std::pow(2., refinement_level);
	    
	    size_x = 10.;
	    size_y = 10.;
	    
	    bottom_y = -5.;
	    top_y = bottom_y + size_y;
      double notch_r = 0.005;
	    
	    GridGenerator::hyper_cube(triatmp, -5, 5);
	    
      
      for (unsigned int i = 0; i< refinement_level; i++)
      {
        for (const auto &cell : triatmp.active_cell_iterators())
          cell->set_refine_flag(RefinementCase<2>::cut_xy);
        //GridGenerator::extrude_triangulation(triatmp, 2, 20, triangulation, true);
        triatmp.execute_coarsening_and_refinement();
      }
      
      std::set<typename Triangulation<2>::active_cell_iterator> cells_to_remove;
      
      for (const auto &cell : triatmp.active_cell_iterators())
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          if ((cell->center()[1] < (cell->diameter()/2.)) && cell->center()[1] > 0. && cell->center()[0] < .0)
            cells_to_remove.insert(cell);
      
      // Just one line of cells
      //    if (std::fabs(cell->vertex(v)[1]) < notch_r && cell->vertex(v)[0] < .0) OR
      // (cell->center()[1] < (cell->diameter()/2.)) && cell->center()[1] > 0. && cell->center()[0] < .0 
      //      cells_to_remove.insert(cell);
      
      GridGenerator::create_triangulation_with_removed_cells(triatmp,
                                                             cells_to_remove,
                                                             triatmp_notch);
      if (dim == 3)
      GridGenerator::extrude_triangulation(triatmp_notch,
	                                           (n_subdivisions_z + 1),
	                                           size_z,
	                                           triangulation);
      else
        triangulation.copy_triangulation(triatmp_notch);

      for (const auto &cell : triangulation.active_cell_iterators())
        for (const auto &face : cell->face_iterators())
        {
          
          if (std::fabs(face->center()[1] - top_y) < 1e-12)
            face->set_boundary_id(11);
          //if (std::fabs(face->center()[0] + 5) < 1e-12)
          //  face->set_boundary_id(12);//5.
          if (std::fabs(face->center()[1] - bottom_y) < 1e-12)
            face->set_boundary_id(13);
        }
      
    }
    
    if (mesh_type == "double_diagonal_notch")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 5. / std::pow(2., refinement_level);
	    
	    size_x = 7.23; // Approximate 10 - 1 - 2.5/sqrt(2)
	    size_y = 10.;
	    
	    bottom_y = 0.;
	    top_y = bottom_y + size_y;
	    
	    grid_2_corner_notches<dim>(triangulation,
	                               (n_subdivisions_z + 1),
	                               size_z);
      
      for (const auto &cell : triangulation.active_cell_iterators())
        for (const auto &face : cell->face_iterators())
        {
          
        	if (std::fabs(face->center()[0] - 0.) < 1e-12)
        		face->set_boundary_id(11);
          //
          if (std::fabs(face->center()[0] - 10.) < 1e-12)
        		face->set_boundary_id(13);
          //
          //if (std::fabs(face->center()[1] + 5.) < 1e-12)
        	//	face->set_boundary_id(13);

          if (std::fabs(face->center()[1] - top_y) < 1e-12)
            face->set_boundary_id(11);
          //if (std::fabs(face->center()[0] + 5) < 1e-12)
          //  face->set_boundary_id(12);//5.
          if (std::fabs(face->center()[1] - bottom_y) < 1e-12)
            face->set_boundary_id(13);
        }
      
      for (unsigned int i = 0; i< refinement_level; i++)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_refine_flag(RefinementCase<dim>::cut_xy);
        //GridGenerator::extrude_triangulation(triatmp, 2, 20, triangulation, true);
        triangulation.execute_coarsening_and_refinement();
      }
    }
    
    if (mesh_type == "CT")
    {
      if (dim == 3)
	      size_z = n_subdivisions_z * 9. / std::pow(2., refinement_level);
	    
	    size_x = 63.8;
	    size_y = 60.96;
	    
	    bottom_y = -30.48;
	    top_y = bottom_y + size_y;
	    
	    grid_CT<dim>(triangulation,
                   (n_subdivisions_z + 1),
                   size_z);
    }
    
    if (mesh_type == "custom_mesh") // Only reads .msh files
    {
      /// Reads mesh from file and creates the custom mesh,
      /// currently only works for gmsh .msh files
      
	    GridIn<dim> gridin;
	    gridin.attach_triangulation(triangulation);
	    std::ifstream msh_f(mesh_filename);
	    gridin.read_msh(msh_f);
	    
	    if (enable_bounding_box)
	      calc_bounding_box();
    }
    
    ////////////////////// XY refinement
    if (mesh_type != "single_notch" &&
        mesh_type != "double_diagonal_notch" &&
        mesh_type != "CT" &&
        mesh_type != "custom_mesh")
    {
        
    	for (const auto &cell : triangulation.active_cell_iterators())
        for (const auto &face : cell->face_iterators())
        {
          
        	//if (std::fabs(face->center()[1] - 5.) < 1e-12)
        	//	face->set_boundary_id(11);
          //
          //if (std::fabs(face->center()[0] - 5) < 1e-12)
        	//	face->set_boundary_id(12);
          //
          //if (std::fabs(face->center()[1] + 5.) < 1e-12)
        	//	face->set_boundary_id(13);

          if (std::fabs(face->center()[1] - top_y) < 1e-12)
            face->set_boundary_id(11);
          //if (std::fabs(face->center()[0] + 5) < 1e-12)
          //  face->set_boundary_id(12);//5.
          if (std::fabs(face->center()[1] - bottom_y) < 1e-12)
            face->set_boundary_id(13);
        }
      
      for (unsigned int i = 0; i< refinement_level; i++)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_refine_flag(RefinementCase<dim>::cut_xy);
        //GridGenerator::extrude_triangulation(triatmp, 2, 20, triangulation, true);
        triangulation.execute_coarsening_and_refinement();
      }
      
    }
    ///////////////////////
     
  	//triangulation.refine_global(2);
  	//triangulation.refine_global(refinement_level);
  	
  	double suggested_L_dam = L_dam;
  	
  	for (const auto &cell : triangulation.active_cell_iterators())
  	  if (cell->is_locally_owned())
  	  {
  	    double cell_diameter = cell->diameter();
  	    
  	    if (cell_diameter > suggested_L_dam)
          suggested_L_dam = cell_diameter;
  	  }
  	
  	suggested_L_dam= Utilities::MPI::max(suggested_L_dam, mpi_communicator);
  	
  	if (suggested_L_dam > L_dam)
  	  pcout << "\t Prescribed L_dam: "<< suggested_L_dam <<"\tInputted L_dam: "<< L_dam <<"\n";
  }



// Setup Systems
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::setup_system(const int flg)
  {
       
    if(!flg)// This part runs everytime system_setup is called and iteration is 0
    {
            
      {
        dof_handler_disp.distribute_dofs(fe_disp);
        // Renumbering
        //DoFRenumbering::component_wise(dof_handler_u, stokes_sub_blocks);
        
        locally_owned_dofs_disp = dof_handler_disp.locally_owned_dofs();
        locally_relevant_dofs_disp =
          DoFTools::extract_locally_relevant_dofs(dof_handler_disp);
          
        hanging_node_constraints_disp.clear();
        
        DoFTools::make_hanging_node_constraints(dof_handler_disp,
                                                hanging_node_constraints_disp);
        
        
        //{// Fixed and sliding node constraint // Change when using other geometries
        //  // Fixed node 
        //  Point<dim> fixed_node_unit_cell(0., 0.);
        //  for (const auto &cell : dof_handler_disp.active_cell_iterators()) 
        //    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        //      if (cell->is_locally_owned())// Check if the vertex is the target vertex
        //        if (fixed_node_unit_cell.distance (cell->vertex(i)) < 1e-2 * cell->diameter())
        //         {
        //           unsigned int dof1=cell->vertex_dof_index(i,0);
        //           unsigned int dof2=cell->vertex_dof_index(i,1);
        //        
        //           hanging_node_constraints_disp.add_line(dof1);
        //           hanging_node_constraints_disp.add_line(dof2);
        //         }
        //
        //   //Sliding node
        //   Point<dim> slide_node_unit_cell(60., 0.);
        //   for (const auto &cell : dof_handler_disp.active_cell_iterators()) 
        //     for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        //       if (cell->is_locally_owned())// Check if the vertex is the target vertex
        //         if (slide_node_unit_cell.distance (cell->vertex(i)) < 1e-2 * cell->diameter())
        //         {
        //           unsigned int dof1=cell->vertex_dof_index(i,1);
        //        
        //           hanging_node_constraints_disp.add_line(dof1);
        //         }
        // }


        //COMMENTED FOR I-BAR
        // Point<dim> left_bottom_back(-size_x/2.0,-size_y/2.0,-size_z/2.0); 
        // Point<dim> right_bottom_back(-size_x/2.0,-size_y/2.0,size_z/2.0); 
        // Point<dim> left_bottom_front(size_x/2.0,-size_y/2.0,-size_z/2.0); 

        // unsigned int vertices_per_cell=GeometryInfo<dim>::vertices_per_cell;

        // // Loop over each locally owned cell
        // //typename DoFHandler<dim>::active_cell_iterator cell= dof_handler->begin_active(), endc = dof_handler->end();

        // //for (; cell!=endc; ++cell)
        // //

        // for (const auto &cell : dof_handler_disp.active_cell_iterators())
        // {  
        //   for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        //   {
        //    if (cell->is_locally_owned())
        //    {
        //       // Check if the vertex is the target vertex
        //       if (left_bottom_back.distance (cell->vertex(i)) < 1e-2 * cell->diameter())
        //        {
        //         // Loop through the list of components with rigid body modes and add an inhomogeneous constraint for each
        //         //for (unsigned int component_num = 0; component_num < rigidBodyModeComponents.size(); component_num++){
        //           //std::cout<<"\n Found...."; 
        //           unsigned int dof1=cell->vertex_dof_index(i,0);
        //           unsigned int dof2=cell->vertex_dof_index(i,1);
        //           unsigned int dof3=cell->vertex_dof_index(i,2);
                  
        //           hanging_node_constraints_disp.add_line(dof1);
        //           hanging_node_constraints_disp.add_line(dof2);
        //           hanging_node_constraints_disp.add_line(dof3);
        //           //constraints_rbc.set_inhomogeneity(dof1,0.0);
        //           //solution(dof1) = 0.0;
        //         }
        //     }
        //   }
        // }


        // for (const auto &cell : dof_handler_disp.active_cell_iterators())
        // {  
        //   for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        //   {
        //    if (cell->is_locally_owned())
        //    {
        //       // Check if the vertex is the target vertex
        //       if (right_bottom_back.distance (cell->vertex(i)) < 1e-2 * cell->diameter())
        //        {
        //         // Loop through the list of components with rigid body modes and add an inhomogeneous constraint for each
        //         //for (unsigned int component_num = 0; component_num < rigidBodyModeComponents.size(); component_num++){
        //           //std::cout<<"\n Found...."; 
        //           unsigned int dof1=cell->vertex_dof_index(i,0);
        //           unsigned int dof2=cell->vertex_dof_index(i,1);
                  
        //           hanging_node_constraints_disp.add_line(dof1);
        //           hanging_node_constraints_disp.add_line(dof2);
        //           //constraints_rbc.set_inhomogeneity(dof1,0.0);
        //           //solution(dof1) = 0.0;
        //         }
        //     }
        //   }
        // }

        // for (const auto &cell : dof_handler_disp.active_cell_iterators())
        // {  
        //   for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        //   {
        //    if (cell->is_locally_owned())
        //    {
        //       // Check if the vertex is the target vertex
        //       if (left_bottom_front.distance (cell->vertex(i)) < 1e-2 * cell->diameter())
        //        {
        //         // Loop through the list of components with rigid body modes and add an inhomogeneous constraint for each
        //         //for (unsigned int component_num = 0; component_num < rigidBodyModeComponents.size(); component_num++){
        //           //std::cout<<"\n Found...."; 
        //           unsigned int dof1=cell->vertex_dof_index(i,1);
                  
        //           hanging_node_constraints_disp.add_line(dof1);
        //           //constraints_rbc.set_inhomogeneity(dof1,0.0);
        //           //solution(dof1) = 0.0;
        //         }
        //     }
        //   }
        // }


        
        hanging_node_constraints_disp.close();
        
      }
      
      {
        dof_handler_dam.distribute_dofs(fe_dam);
         
        locally_owned_dofs_dam = dof_handler_dam.locally_owned_dofs();
        locally_relevant_dofs_dam =
          DoFTools::extract_locally_relevant_dofs(dof_handler_dam);
        
        hanging_node_constraints_dam.clear();
        
        DoFTools::make_hanging_node_constraints(dof_handler_dam,
                                                hanging_node_constraints_dam);
        
        hanging_node_constraints_dam.close();
      }
      
      disp.reinit(dof_handler_disp.n_dofs());
      disp_NR.reinit(dof_handler_disp.n_dofs());
      
      disp = 0;//Final displacement after everything converges
		  disp_NR = 0; //Displacement for each Staggered solution after NR blocks for disp has converged
		  dam.reinit(dof_handler_dam.n_dofs());
      dam_NR.reinit(dof_handler_dam.n_dofs());
      
      dam = 1;//Final damage after everything converges
		  dam_NR = 1; //damage for each Staggered solution after NR blocks for damage has converged
		  
		  
		  // This is notch condition, taking damage = 0 at notch dofs
		  /*
		  if (mesh_type == "single_notch")
      {
		    auto cell = triangulation.begin();
		    double cell_diameter= cell->diameter();
		    MappingQ<dim, dim> mapping(1);
		    std::vector<Point<dim>> support_points(dof_handler_dam.n_dofs());
		    
		    
		    DoFTools::map_dofs_to_support_points<dim, dim>(mapping,
                                                       dof_handler_dam,
                                                       support_points);
		    for (unsigned int i = 0; i < dof_handler_dam.n_dofs(); ++i)
          if(locally_owned_dofs_dam.is_element(i) &&
             support_points[i](1)<=cell_diameter && 
             support_points[i](1)>=0. &&
             support_points[i](0)<1e-12)
          {
            dam(i) = 0.00;
            dam_NR(i) = 0.00; 
          }
      }*/
		  
		  
		  if ((mesh_type == "custom_mesh") && (mesh_filename == "inputFiles/Notchsmall"))
      {
		    
		    MappingQ<dim, dim> mapping(1);
		    
		    
		    std::vector<bool> boundary_dofs(dof_handler_dam.n_dofs());
		    DoFTools::extract_boundary_dofs(dof_handler_dam,
                                    {},//ComponentMask(),
                                    boundary_dofs,
                                    {});
		    
		    
		    std::vector<Point<dim>> support_points(dof_handler_dam.n_dofs());
		    DoFTools::map_dofs_to_support_points<dim, dim>(mapping,
                                                       dof_handler_dam,
                                                       support_points);
		    
		    for (unsigned int i = 0; i < dof_handler_dam.n_dofs(); ++i)
          if (boundary_dofs[i] == true)
            if(locally_owned_dofs_dam.is_element(i) &&
               std::abs(support_points[i](1))<=0.25)
            {
              dam(i) = 0.98;
              dam_NR(i) = 0.98; 
            }
      }
		  
    }

    // This part runs everytime system_setup is called
    DynamicSparsityPattern sparsity_pattern_dam(locally_relevant_dofs_dam);
    DoFTools::make_sparsity_pattern(dof_handler_dam,
                                    sparsity_pattern_dam,
                                    hanging_node_constraints_dam,
                                    /*keep constrained dofs*/ false);
    SparsityTools::distribute_sparsity_pattern(sparsity_pattern_dam,
                                               locally_owned_dofs_dam,
                                               mpi_communicator,
                                               locally_relevant_dofs_dam);
    system_matrix_dam.reinit(locally_owned_dofs_dam,
                         locally_owned_dofs_dam,
                         sparsity_pattern_dam,
                         mpi_communicator);
    system_rhs_dam.reinit(locally_owned_dofs_dam, mpi_communicator);
    
    d_dam.reinit(dof_handler_dam.n_dofs());
    
    DynamicSparsityPattern sparsity_pattern_disp(locally_relevant_dofs_disp);
    DoFTools::make_sparsity_pattern(dof_handler_disp,
                                    sparsity_pattern_disp,
                                    hanging_node_constraints_disp,
                                    /*keep constrained dofs*/ false);
    SparsityTools::distribute_sparsity_pattern(sparsity_pattern_disp,
                                               locally_owned_dofs_disp,
                                               mpi_communicator,
                                               locally_relevant_dofs_disp);
    system_matrix_disp.reinit(locally_owned_dofs_disp,
                         locally_owned_dofs_disp,
                         sparsity_pattern_disp,
                         mpi_communicator);
    system_rhs_disp.reinit(locally_owned_dofs_disp, mpi_communicator);
    system_rxn.reinit(locally_owned_dofs_disp, mpi_communicator);
    d_disp.reinit(dof_handler_disp.n_dofs());
  }

  // Setup Qp history //Problem here?
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::setup_Qp_history()
  {
    
    triangulation.clear_user_data();
    
    {
      std::vector<PointHistory<spacedim>> tmp;  // PointHistory may not require dim
      quadrature_point_history.swap(tmp);   
    }
    quadrature_point_history.resize(
      triangulation.n_locally_owned_active_cells() * quadrature_formula.size());
    
    
    unsigned int history_index = 0;
    for (auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell->set_user_pointer(&quadrature_point_history[history_index]);
          history_index += quadrature_formula.size();
        }

    Assert(history_index == quadrature_point_history.size(),
           ExcInternalError());
    
    for (auto &cell : dof_handler_disp.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        
        // Next, get a pointer to the quadrature point history data local to
        // the present cell, and, as a defensive measure, make sure that
        // this pointer is within the bounds of the global array:
        
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());
        
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          
          Tensor<4, spacedim> tmp_Cijkl = 
            calc_strss_strain_tensor_constant (local_quadrature_points_history[q].eps_e_NR,
                                               local_quadrature_points_history[q].g_hat);// correct digits but no decimal, although it shows it to be double, eps double, cout prints 6 digits only, precision problems can occur during printing, don't get confused,
          //pcout<<typeid(MU).name()<<std::endl;
          //std::cout.precision(10);
          //pcout<<local_quadrature_points_history[q].eps_e[0]<<std::endl;
          local_quadrature_points_history[q].Psi_plus_max_NR =
            Psi_plus(local_quadrature_points_history[q].eps_e_NR);
          local_quadrature_points_history[q].Psi_plus_max =
            local_quadrature_points_history[q].Psi_plus_max_NR;
          
                      
          local_quadrature_points_history[q].strss_NR = 
            calc_strss_from_strain (tmp_Cijkl,
                                    local_quadrature_points_history[q].eps_e_NR);
          
          //local_quadrature_points_history[q].Qindex_u = q;
        }
      }
      
    // damage Qp
    
    
    FEValues<dim> fe_values(fe_dam,
                            quadrature_formula,
                            update_values |
                              update_quadrature_points );

    //const unsigned int n_q_points    = quadrature_formula.size();
    
    std::vector<double> val_dam(quadrature_formula.size());
    
    for (auto &cell : dof_handler_dam.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        
        fe_values.reinit(cell);
        fe_values.get_function_values(dam_NR, val_dam);
        
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());
        
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          // Here there should be plastic integration part
          if (enable_elastoplastic_brittle_fracture || enable_brittle_LEFM)
            local_quadrature_points_history[q].acc_eps_p_NR = Ep_crit;
          
          local_quadrature_points_history[q].acc_eps_p =      // This is fine
            local_quadrature_points_history[q].acc_eps_p_NR;
          
          local_quadrature_points_history[q].g_hat= std::pow(val_dam[q], 2. *
                     local_quadrature_points_history[q].acc_eps_p/Ep_crit) + Eta;// correct
          //local_quadrature_points_history[q].Qindex_s = q;
        }
      }
      
    /*
    for (auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());
        
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          pcout << "Error" << local_quadrature_points_history[q].Qindex_s<<"\t" 
                            <<local_quadrature_points_history[q].Qindex_u<<"\n";
          if (local_quadrature_points_history[q].Qindex_s !=
               local_quadrature_points_history[q].Qindex_u)
          {  pcout << "Error" << local_quadrature_points_history[q].Qindex_s<<"\t" <<local_quadrature_points_history[q].Qindex_u<<"\n";
          }
        }
      }
        
    */  
  }


  template <int dim, int spacedim>
  void PFC<dim, spacedim>::assemble_system_disp( const int iteration)
  {
    
    Tensor <1, dim> vel;
    vel[1] = velocity;
    
    // Apply BC on disp, BC may be outside staggered
    system_rhs_disp    = 0;
    system_matrix_disp = 0;

    FEValues<dim> fe_values(fe_disp,
                            quadrature_formula,
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe_disp.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); //(Ai X Bk)
    Vector<double>     cell_residual(dofs_per_cell); // (Ai)

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    //BodyForce<dim>              body_force;
    //std::vector<Vector<double>> body_force_values(n_q_points,
    //                                              Vector<double>(dim));

   std::vector<std::vector<Tensor<1, dim>>> disp_grad(
    quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim));
  
    for (const auto &cell : dof_handler_disp.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          
          
          cell_matrix = 0;
          cell_residual    = 0;

          fe_values.reinit(cell);
          fe_values.get_function_gradients(disp_NR, disp_grad);
          
          const PointHistory<spacedim> *local_quadrature_points_history =
            reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
          Assert(local_quadrature_points_history >=
                   &quadrature_point_history.front(),
                 ExcInternalError());
          Assert(local_quadrature_points_history <=
                   &quadrature_point_history.back(),
                 ExcInternalError());
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            //ssFile << "\n\t Qp:" << std::endl;
            
            const Tensor<2, dim> tmp_strain = get_strain(disp_grad[q]);
            Tensor<2, spacedim> strain;
            
              strain = tensordimtospacedim<dim, spacedim> ( tmp_strain);
            
            double yield_fn;
            
            Tensor<2, spacedim> strss_pc, eps_e_pc, eps_p_pc;
            Tensor<4, spacedim> Cijkl_pc;
            
            double alpha_pc, del_lambda_pc;
            
            if ( !enable_brittle_LEFM )
              predictor_corrector ( strain,
                                    local_quadrature_points_history[q].eps_e,
                                    local_quadrature_points_history[q].eps_p,
                                    local_quadrature_points_history[q].alpha,
                                    local_quadrature_points_history[q].g_hat,
                                    yield_fn,
                                    del_lambda_pc,
                                    alpha_pc,
                                    strss_pc,
                                    eps_e_pc,
                                    eps_p_pc,
                                    Cijkl_pc);
           
           if ( enable_brittle_LEFM )
           {
             eps_e_pc = strain;
             Cijkl_pc = calc_strss_strain_tensor_constant(
                         eps_e_pc,
                         local_quadrature_points_history[q].g_hat);
             strss_pc = calc_strss_from_strain (Cijkl_pc, eps_e_pc);
           }
            
            for (unsigned int m = 0; m < dofs_per_cell; ++m)
            { 
              const unsigned int component_m = fe_disp.system_to_component_index(m).first;
              
              for (unsigned int n = 0; n < dofs_per_cell; ++n)
              {
                
                const unsigned int component_n = fe_disp.system_to_component_index(n).first;
                /*const Tensor<2, dim>
                eps_phi_m = get_strain(fe_values, m, q);*/
                const Tensor<2, dim> eps_phi_n = get_strain(fe_values, n, q);
                const Tensor<2, spacedim> eps_n = 
                      tensordimtospacedim<dim, spacedim> (get_strain(fe_values, n, q));
                const Tensor<2, spacedim> d_sigma_n = double_contract(Cijkl_pc, eps_n);
                //std::cout << eps_phi_m << eps_phi_n  ;
                for(unsigned int J=0;J<dim;J++)
                  //for(unsigned int L=0;L<dim;L++)
                  {
                  
                //std::cout<<local_quadrature_points_history[q].Cijkl<<std::endl;
                //cell_matrix.print_formatted(ssFile);// Not working
                  /* cell_matrix(m, n) += (fe_values.shape_grad(m, q)[J]  *          
                                        Cijkl_pc[component_m][J][component_n][L]*
                                        fe_values.shape_grad(n, q)[L]
                                        ) * 
                                       fe_values.JxW(q); */ 
                  cell_matrix(m, n) += (fe_values.shape_grad(m, q)[J] *
                                        d_sigma_n[component_m][J]
                                       ) * 
                                      fe_values.JxW(q); 
                  
                  }
                
                //std::cout<<local_quadrature_points_history[q].Cijkl<<std::endl;
                
                
                }
            } 
         
            for (unsigned int m = 0; m < dofs_per_cell; ++m)
            {
              //const unsigned int component_i =
              //  fe_u.system_to_component_index(i).first;
              const unsigned int component_m =
                fe_disp.system_to_component_index(m).first;

                for(int J=0;J<dim;J++)//if(component_i==j)
                //for(unsigned int O=0;O<dim;O++)
                {
              
                  cell_residual(m) -= 
                    (//body_force_values[q](component_i) *
                      // fe_values.shape_value(i, q) 
                       fe_values.shape_grad(m, q)[J] *
                       strss_pc[component_m][J]) *
                       fe_values.JxW(q);
                  //cell_residual(m) -= 
                  // (//body_force_values[q](component_i) *
                  //   // fe_values.shape_value(i, q) 
                  //  get_strain(fe_values, m, q)[O][J] *
                  //  strss_pc[O][J]) *
                  //  fe_values.JxW(q);

              
                }
            }
          }

          cell->get_dof_indices(local_dof_indices);

          hanging_node_constraints_disp.distribute_local_to_global(cell_matrix,
                                                              cell_residual,
                                                              local_dof_indices,
                                                              system_matrix_disp,
                                                              system_rhs_disp);
        }

    system_matrix_disp.compress(VectorOperation::add);
    system_rhs_disp.compress(VectorOperation::add);

    
    //pcout << "Is Symmetric? "<< system_matrix_disp.is_symmetric() <<std::endl;

    const FEValuesExtractors::Scalar          x_component(0);
    const FEValuesExtractors::Scalar          y_component(1);
    std::map<types::global_dof_index, double> boundary_values;
    // interpolate_boundary_values adds but does not remove
    /*VectorTools::interpolate_boundary_values(dof_handler_disp,
                                             13,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values,
                                             fe_disp.component_mask(y_component));*/
    VectorTools::interpolate_boundary_values(dof_handler_disp,
                                             13,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values);
    //pcout<< boundary_values.size()<<std::endl;
    
    if (!iteration)
    {
      /*VectorTools::interpolate_boundary_values(dof_handler_disp,
                                               11,
                                               IncrementalBoundaryValues<dim>(dt, vel),
                                               boundary_values,
                                             fe_disp.component_mask(y_component));*/
      
      VectorTools::interpolate_boundary_values(dof_handler_disp,
                                               11,
                                               IncrementalBoundaryValues<dim>(dt, vel),
                                               boundary_values);
      
      VectorTools::interpolate_boundary_values(dof_handler_disp,
                                               12,
                                               IncrementalBoundaryValues<dim>(dt, vel),
                                               boundary_values);
      /*pcout << boundary_values.size(); // Point test
      for (const auto & [key, value] : boundary_values)
        std::cout << '[' << key << "] = " << value << "; ";
      
      boundary_values.erase(6);
      
      for (const auto & [key, value] : boundary_values)
        std::cout << '[' << key << "] = " << value << "; ";*/
      
    }
    else
    {
      /*VectorTools::interpolate_boundary_values(dof_handler_disp,
                                             11,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values,
                                             fe_disp.component_mask(y_component));*/
      VectorTools::interpolate_boundary_values(dof_handler_disp,
                                             11,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_disp,
                                             12,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values);
      //boundary_values.erase(6); // Point test
    }
    
    PETScWrappers::MPI::Vector tmp(locally_owned_dofs_disp, mpi_communicator);
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix_disp, tmp, system_rhs_disp, false);
    d_disp = tmp;
  }
  
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::assemble_system_dam()
  {   
    system_rhs_dam    = 0;
    system_matrix_dam = 0;

    FEValues<dim> fe_values(fe_dam,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe_dam.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_residual(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double> val_dam(quadrature_formula.size());
    std::vector<Tensor<1, dim>> grad_dam(quadrature_formula.size());

//        BodyPotential<dim>              body_force;
//        std::vector<Vector<double>> body_force_values(n_q_points,
//                                                      Vector<double>(dim));
      //std::vector<std::vector<Tensor<1, dim>>> solution_grad(
  //quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim));

    
    for (const auto &cell : dof_handler_dam.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_residual    = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(dam_NR, val_dam);
        fe_values.get_function_gradients(dam_NR, grad_dam);
        
        const PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {  
          double p_k = (local_quadrature_points_history[q_point].acc_eps_p +
            
              local_quadrature_points_history[q_point].alpha_NR -
              local_quadrature_points_history[q_point].alpha)
                / Ep_crit;
          
          if (enable_elastoplastic_brittle_fracture || enable_brittle_LEFM)
            p_k = 1.;
          
          double Psi_plus_max = local_quadrature_points_history[q_point].Psi_plus_max_NR;
          
          //ssFile << "\t pk:"<< p_k<< std::endl;
          
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               
              double grad_product(0.);
              for(unsigned int K=0;K<dim;K++)
                grad_product +=fe_values.shape_grad(i, q_point)[K] *
                                    fe_values.shape_grad(j, q_point)[K] ;
              
              
              cell_matrix(i, j) += (grad_product//
                                    
                                    +
                                    
                                    fe_values.shape_value(i, q_point) * 
                                    
                                    ((1. / (4.*L_dam*L_dam)) +
                                    (std::pow (val_dam[q_point], (2. * p_k - 2.)) / (Gc * L_dam))
                                    * (2. * p_k - 1.) * p_k * Psi_plus_max) * 
                                    
                                    fe_values.shape_value(j, q_point)) *
                                             
                                   fe_values.JxW(q_point); //
            }

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            
            double grad_dam_product(0.);
            
            for(unsigned int K=0;K<dim;K++)
              grad_dam_product += fe_values.shape_grad(i, q_point)[K] * grad_dam[q_point][K];
            
            cell_residual(i) -=  (grad_dam_product -
            
                                  fe_values.shape_value(i, q_point) * 
                                  (
                                   ( (1. - val_dam[q_point])/ (4. * L_dam * L_dam)) -
                                   ( std::pow (val_dam[q_point], (2. * p_k - 1.) ) / (Gc * L_dam))
                                   * p_k * Psi_plus_max
                                   )
                                   ) *
                                  fe_values.JxW(q_point);
          }
        }

        #ifdef PFC_DEBUG
        if ((cell->center()[1] < 7.) && (cell->center()[1] > -7.))
        {
          outfile << "\n\tCell Id: \t"
                  << cell->id();
                    /*<< "\t Cell center: \t"
                    <<cell->center();*/
          double max_val = cell_residual.linfty_norm();
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            outfile << "\n\t DOF index: \t"
                    << i
                    << "\t Cell ratio: \t"
                    <<(cell_residual(i)/max_val)
                    << "\t Cell residue: \t"
                    <<cell_residual(i);
        }
        #endif
        cell->get_dof_indices(local_dof_indices);

        hanging_node_constraints_dam.distribute_local_to_global(cell_matrix,
                                                            cell_residual,
                                                            local_dof_indices,
                                                            system_matrix_dam,
                                                            system_rhs_dam);
      }

    system_matrix_dam.compress(VectorOperation::add);
    system_rhs_dam.compress(VectorOperation::add);


    //pcout << "Is Symmetric? "<< system_matrix_dam.is_symmetric() <<std::endl;
    
    
    //std::map<types::global_dof_index, double> boundary_values;
    
    // Only Neumann no dirichlet BC
   
    //PETScWrappers::MPI::Vector tmp(locally_owned_dofs_dam, mpi_communicator);
    /*MatrixTools::apply_boundary_values(
      boundary_values, system_matrix_dam, tmp, system_rhs_dam, false);*/
    //d_dam = tmp;
    
  }


// Solve one iteration

  template <int dim, int spacedim>
  void PFC<dim, spacedim>::solve_NR_disp()
  {
    
    
    
    //pcout << " norm of residue u is " << system_rhs_u.l2_norm() << std::endl;

    const unsigned int n_iterations = solve_linear_problem_disp();

    //pcout << "    Solver converged in " << n_iterations << " iterations."
    //      << std::endl;
    disp_NR.add(1.0, d_disp);

    //pcout << "    Updating quadrature point data..." << std::flush;
    //update_quadrature_point_history(); // Comes at end of NR
    //pcout << std::endl;
      
  }

  template <int dim, int spacedim>
  void PFC<dim, spacedim>::solve_NR_dam()
  {
    
    
    //pcout << " norm of residue s is " << system_rhs_s.l2_norm() << std::endl;

    const unsigned int n_iterations = solve_linear_problem_dam();

    //pcout << "    Solver converged in " << n_iterations << " iterations."
    //      << std::endl;
    dam_NR.add(1.0, d_dam);
    
    #ifdef PFC_DEBUG
    //outfile << "\n\tDam NR: \n"<< dam_NR<< "\n\tDam ds: \n"<<d_dam;
    #endif
    
    for (unsigned int i = 0; i< dof_handler_dam.n_dofs(); i++)
    {
      dam_NR[i] = ((dam_NR[i] < 0.) ? 0. : dam_NR[i]);
      dam_NR[i] = ((dam_NR[i] > dam[i]) ? dam[i] : dam_NR[i]);
    }
    //pcout << "    Updating quadrature point data..." << std::flush;
    //update_quadrature_point_history(); // Comes at end of NR
    //pcout << std::endl;
    
  }



  // Linear problem solver for disp
  template <int dim, int spacedim>
  unsigned int PFC<dim, spacedim>::solve_linear_problem_disp()
  {
   
    PETScWrappers::MPI::Vector distributed_incremental_displacement(
      locally_owned_dofs_disp, mpi_communicator);
    distributed_incremental_displacement = d_disp;

    SolverControl solver_control(dof_handler_disp.n_dofs(),
                                 lin_solve_disp_tol * system_rhs_disp.l2_norm());
    
    #ifdef PFC_DEBUG
      //solver_control.enable_history_data();
    #endif

    //PETScWrappers::SolverCG cg(solver_control);

    //PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix_disp);

    //cg.solve(system_matrix_disp,
     //        distributed_incremental_displacement,
      //       system_rhs_disp,
      //       preconditioner);
    
    
    // Create the LU solver
    PETScWrappers::SparseDirectMUMPS solver(solver_control);

    // Solve the system using LU
    solver.solve(system_matrix_disp, distributed_incremental_displacement, system_rhs_disp);
    
    #ifdef PFC_DEBUG
    /*
    std::vector< double > history_data = solver_control.get_history_data();
    
    outfile << "\n CG step Residues for Displacement:\n";
    
    std::vector<double>::iterator res;
    for (res=history_data.begin();res!=history_data.end();res++)
      outfile << *res << "\t";
    */
    #endif
    
    d_disp = distributed_incremental_displacement;

    hanging_node_constraints_disp.distribute(d_disp);

    return solver_control.last_step();
       
  }
  
  // for dam
  template <int dim, int spacedim>
  unsigned int PFC<dim, spacedim>::solve_linear_problem_dam()
  {
    
    PETScWrappers::MPI::Vector distributed_incremental_damage(
      locally_owned_dofs_dam, mpi_communicator);
    distributed_incremental_damage = d_dam;

    SolverControl solver_control(dof_handler_dam.n_dofs(),
                                 lin_solve_dam_tol * system_rhs_dam.l2_norm());
    
    
    #ifdef PFC_DEBUG
      solver_control.enable_history_data();
    #endif
    
    /*PETScWrappers::SolverCG cg(solver_control);
    //SolverBiCG

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix_dam);
    //PETScWrappers::PreconditionGMRES preconditioner(system_matrix_dam);
    
    cg.solve(system_matrix_dam,
             distributed_incremental_damage,
             system_rhs_dam,
             preconditioner);*/
    
    // Create the LU solver
    PETScWrappers::SparseDirectMUMPS solver(solver_control);

    // Solve the system using LU
    solver.solve(system_matrix_dam, distributed_incremental_damage, system_rhs_dam);
    
    #ifdef PFC_DEBUG
    std::vector< double > history_data = solver_control.get_history_data();
    
    outfile << "\n CG step Residues for Damage:\n";

    std::vector<double>::iterator res;
    for (res=history_data.begin();res!=history_data.end();res++)
      outfile << *res << "\t";
    #endif
    
    d_dam = distributed_incremental_damage;

    hanging_node_constraints_dam.distribute(d_dam);

    return solver_control.last_step();
     
  }
  
  
  // Update Qp history
  
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::update_Qp_history_disp_NR() // Update after disp_NR  convergence: eps_e, eps_p, strss, alpha, Psi_plus_max
  {

    // For strains and all
    FEValues<dim> fe_values(fe_disp,
                  quadrature_formula,
                  update_values | update_gradients |
                  update_quadrature_points | update_JxW_values);
    
    std::vector<std::vector<Tensor<1, dim>>> disp_grad(
      quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim));
    
    for (auto &cell : dof_handler_disp.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        
        //outfile<<"\n \n Printing Displacement NR state variables at cell:\t"<<cell->id(); 

        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());

        fe_values.reinit(cell);
        fe_values.get_function_gradients(disp_NR, disp_grad);
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          //outfile << "\n\t Quadrature point:\t" << q;
          // compute strain from u
          // compute delta strain from eps_e 0 and eps_p 0 
          const Tensor<2, dim> tmp_strain = get_strain(disp_grad[q]);
          Tensor<2, spacedim> strain;
          //std::cout<<tmp_strain<<std::endl;
          
            strain = tensordimtospacedim<dim, spacedim> ( tmp_strain);
          
          //pcout << del_strain << std::endl;
          
          //std::cout<<"Cijkl"<<tmp_Cijkl<<std::endl;
          //std::cout<<"Strain"<<local_quadrature_points_history[q].eps_e[1]<<std::endl;
          //std::cout<<"Stress"<<local_quadrature_points_history[q].strss[1]<<std::endl;
          /// NOW, testing predictor corrector
          double yield_fn;
          
          Tensor<2, spacedim> strss_pc, eps_e_pc, eps_p_pc;
          Tensor<4, spacedim> Cijkl_pc;
          double alpha_pc, del_lambda_pc;
          
          if ( !enable_brittle_LEFM )
          {
            predictor_corrector ( strain,
                                  local_quadrature_points_history[q].eps_e,
                                  local_quadrature_points_history[q].eps_p,
                                  local_quadrature_points_history[q].alpha,
                                  local_quadrature_points_history[q].g_hat,
                                  yield_fn,
                                  del_lambda_pc,
                                  alpha_pc,
                                  strss_pc,
                                  eps_e_pc,
                                  eps_p_pc,
                                  Cijkl_pc);
            
            local_quadrature_points_history[q].del_lambda = del_lambda_pc;
            local_quadrature_points_history[q].alpha_NR = alpha_pc;
            local_quadrature_points_history[q].strss_NR = strss_pc;
            local_quadrature_points_history[q].eps_e_NR = eps_e_pc;
            local_quadrature_points_history[q].eps_p_NR = eps_p_pc;
          }
          
          if ( enable_brittle_LEFM )
          {
            eps_e_pc = strain;
            Cijkl_pc = calc_strss_strain_tensor_constant(
                        eps_e_pc,
                        local_quadrature_points_history[q].g_hat);
            strss_pc = calc_strss_from_strain (Cijkl_pc, eps_e_pc);
            
            local_quadrature_points_history[q].strss_NR = strss_pc;
            local_quadrature_points_history[q].eps_e_NR = eps_e_pc;
          }
          //local_quadrature_points_history[q].Psi_plus_max_NR = Psi_plus(eps_e_pc);
          
          local_quadrature_points_history[q].Psi_plus_max_NR = 
            std::max(local_quadrature_points_history[q].Psi_plus_max,
                     Psi_plus(eps_e_pc));
          
          //outfile<<"\n\t eps_e:\t"<<eps_e_pc; 
          //outfile<<"\n\t eps_p:\t"<<eps_p_pc; 
          //outfile<<"\n\t alpha:\t"<<alpha_pc; 
          //outfile<<"\n\t Psi_plus_max:\t"<<Psi_plus(eps_e_pc); 
 /*         
          double tmp_Yd = yield_function(strss_pc,alpha_pc);
 */         
          
          
             
           //plastic_strain << ((tmp_strain.norm() <= eps_p_pc.norm()) ? 1 : 0) << "\t";
        }
      }
    
  }
  
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::update_Qp_history_dam_NR() // Update after dam_NR convergence: acc_eps_p, g_hat
  {
    
    
      FEValues<dim> fe_values(fe_dam,
                    quadrature_formula,
                    update_values | update_gradients |
                    update_quadrature_points | update_JxW_values);
      
      std::vector<double> val_dam(quadrature_formula.size());
      
      for (auto &cell : dof_handler_dam.active_cell_iterators())
        if (cell->is_locally_owned())
        {
          
          //outfile<<"\n \n Printing Damage NR state variables at cell:"<<cell->id()<<std::endl; 

           PointHistory<spacedim> *local_quadrature_points_history =
            reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
          Assert(local_quadrature_points_history >=
                   &quadrature_point_history.front(),
                 ExcInternalError());
          Assert(local_quadrature_points_history <=
                   &quadrature_point_history.back(),
                 ExcInternalError());


          fe_values.reinit(cell);
          fe_values.get_function_values(dam_NR, val_dam);
          for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
          {

            //outfile << "\n\t Quadrature point:\t" << q ;

      // Update p before g_hat after NR s
            //std::cout<< local_quadrature_points_history[q].p<<std::endl;
            local_quadrature_points_history[q].acc_eps_p_NR = 
              local_quadrature_points_history[q].acc_eps_p + 
                
                (local_quadrature_points_history[q].alpha_NR -
                  local_quadrature_points_history[q].alpha);
                      
            if (enable_elastoplastic_brittle_fracture || enable_brittle_LEFM)
              local_quadrature_points_history[q].acc_eps_p_NR = Ep_crit;
            // Update g_hat before, during and after NR s
            local_quadrature_points_history[q].g_hat =
              std::pow(val_dam[q], 2 * local_quadrature_points_history[q].acc_eps_p_NR / Ep_crit) + Eta;

            //outfile<<"\n\t acc_eps_p_NR:\t"<<local_quadrature_points_history[q].acc_eps_p_NR; 
            //outfile<<"\n\t g_hat:\t"<<local_quadrature_points_history[q].g_hat; 

          }
        }
    
  }
  
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::accept_NR_solution() // Accept the solution after Energy convergence
  {
    
    disp = disp_NR;
    dam = dam_NR;
  }
  
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::accept_NR_state() // Accept the State vars after energy convergence
  {
    
    for (auto &cell : dof_handler_dam.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());


        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          local_quadrature_points_history[q].strss =
            local_quadrature_points_history[q].strss_NR;
          
          local_quadrature_points_history[q].eps_e =
            local_quadrature_points_history[q].eps_e_NR;
          
          local_quadrature_points_history[q].eps_p =
            local_quadrature_points_history[q].eps_p_NR;
                        
          local_quadrature_points_history[q].alpha =
            local_quadrature_points_history[q].alpha_NR;
          
          local_quadrature_points_history[q].Psi_plus_max = 
            local_quadrature_points_history[q].Psi_plus_max_NR;
          
          local_quadrature_points_history[q].acc_eps_p = 
            local_quadrature_points_history[q].acc_eps_p_NR;
        }
      }
  }
  
  // Temporary, may or may not be used
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::reset_vark() // reset all k variables except g_hat and Psi_plus_max
  {
    
    disp_NR = disp;
    dam_NR = dam;
    
    for (auto &cell : dof_handler_dam.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());


        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
 /*         local_quadrature_points_history[q].strss =
            local_quadrature_points_history[q].strss_NR;
          
          local_quadrature_points_history[q].eps_e =
            local_quadrature_points_history[q].eps_e_NR;
          
          local_quadrature_points_history[q].eps_p =
            local_quadrature_points_history[q].eps_p_NR;
                        
          local_quadrature_points_history[q].alpha =
            local_quadrature_points_history[q].alpha_NR;
          
          local_quadrature_points_history[q].Psi_plus_max = 
            local_quadrature_points_history[q].Psi_plus_max_NR;
          
          local_quadrature_points_history[q].p = 
            local_quadrature_points_history[q].p_k;*/
        }
      }
    
  }
    
  // Output results
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::output_results(const unsigned int cycle)
  {
    //DataOut<dim> data_out_disp, data_out_dam;
    DataOut<dim> data_out;
    
    data_out.attach_dof_handler(dof_handler_disp);
    
    if(cycle % outputFrequency ==0)
    {
      /*
      std::vector<std::string> solution_names;
      
      switch (dim)
        {
          case 1:
            solution_names.emplace_back("delta_x");
            break;
          case 2:
            solution_names.emplace_back("delta_x");
            solution_names.emplace_back("delta_y");
            break;
          case 3:
            solution_names.emplace_back("delta_x");
            solution_names.emplace_back("delta_y");
            solution_names.emplace_back("delta_z");
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
      data_out_u.add_data_vector(u, solution_names);
      */
      // Does it interpret? Yes it does. STEP 22
      std::vector<std::string> disp_name(dim, "Displacement");
      
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
          
      data_out.add_data_vector(disp,
                               disp_name,
                               DataOut<dim>::type_dof_data,
                               data_component_interpretation);
      /*data_out_u.add_data_vector(u,
                               "Displacement",
                               DataOut<dim>::type_dof_data);*/// does the same thing
      
      Vector<double> norm_of_stress(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              Tensor<2, spacedim> accumulated_stress;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_stress +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .strss;

              norm_of_stress(cell->active_cell_index()) =
                (accumulated_stress / quadrature_formula.size()).norm();
            }
          else
            norm_of_stress(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(norm_of_stress,
                                    "Norm_of_stress: $||\\mathbf{\\sigma}||$");
      
      Vector<double> norm_of_ep(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              Tensor<2, spacedim> accumulated_ep;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_ep +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .eps_p;

              norm_of_ep(cell->active_cell_index()) =
                (accumulated_ep / quadrature_formula.size()).norm();
            }
          else
            norm_of_ep(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(norm_of_ep,
                                    "Norm_of_eps_p: $||\\mathbf{\\varepsilon}^{p}||$");
      
      Vector<double> norm_of_ee(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              Tensor<2, spacedim> accumulated_ee;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_ee +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .eps_e;

              norm_of_ee(cell->active_cell_index()) =
                (accumulated_ee / quadrature_formula.size()).norm();
            }
          else
            norm_of_ee(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(norm_of_ee,
                                    "Norm_of_eps_e: $||\\mathbf{\\varepsilon}^{e}||$");
      
      Vector<double> History(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              double accumulated_his=0;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_his +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .Psi_plus_max;

              History(cell->active_cell_index()) =
                (accumulated_his / quadrature_formula.size());
            }
          else
            History(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(History,
                                    "Max_Elastic_Energy_Density_History: $\\mathcal{H}_{e}$");
  /*    
      Vector<double> Tracey(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              double accumulated_tr=0;//here is the problem
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_tr +=trace(
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .eps_e);

              Tracey(cell->active_cell_index()) =
                (accumulated_tr / quadrature_formula.size());
            }
          else
            Tracey(cell->active_cell_index()) = -1e+20;
      }
      data_out_u.add_data_vector(Tracey, "Tracey");
  */  
      
      Vector<double> Fn(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              double Fn_yield=0;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                Fn_yield += yield_function(
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .strss
                    ,reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .alpha );

              Fn(cell->active_cell_index()) =
                (Fn_yield / quadrature_formula.size());
            }
          else
            Fn(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(Fn, "Yield_fn: $f(\\mathbf{\\sigma}_{dev},\\alpha})$");
      
      Vector<double> alpha(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              double accumulated_alpha=0;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_alpha +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .alpha;

              alpha(cell->active_cell_index()) =
                (accumulated_alpha / quadrature_formula.size());
            }
          else
            alpha(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(alpha, "Alpha: $\\alpha$");
      
      Vector<double> accumulated_plasticity(triangulation.n_active_cells());
      {
        for (auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              double accumulated_plastic_strain_var=0;
              for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                accumulated_plastic_strain_var +=
                  reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer())[q]
                    .acc_eps_p;

              accumulated_plasticity(cell->active_cell_index()) =
                (accumulated_plastic_strain_var / quadrature_formula.size());
            }
          else
            accumulated_plasticity(cell->active_cell_index()) = -1e+20;
      }
      data_out.add_data_vector(accumulated_plasticity,
                                  "Accumulated_plasticity: $||\\varepsilon^{p}_{eq}||$");
      
      
      data_out.add_data_vector(dof_handler_dam, dam, "Damage: $s$");

  // Partitioning ?Where is the problem /// This is fine its the tensor part
      
      std::vector<types::subdomain_id> partition_int(
        triangulation.n_active_cells());
      GridTools::get_subdomain_association(triangulation, partition_int);
      const Vector<double> partitioning(partition_int.begin(),
                                        partition_int.end());
      // u
      data_out.add_data_vector(partitioning, "partitioning");
      data_out.build_patches();

      //data_out_dam.add_data_vector(partitioning, "partitioning");
      //data_out_dam.build_patches();
      
      const std::string pvtu_filename_disp = data_out.write_vtu_with_pvtu_record(
        "./results/", "solution", timestep_no, mpi_communicator, 4);

      if (this_mpi_process == 0)
        {
          static std::vector<std::pair<double, std::string>> times_and_names;
          times_and_names.emplace_back(time, pvtu_filename_disp);
          std::ofstream pvd_output("results/solution.pvd");
          DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }
      // s
      

     /* 

      const std::string pvtu_filename_dam = data_out_dam.write_vtu_with_pvtu_record(
        "./results/", "solution_s", timestep_no, mpi_communicator, 4);

      if (this_mpi_process == 0)
        {
          static std::vector<std::pair<double, std::string>> times_and_names;
          times_and_names.emplace_back(time, pvtu_filename_dam);
          std::ofstream pvd_output("results/solution_s.pvd");
          DataOutBase::write_pvd_record(pvd_output, times_and_names);
        }*/
    }
        
  }
  
  //Q point output
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::output_Qp_results(const unsigned int cycle)
  {
    ///////////////////////////////////////// Strss field Not working properly
    

    FEValues<dim> fe_values(fe_disp,
                    quadrature_formula,
                    update_values | update_gradients |
                      update_quadrature_points | update_JxW_values);

    //const unsigned int dofs_per_cell = fe_disp.n_dofs_per_cell();
    //const unsigned int n_q_points    = quadrature_formula.size();
    
    FE_Q<dim>     tensor_fe (1);
    DoFHandler<dim> tensor_dof_handler (triangulation);
    tensor_dof_handler.distribute_dofs (tensor_fe);

    std::vector< std::vector< Vector<double> > >
                 strs_field (dim, std::vector< Vector<double> >(dim)),
                 local_strs_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_strs_cell_values (dim, std::vector< Vector<double> >(dim)); 
    // strss ZZ field separate
    Vector<double> strsZZ_field, local_strsZZ_values_at_qpoints, local_strsZZ_cell_values;

    std::vector< std::vector< Vector<double> > >
                 epse_field (dim, std::vector< Vector<double> >(dim)),
                 local_epse_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_epse_cell_values (dim, std::vector< Vector<double> >(dim)); 
   
    Vector<double> epseZZ_field, local_epseZZ_values_at_qpoints, local_epseZZ_cell_values;
    
    std::vector< std::vector< Vector<double> > >
                 epsp_field (dim, std::vector< Vector<double> >(dim)),
                 local_epsp_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_epsp_cell_values (dim, std::vector< Vector<double> >(dim)); 
    
    Vector<double> epspZZ_field, local_epspZZ_values_at_qpoints, local_epspZZ_cell_values;
    
    Vector<double> alpha_field, local_alpha_values_at_qpoints, local_alpha_cell_values;
    Vector<double> his_field, local_his_values_at_qpoints, local_his_cell_values;
    Vector<double> p_field, local_p_values_at_qpoints, local_p_cell_values;
    Vector<double> g_field, local_g_values_at_qpoints, local_g_cell_values;
    
    if(cycle % outputFrequency ==0)
    {
      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
        {
          epse_field[i][j].reinit(tensor_dof_handler.n_dofs());
          local_epse_values_at_qpoints[i][j].reinit(quadrature_formula.size());
          local_epse_cell_values[i][j].reinit(tensor_fe.dofs_per_cell);
          
          epsp_field[i][j].reinit(tensor_dof_handler.n_dofs());
          local_epsp_values_at_qpoints[i][j].reinit(quadrature_formula.size());
          local_epsp_cell_values[i][j].reinit(tensor_fe.dofs_per_cell);
          
          strs_field[i][j].reinit(tensor_dof_handler.n_dofs());
          local_strs_values_at_qpoints[i][j].reinit(quadrature_formula.size());
          local_strs_cell_values[i][j].reinit(tensor_fe.dofs_per_cell);

        }
      
      strsZZ_field.reinit(tensor_dof_handler.n_dofs()); 
      local_strsZZ_values_at_qpoints.reinit(quadrature_formula.size());
      local_strsZZ_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      epseZZ_field.reinit(tensor_dof_handler.n_dofs()); 
      local_epseZZ_values_at_qpoints.reinit(quadrature_formula.size());
      local_epseZZ_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      epspZZ_field.reinit(tensor_dof_handler.n_dofs()); 
      local_epspZZ_values_at_qpoints.reinit(quadrature_formula.size());
      local_epspZZ_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      alpha_field.reinit(tensor_dof_handler.n_dofs()); 
      local_alpha_values_at_qpoints.reinit(quadrature_formula.size());
      local_alpha_cell_values.reinit(tensor_fe.dofs_per_cell);
        
      his_field.reinit(tensor_dof_handler.n_dofs()); 
      local_his_values_at_qpoints.reinit(quadrature_formula.size());
      local_his_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      p_field.reinit(tensor_dof_handler.n_dofs()); 
      local_p_values_at_qpoints.reinit(quadrature_formula.size());
      local_p_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      g_field.reinit(tensor_dof_handler.n_dofs()); 
      local_g_values_at_qpoints.reinit(quadrature_formula.size());
      local_g_cell_values.reinit(tensor_fe.dofs_per_cell);
      
      FullMatrix<double> qpoint_to_dof_matrix (tensor_fe.dofs_per_cell,
                                               quadrature_formula.size());
      FETools::compute_projection_from_quadrature_points_matrix
                (tensor_fe,
                 quadrature_formula, quadrature_formula,
                 qpoint_to_dof_matrix);
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_disp.begin_active(),
                                                     endc = dof_handler_disp.end(),
                                                     dg_cell =
                                                      tensor_dof_handler.begin_active();
      for (; cell!=endc; ++cell, ++dg_cell)
        if(cell->is_locally_owned())
        {
          
          fe_values.reinit(cell);
          PointHistory<spacedim> *local_quadrature_points_history
            = reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
          Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
                  ExcInternalError());
          Assert (local_quadrature_points_history < &quadrature_point_history.back(),
                  ExcInternalError());
          
          fe_values.reinit(cell);
          //fe_values.get_function_gradients(solution,
          //                             solution_grad);

          // Then loop over the quadrature points of this cell:
          for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
            {

          
              for (unsigned int i=0; i<dim; i++)
                for (unsigned int j=0; j<dim; j++)
                {

                  local_epse_values_at_qpoints[i][j](q)
                      = local_quadrature_points_history[q].eps_e[i][j];
                  local_epsp_values_at_qpoints[i][j](q)
                      = local_quadrature_points_history[q].eps_p[i][j];
                  local_strs_values_at_qpoints[i][j](q)
                      = local_quadrature_points_history[q].strss[i][j]; 

                }
              
              local_strsZZ_values_at_qpoints(q)=local_quadrature_points_history[q].strss[2][2];
              local_epseZZ_values_at_qpoints(q)=local_quadrature_points_history[q].eps_e[2][2];
              local_epspZZ_values_at_qpoints(q)=local_quadrature_points_history[q].eps_p[2][2];
              
              local_alpha_values_at_qpoints(q)=local_quadrature_points_history[q].alpha;
              local_his_values_at_qpoints(q)=local_quadrature_points_history[q].Psi_plus_max;
              local_p_values_at_qpoints(q)=local_quadrature_points_history[q].acc_eps_p;
              local_g_values_at_qpoints(q)=local_quadrature_points_history[q].g_hat;

            }

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
            {
              for (unsigned int i=0; i<dim; i++)
                for (unsigned int j=0; j<dim; j++)
                  {

                    qpoint_to_dof_matrix.vmult (local_epse_cell_values[i][j],
                                                local_epse_values_at_qpoints[i][j]);

                    qpoint_to_dof_matrix.vmult (local_epsp_cell_values[i][j],
                                                local_epsp_values_at_qpoints[i][j]);

                    qpoint_to_dof_matrix.vmult (local_strs_cell_values[i][j],
                                                local_strs_values_at_qpoints[i][j]);

                    dg_cell->set_dof_values (local_epse_cell_values[i][j],
                                             epse_field[i][j]);
                                             
                    dg_cell->set_dof_values (local_epsp_cell_values[i][j],
                                             epsp_field[i][j]);
                                             
                    dg_cell->set_dof_values (local_strs_cell_values[i][j],
                                             strs_field[i][j]);

                  }
                  
              qpoint_to_dof_matrix.vmult (local_strsZZ_cell_values,
                                                local_strsZZ_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_epseZZ_cell_values,
                                                local_epseZZ_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_epspZZ_cell_values,
                                                local_epspZZ_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_alpha_cell_values,
                                                local_alpha_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_his_cell_values,
                                                local_his_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_p_cell_values,
                                                local_p_values_at_qpoints);
              
              qpoint_to_dof_matrix.vmult (local_g_cell_values,
                                                local_g_values_at_qpoints);
              
              
              dg_cell->set_dof_values(local_strsZZ_cell_values,
                                            strsZZ_field);
              
              dg_cell->set_dof_values(local_epseZZ_cell_values,
                                            epseZZ_field);
              
              dg_cell->set_dof_values(local_epspZZ_cell_values,
                                            epspZZ_field);
              
              dg_cell->set_dof_values(local_alpha_cell_values,
                                            alpha_field);
              
              dg_cell->set_dof_values(local_his_cell_values,
                                            his_field);
              
              dg_cell->set_dof_values(local_p_cell_values,
                                            p_field);
              
              dg_cell->set_dof_values(local_g_cell_values,
                                            g_field);

            }
        }
        

      //
      DataOut<dim> data_out_tensor; 

      //
      data_out_tensor.attach_dof_handler(tensor_dof_handler);  /// Fix this
//std::cout <<tensor_dof_handler.n_dofs();
      //

      data_out_tensor.add_data_vector(epse_field[0][0], "epse11");
      data_out_tensor.add_data_vector(epse_field[0][1], "epse12");
      data_out_tensor.add_data_vector(epse_field[1][0], "epse21");
      data_out_tensor.add_data_vector(epse_field[1][1], "epse22");
      data_out_tensor.add_data_vector(epseZZ_field, "epse33");
        
      data_out_tensor.add_data_vector(epsp_field[0][0], "epsp11");
      data_out_tensor.add_data_vector(epsp_field[0][1], "epsp12");
      data_out_tensor.add_data_vector(epsp_field[1][0], "epsp21");
      data_out_tensor.add_data_vector(epsp_field[1][1], "epsp22");
      data_out_tensor.add_data_vector(epspZZ_field, "epsp33");
        
      data_out_tensor.add_data_vector(strs_field[0][0], "strss11");
      data_out_tensor.add_data_vector(strs_field[0][1], "strss12");
      data_out_tensor.add_data_vector(strs_field[1][0], "strss21");
      data_out_tensor.add_data_vector(strs_field[1][1], "strss22");
      data_out_tensor.add_data_vector(strsZZ_field, "strss33");

    
      data_out_tensor.add_data_vector(alpha_field, "Alpha");
      data_out_tensor.add_data_vector(his_field, "Max_Elastic_Energy_Density_History");
      data_out_tensor.add_data_vector(p_field, "Accumulated_plasticity");
      data_out_tensor.add_data_vector(g_field, "g_hat");

      
      //
      data_out_tensor.build_patches(); 
      
      const std::string pvtu_filename_tensor = data_out_tensor.write_vtu_with_pvtu_record(
      "./results/", "tensor", timestep_no, mpi_communicator, 4);
      
      if (this_mpi_process == 0)
      {
        static std::vector<std::pair<double, std::string>> times_and_names_tensor;
        times_and_names_tensor.emplace_back(time, pvtu_filename_tensor);
        std::ofstream pvd_output_tensor("results/tensor.pvd");
        DataOutBase::write_pvd_record(pvd_output_tensor, times_and_names_tensor);
      }

    }  
  }
  
  // Stress strain
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::output_stress_strain_plot()
  {
    
    FEValues<dim> fe_values(fe_disp,
                            quadrature_formula,
                            update_gradients |
                            update_quadrature_points | update_JxW_values);
    
    const unsigned int dofs_per_cell = fe_disp.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
    unsigned int total_cells(0), total_dofs(0);
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    
    std::vector<bool> boundary_dofs_y(dof_handler_disp.n_dofs()),
                      boundary_dofs_x(dof_handler_disp.n_dofs()),
                      boundary_dofs_blk(dof_handler_disp.n_dofs());// Depends on BC
    
    FEValuesExtractors::Scalar                x_component(0);
    FEValuesExtractors::Scalar                y_component(1);
    
    DoFTools::extract_boundary_dofs(dof_handler_disp,
                                    fe_disp.component_mask(y_component),//ComponentMask(),
                                    boundary_dofs_y,
                                    {11});
    
    DoFTools::extract_boundary_dofs(dof_handler_disp,
                                    fe_disp.component_mask(x_component),//ComponentMask(),
                                    boundary_dofs_x,
                                    {11});
    
    /*DoFTools::extract_boundary_dofs(dof_handler_disp,
                                    fe_disp.component_mask(y_component),//ComponentMask(),
                                    boundary_dofs_blk,
                                    {11});*/
    DoFTools::extract_boundary_dofs(dof_handler_disp,
                                    fe_disp.component_mask(x_component),//ComponentMask(),
                                    boundary_dofs_blk,
                                    {12});
    
    //std::cout << boundary_dofs_blk.size() <<std::endl;
    //for (unsigned int i = 0; i < boundary_dofs.size(); ++i)
    //  std::cout << i << "\t"<< boundary_dofs[i] <<std::endl;
    
    const IndexSet bd = DoFTools::extract_boundary_dofs(dof_handler_disp,
                                                        fe_disp.component_mask(y_component),
                                                        {11});
    
    Vector<double>     cell_residual(dofs_per_cell);
    system_rxn = 0; 

    for (auto &cell : dof_handler_disp.active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
    {
      PointHistory<spacedim> *local_quadrature_points_history =
            reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());

      
      cell_residual    = 0;

      fe_values.reinit(cell);
          
      for (unsigned int p = 0; p < dofs_per_cell; ++p)
      {
        const unsigned int component_p =
          fe_disp.system_to_component_index(p).first;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {             
               
              const Tensor<2, spacedim> &stress =
                local_quadrature_points_history[q_point].strss;
              for(int J=0;J<spacedim;J++)                               
              {
                 cell_residual(p) += fe_values.shape_grad(p, q_point)[J] * 
                                     stress[component_p][J] *
                                     fe_values.JxW(q_point); 
              }
              
          }
      }

      //std::cout<<"\n cell_residual:"<<cell_residual; 

      cell->get_dof_indices(local_dof_indices);
      
      hanging_node_constraints_disp.distribute_local_to_global(cell_residual,
                                                      local_dof_indices,
                                                      system_rxn);

            
    }

    system_rxn.compress(VectorOperation::add);

    
    //Reaction forces

    double rxn_y(0.), rxn_x(0.), rxn_blk(0.);
    double tensile_strain(0.), shear_strain(0.) , bulk_strain(0.),
           disp_y(0.)        , disp_x(0.)       , disp_xny(0.)   ;
    int    disp_count_y(0)   , disp_count_x(0)  , disp_count_xny(0) ;
    double tensile_stress(0.), shear_stress(0.) , bulk_stress(0.);
    double avg_eps_p(0.), avg_eps_p_eq(0.), avg_g_hat(0.), avg_damage(0.);
        
    MappingQ<dim, dim> mapping(1);

    FEValues<dim , dim> fe_v(mapping,
                            fe_disp,
                            quadrature_formula,
                            update_values |
                              update_quadrature_points | update_JxW_values);

    std::vector<Point<dim>> support_points(dof_handler_disp.n_dofs());
    DoFTools::map_dofs_to_support_points<dim, dim>(mapping,
                                                   dof_handler_disp,
                                                   support_points);
    
    for (auto &cell : dof_handler_disp.active_cell_iterators())
    if (cell->is_locally_owned())
    {
      PointHistory<dim> *local_quadrature_points_history =
            reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
      
      total_cells++;
      
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        avg_eps_p     += (local_quadrature_points_history[q_point].eps_p).norm();
        avg_eps_p_eq  += local_quadrature_points_history[q_point].acc_eps_p;
        avg_g_hat     += local_quadrature_points_history[q_point].g_hat;
      }
      
    }
    
    avg_damage= dam.l1_norm();
    total_dofs = dof_handler_dam.n_dofs();
    
    for (unsigned int i = 0; i < dof_handler_disp.n_dofs(); ++i)
      if(locally_owned_dofs_disp.is_element(i))
        {
          //if (boundary_dofs_x[i] == true)
          //{
          //  rxn_x+=system_rxn(i);
          //  disp_x+=disp(i); 
          //  disp_count_x++;
          //}
          
          if (boundary_dofs_y[i] == true)
          {
            rxn_y+=system_rxn(i);
            disp_y+=disp(i); 
            disp_count_y++;
            
            //rxn_blk+=system_rxn(i);// for bulk stress calculation
            
          }
          
          //if (boundary_dofs_blk[i] == true)
          //{
          //  rxn_blk+=system_rxn(i);
          //  //disp_xny+= (disp(i) * sample_size); // + dx * dy missing // not correct
          //  //disp_count_xny++;
          //}
          
        }

    //rxn_x= Utilities::MPI::sum(rxn_x, mpi_communicator);
    //disp_x= Utilities::MPI::sum(disp_x, mpi_communicator);
    //disp_count_x = Utilities::MPI::sum(disp_count_x, mpi_communicator);
    
    rxn_y= Utilities::MPI::sum(rxn_y, mpi_communicator);
    disp_y= Utilities::MPI::sum(disp_y, mpi_communicator); 
    disp_count_y = Utilities::MPI::sum(disp_count_y, mpi_communicator);
    
    // Recheck
    //rxn_blk= Utilities::MPI::sum(rxn_blk, mpi_communicator);
    //disp_xny= Utilities::MPI::sum(disp_xny, mpi_communicator); 
    //disp_count_xny = Utilities::MPI::sum(disp_count_xny, mpi_communicator);
    
    avg_eps_p= Utilities::MPI::sum(avg_eps_p, mpi_communicator);
    avg_eps_p_eq= Utilities::MPI::sum(avg_eps_p_eq, mpi_communicator);
    avg_g_hat= Utilities::MPI::sum(avg_g_hat, mpi_communicator);
    avg_damage= Utilities::MPI::sum(avg_damage, mpi_communicator);
    total_cells= Utilities::MPI::sum(total_cells, mpi_communicator);
    total_dofs= Utilities::MPI::sum(total_dofs, mpi_communicator);
    
    avg_eps_p= avg_eps_p/ (total_cells * n_q_points);
    avg_eps_p_eq= avg_eps_p_eq/ (total_cells * n_q_points);
    avg_g_hat= avg_g_hat/ (total_cells * n_q_points);
    avg_damage= avg_damage/ total_dofs;
    
    // Strains
    tensile_strain = disp_y/disp_count_y/size_y;
    //shear_strain   = disp_x/disp_count_x/sample_size;
    //bulk_strain    = disp_xny/disp_count_xny/(sample_size * sample_size); // not correct
    //bulk_strain    = (disp_y * sample_size)
    //                 / disp_count_y 
    //                 / (sample_size * sample_size); // temporary
    
    // Stresses
    // 3D : rxn/pow(sample_size,2);
    tensile_stress = rxn_y/(size_x * size_z); 
    //shear_stress   = rxn_x/(size_x * size_z);
    //bulk_stress    = rxn_blk/(sample_size + sample_size);
    
    if(this_mpi_process == 0)
      ssFile<< std::setprecision(8)
            << time << "\t"
            << (disp_y/disp_count_y) << "\t"
            << rxn_y/(1000 * size_z) << "\t" //(1./ size_z) for scaling the forces per unit lgth
            << avg_eps_p << "\t"
            << avg_eps_p_eq << "\t"
            << avg_eps_p_eq/Ep_crit << "\t"
            << avg_g_hat << "\t"
            << avg_damage << "\t"
            << tensile_strain << "\t"
            << tensile_stress << std::endl;
            //<< shear_strain << "\t"
            //<< shear_stress << "\t"
            //<< bulk_strain << "\t"
            //<< bulk_stress << std::endl; 
        
  }
  
  // Calculate timestep //use before update
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::calc_timestep()
  {
    double deps_p_max(0.), d_dam_max(0.);
    double del = 0.;
    
    //Vector<double> dam_diff; //damage difference
    //dam_diff.reinit(dof_handler_dam.n_dofs()); 
    
    for (auto &cell : dof_handler_disp.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());
      
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          //del = (local_quadrature_points_history[q].eps_p_NR -
          //       local_quadrature_points_history[q].eps_p).norm();
          
          //the scaled one?
          //del = ((local_quadrature_points_history[q].eps_p_NR -
          //       local_quadrature_points_history[q].eps_p).norm() /
          //       (local_quadrature_points_history[q].eps_e_NR -
          //       local_quadrature_points_history[q].eps_e).norm());
          
          del = (local_quadrature_points_history[q].eps_p_NR -
                 local_quadrature_points_history[q].eps_p).norm();
                 // /
                 //(local_quadrature_points_history[q].eps_e_NR -
                 //local_quadrature_points_history[q].eps_e).norm());
          
          
          deps_p_max = std::max(deps_p_max, del);
        }
      
      }
      
    deps_p_max= Utilities::MPI::max(deps_p_max, mpi_communicator);
    
    //check if max change in eps_p is > 0.002. If yes, reduce dt
    
    dt = default_dt * ((deps_p_max > d_eps_p_limit)? (d_eps_p_limit / deps_p_max) : 1.);
    
     // reduce time step after disp=5 -- expected failure !!
    //if(!is_reduced_dt && present_strain>=8.4)  // LOAD_TYPE change
    //{
    //  dt = dt/5.0; 
    //  is_reduced_dt = true; 
    //}
     
    /*
    dam_diff.add(1., dam, -1., dam_NR);
    d_dam_max = dam_diff.linfty_norm();
    d_dam_max= Utilities::MPI::max(d_dam_max, mpi_communicator);
    //
    
    dt = dt * ((d_dam_max > d_dam_max_limit)? (d_dam_max_limit / d_dam_max) : 1.); 
     */
    //pcout << dt<<std::endl;
  }
  
  // Compute Psi total and then energy  
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::compute_Energy()
  {
    
    double Energy, Psi_tot;
    Psi_tot = 0;
    Energy = 0;
    
    FEValues<dim> fe_values(fe_dam,
                        quadrature_formula,
                        update_values | update_gradients |
                        update_quadrature_points | update_JxW_values);
    
    std::vector<double> val_dam(quadrature_formula.size());
    std::vector<Tensor<1, dim>> grad_dam(quadrature_formula.size());
    
    for (auto &cell : dof_handler_dam.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        PointHistory<spacedim> *local_quadrature_points_history =
          reinterpret_cast<PointHistory<spacedim> *>(cell->user_pointer());
        Assert(local_quadrature_points_history >=
                 &quadrature_point_history.front(),
               ExcInternalError());
        Assert(local_quadrature_points_history <=
                 &quadrature_point_history.back(),
               ExcInternalError());


        fe_values.reinit(cell);
        fe_values.get_function_values(dam_NR, val_dam);
        fe_values.get_function_gradients(dam_NR, grad_dam);
        
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
        {
          Psi_tot = Psi_total( local_quadrature_points_history[q].eps_e_NR,
                             grad_dam[q],
                             val_dam[q],
                             local_quadrature_points_history[q].g_hat,
                             local_quadrature_points_history[q].alpha_NR);
          
          Energy += Psi_tot * fe_values.JxW(q);
        } 
      }
    Energy= Utilities::MPI::sum(Energy, mpi_communicator);
    return Energy;
  }
  
  // Compute trace
  template <int dim>
  double trace_rank2(const Tensor<2,dim> &mat)
  {
  	double tr= 0;
  	for(unsigned int i = 0; i <dim; i++)
  		tr +=mat[i][i];
  	
  	return tr;
  }
  
  template <int dim>
  double trace_rank2(const SymmetricTensor<2,dim> &mat)
  {
  	double tr= 0;
  	for(unsigned int i = 0; i <dim; i++)
  		tr +=mat[i][i];
  	
  	return tr;
  }
  
  template <int spacedim>
  double trace(const Tensor<2,spacedim> &mat)
  {
  	double tr= 0;
  	for(unsigned int i = 0; i <spacedim; i++)
  		tr +=mat[i][i];
  	
  	return tr;
  }
  
  // Compute deviatoric stress
  template <int spacedim>
  Tensor<2,spacedim> deviatoric(const Tensor<2,spacedim> &ten)
  {
  	Tensor<2,spacedim> dev_ten;
  	double tr = trace(ten);
  	  	  	
  	for(unsigned int i=0; i<spacedim; i++)
      for(unsigned int j=0; j<spacedim; j++)
        dev_ten[i][j] = ten[i][j] - (1.0/double(spacedim))* ((i==j)?tr:0.0);
  	
    return dev_ten;
  }
  
  
/*  template <int dim>///////? deviatoric not working deviator working better why?
  Tensor<2,dim> deviatoric(const Tensor<2,dim> &ten)
  {
  	Tensor<2,dim> dev_ten;
  	double trace = trace_rank2(ten);
  	  	  	
  	for(unsigned int i=0; i<dim; i++)
      for(unsigned int j=0; j<dim; j++)
        dev_ten[i][j] = ten[i][j] - (1.0/double(dim))* ((i==j)?trace:0.0); ///dim was problem
  	
    return dev_ten;
  }*/
  // Calculate stress strain tensor constant
  template <int dim, int spacedim>
  Tensor<4,spacedim> PFC<dim, spacedim>::calc_strss_strain_tensor_constant(
                                                const Tensor<2, spacedim> &eps_e,
                                                const double &g_hat)
  {
    Tensor<4, spacedim> Cijkl;
    double tr = trace(eps_e);
    //std::cout << g_hat;
    if (tr < 0)
    {
      for(unsigned int i=0; i <spacedim; i++)
        for(unsigned int j=0; j <spacedim; j++)
          for(unsigned int k=0; k <spacedim; k++)
            for(unsigned int l=0; l <spacedim; l++)
              Cijkl[i][j][k][l] = 2.0*MU*g_hat*(0.5*(((i==k) && (j==l) ? 1. : 0.)
                                  + ((i==l) && (j==k) ? 1. : 0.))
                                  -(1./double(spacedim))*((i==j) && (k==l)? 1. : 0.))
                                  + Kn*((i==j) && (k==l) ? 1. : 0.);
    }
    else 
    {
      for(unsigned int i=0; i <spacedim; i++)
        for(unsigned int j=0; j <spacedim; j++)
          for(unsigned int k=0; k <spacedim; k++)
            for(unsigned int l=0; l <spacedim; l++)
              Cijkl[i][j][k][l] = g_hat*(2.*MU*(0.5*(((i==k) && (j==l) ? 1. : 0.)
                                  + ((i==l) && (j==k) ? 1. : 0.))
                                  -(1./double(spacedim))*((i==j) && (k==l) ? 1. : 0.))
                                  + Kn*((i==j) && (k==l) ? 1.: 0.));
    }
    return Cijkl;
  }
  // Calculate stress from strain
  template <int dim, int spacedim>
  Tensor<2, spacedim> PFC<dim, spacedim>::calc_strss_from_strain(const Tensor<4, spacedim> &elasticity_tensor, const Tensor<2, spacedim> &eps_e)
  {
    // Cijkl * eps_e
    Tensor<2, spacedim> stress;
    
    for(unsigned int i=0; i <spacedim; i++)
      for(unsigned int j=0; j <spacedim; j++)
        for(unsigned int k=0; k <spacedim; k++)
          for(unsigned int l=0; l <spacedim; l++)
            stress[i][j] += elasticity_tensor[i][j][k][l] * eps_e[k][l];
    
    
    return stress;
  }

  // Psi plus
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::Psi_plus(const Tensor<2, spacedim> &eps_e)
  {
    double Psi_p;
    double tr = trace(eps_e);
    Tensor<2, spacedim> dev_eps_e = deviatoric(eps_e);
    tr = (tr < 0) ? 0 : tr;
    Psi_p = 0.5*Kn*(tr*tr) + MU * pow((dev_eps_e.norm()),2.);    
    return Psi_p;
  }
  // Psi minus
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::Psi_minus(const Tensor<2, spacedim> &eps_e)
  {
    double Psi_m;
    double tr = trace(eps_e);
    tr = (tr < 0) ? tr : 0;
    Psi_m = 0.5*Kn*(tr*tr);
    return Psi_m;
  }
  // Psi total
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::Psi_total(const Tensor<2, spacedim> &eps_e,
     const Tensor<1,dim> &grad_dam,
     const double &val_dam,
     const double &g_hat,
     const double &alpha)
  {
    double Psi_tot;
    double norm_dam = grad_dam.norm();
    Psi_tot = g_hat*Psi_plus(eps_e) + Psi_minus(eps_e) 
              + (Ystrss*alpha + 0.5*H_pl*alpha*alpha)
              + Gc*(((1-val_dam)*(1-val_dam))/(4*L_dam)+L_dam*norm_dam*norm_dam);
    return Psi_tot;
  }
  // Compute yield function
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::yield_function(const Tensor<2, spacedim> &strss, const double &alpha)
  {
  	double yield_fn;
  	Tensor<2, spacedim> dev_strss = deviatoric(strss);
  	//double sec_inv = -1.*second_invariant(dev_strss);
  	
  	//if (sec_inv < 0)
  	//std::cout << sec_inv << "\t";
  	
  	yield_fn = sqrt(3.0/2.0)*dev_strss.norm() - (Ystrss + H_pl*alpha);
  	//yield_fn = sqrt(3.* sec_inv) - (Ystrss + H_pl*alpha);
  	return yield_fn;
  }
  // Compute Delta lambda
  template <int dim, int spacedim>
  double PFC<dim, spacedim>::delta_lambda (const double &yield_fn, const double &g_hat)
  {
    double del_lambda;
    del_lambda = yield_fn/(H_pl+3.*g_hat*MU);
    
    return del_lambda;
  }
  // Predictor Corrector for u
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::predictor_corrector (const Tensor<2, spacedim> &strain,
                                                const Tensor<2, spacedim> &eps_e_N,
                                                const Tensor<2, spacedim> &eps_p_N,
                                                const double &alpha_N,
                                                const double &const_g_hat,
                                                double &yield_fn,
                                                double &del_lambda,
                                                double &alpha,
                                                Tensor<2, spacedim> &strss,
                                                Tensor<2, spacedim> &eps_e,
                                                Tensor<2, spacedim> &eps_p,
                                                Tensor<4, spacedim> &Cijkl)
  {
    //std::cout<<  yield_fn;
    //const Tensor<2, spacedim> Deps = strain
    //const Tensor<2, spacedim> eps_e_tr = eps_e_N + Deps;
    const Tensor<2, spacedim> eps_e_tr = strain - eps_p_N;
    const Tensor<2, spacedim> eps_p_tr = eps_p_N;
    
    //2 steps for stress calculation should be consistent, Cijkl and then strss
    const Tensor<4, spacedim> tmp_Cijkl = calc_strss_strain_tensor_constant(eps_e_tr,
                                                                            const_g_hat);
    
    const Tensor<2, spacedim> strss_tr = calc_strss_from_strain (tmp_Cijkl, eps_e_tr);
    
    const double alpha_tr = alpha_N;
    yield_fn = yield_function (strss_tr, alpha_tr);
    
    if (yield_fn <=0)
    {
      del_lambda = 0;
      eps_e = eps_e_tr;
      eps_p = eps_p_tr;
      alpha = alpha_tr;
      
      strss = strss_tr;// e_N, e, e_tr
      Cijkl = tmp_Cijkl; //e_N or just e? so Cijkl=Cijkltmp
      //ssFile <<Cijkl<<std::endl;
    }
    else
    {
      del_lambda = delta_lambda(yield_fn, const_g_hat); //const_g_hat
      Tensor<2, spacedim> dev_strss_tr = deviatoric(strss_tr);
      double              strss_norm_tr = dev_strss_tr.norm();
      Tensor<2, spacedim> strss_normal_tr = dev_strss_tr / strss_norm_tr;
      Tensor<2, spacedim> eps_e_dev, eps_e_vol;
      Tensor<4, spacedim> Cijkl_elas, Cijkl_elas_NR, Cijkl_tr, Cijkl_inter, Cijkl_inter1;
      
      // Tensor scalar multiplication from left or right is same
      eps_p = eps_p_tr + del_lambda * sqrt (3.0 / 2.0) * strss_normal_tr;
      alpha = alpha_tr + del_lambda;
      //ssFile <<del_lambda<<std::endl;// Somewhere here is the problem ?
      eps_e_dev = deviatoric(eps_e_tr);
      eps_e_vol = eps_e_tr - eps_e_dev;
      eps_e_dev = eps_e_dev - del_lambda * sqrt (3.0 / 2.0) * strss_normal_tr;
      eps_e = eps_e_vol + eps_e_dev;
      
      Cijkl_elas_NR = calc_strss_strain_tensor_constant(eps_e_N,
                                                        const_g_hat);
      Cijkl_elas = calc_strss_strain_tensor_constant(eps_e,
                                                     const_g_hat);//?
      
      strss = calc_strss_from_strain(Cijkl_elas, eps_e);//?
      //
      int Copt = 2; //2
      if (Copt == 0)
        Cijkl = Cijkl_elas_NR;
      if (Copt == 1)
      {
        Cijkl_tr = calc_strss_strain_tensor_constant(eps_e_tr,
                                                     const_g_hat);
        // Cijqr(del_qs del_rt - factor *()_qrmn Cmnop del_os del_pt)
        // factor = sqrt(1.5)/(h + 3 mu g_hat)
        // ()_qrmn = f del_qm del_rn/ strss_norm_tr +
        //            () dev_strss_tr_qr * dev_strss_tr_mn/2  -
        //             f del_qr del_mn/spacedim/strss_norm_tr
        // () = sigma_y + h alpha_tr
        
        for(unsigned int i=0; i <spacedim; i++)
          for(unsigned int j=0; j <spacedim; j++)
            for(unsigned int k=0; k <spacedim; k++)
              for(unsigned int l=0; l <spacedim; l++)
                Cijkl_inter1[i][j][k][l] = yield_fn*((i==k) && (j==l) ? 1.: 0.)/strss_norm_tr +
                                           .5* (Ystrss + H_pl*alpha_tr) * dev_strss_tr[i][j] * dev_strss_tr[k][l] / std::pow(strss_norm_tr, 3.) -
                                           yield_fn*((i==j) && (k==l) ? 1.: 0.)/(double(spacedim)*strss_norm_tr); 
        
        Cijkl_inter1 *= (std::sqrt(1.5)/(H_pl+3.*const_g_hat*MU)); 
        
        Cijkl_inter = double_contract(Cijkl_inter1, Cijkl_tr); 
        
        for(unsigned int i=0; i <spacedim; i++)
          for(unsigned int j=0; j <spacedim; j++)
            for(unsigned int k=0; k <spacedim; k++)
              for(unsigned int l=0; l <spacedim; l++)
                Cijkl_inter[i][j][k][l] -=  ((i==k) && (j==l) ? 1.: 0.);
       
        Cijkl_inter *= -1.;
         
        Cijkl = double_contract(Cijkl_elas_NR, Cijkl_inter);//?
      }
      
      if (Copt == 2)
      {
        // Cijqr(del_qs del_rt - factor *()_qrmn Cmnop del_os del_pt)
        // factor = sqrt(1.5)/(h + 3 mu g_hat)
        // ()_qrmn = f del_qm del_rn/ strss_norm_tr +
        //            () dev_strss_tr_qr * dev_strss_tr_mn/2  -
        //             f del_qr del_mn/spacedim/strss_norm_tr
        // () = sigma_y + h alpha_tr
        
        for(unsigned int i=0; i <spacedim; i++)
          for(unsigned int j=0; j <spacedim; j++)
            for(unsigned int k=0; k <spacedim; k++)
              for(unsigned int l=0; l <spacedim; l++)
                Cijkl_inter[i][j][k][l] = yield_fn*((i==k) && (j==l) ? 1.: 0.)/strss_norm_tr +
                                           .5* (Ystrss + H_pl*alpha_tr) * dev_strss_tr[i][j] * dev_strss_tr[k][l] / std::pow(strss_norm_tr, 3.) -
                                           yield_fn*((i==j) && (k==l) ? 1.: 0.)/(double(spacedim)*strss_norm_tr); 
        
        Cijkl_inter *= (2.*const_g_hat*MU*std::sqrt(1.5)/(H_pl+3.*const_g_hat*MU)); 
        
        //Cijkl_inter = double_contract(Cijkl_inter1, Cijkl_tr); 
        
        for(unsigned int i=0; i <spacedim; i++)
          for(unsigned int j=0; j <spacedim; j++)
            for(unsigned int k=0; k <spacedim; k++)
              for(unsigned int l=0; l <spacedim; l++)
                Cijkl_inter[i][j][k][l] -=  ((i==k) && (j==l) ? 1.: 0.);
       
        Cijkl_inter *= -1.;
         
        Cijkl = double_contract(Cijkl_elas, Cijkl_inter);//?
      }
      
      if (Copt == 3)
      {
        // from the website
        
        Cijkl_inter = dyadic_product(strss_normal_tr, strss_normal_tr);
        Cijkl_inter *= (2 * MU * const_g_hat * const_g_hat *
                        ((3 * MU /(H_pl+3.*const_g_hat*MU))- 
                        std::sqrt(6)* MU * del_lambda/strss_norm_tr));
        
        for(unsigned int i=0; i <spacedim; i++)
          for(unsigned int j=0; j <spacedim; j++)
            for(unsigned int k=0; k <spacedim; k++)
              for(unsigned int l=0; l <spacedim; l++)
                Cijkl_inter1[i][j][k][l] = 0.5*(((i==k) && (j==l) ? 1. : 0.)
                                           + ((i==l) && (j==k) ? 1. : 0.))
                                           -(1./double(spacedim))*((i==j) && (k==l)? 1. : 0.);
        
        Cijkl_inter1 *= (2 * MU * std::sqrt(6)* MU * del_lambda * const_g_hat/strss_norm_tr);
        
        
         
        Cijkl = Cijkl_elas_NR -  Cijkl_inter - Cijkl_inter1;//?
      }
    }
    
  }
  
  // RUN
  template <int dim, int spacedim>
  void PFC<dim, spacedim>::run()
  {
 /*   do_initial_timestep();
*/	
	  double system_rhs_disp_norm0(0.),
	         system_rhs_dam_norm0(0.),
	         system_rhs_norm1(0.);
	  
	  double Energy0(0.),
	         Energy1(0.),
	         Energy_err(0.);// Energy and relative energy error
	  
	  //int flag_E, flag_U, flag_S;// use incase not able to capture error
	  int flag_BC = 0;
	  	  
	  if (timestep_no == 0)
	  {	
		  pcout<<"\n t=0 \n .. Making grid..."<<std::endl;

      make_grid();

		  pcout<<"\n Setting up system ..." <<std::endl; 
      setup_system(0);
		  
      pcout<<"\n Setting up Qp history ..."<<std::endl; 
      setup_Qp_history();
		  
		  if(this_mpi_process == 0)
		    ssFile << "Time (s)\t" 
               << "Displacement (mm)\t" 
               << "Force (kN/mm)\t"
               << "Avg. eps_p\t"
               << "Avg. eps_p_eq\t"
               << "Avg. p\t"
               << "Avg. g_hat\t"
               << "Avg. damage\t"
               << "Strain\t" 
		           << "Stress" << std::endl;
		           //<< "Shear Strain\t"
		           //<< "Shear Stress\t" 
		           //<< "Bulk Strain\t"
		           //<< "Bulk Stress" << std::endl;
		  
		  Energy0 = compute_Energy();
		  
		  output_results(timestep_no);
		  output_Qp_results(timestep_no);
		  output_stress_strain_plot();
		  
	  }
	//output_results();// Check
	
	//print_mesh_info(triangulation,"test.vtu");
		
		//delete it later
		/*{
		const std::string pl = ("debug/plasticcompare-" + Utilities::int_to_string(this_mpi_process, 3) +
      	 ".txt");
    plastic_strain.open(pl);
		}*/
		
		//#ifdef PFC_DEBUG_STAGGER
    //  const std::string stag_list = "debug/stagger-list.txt";
    //  debug_stagger_list.open(stag_list);
    //#endif
		
    while (time < end_time) 
    {
    	
    	
    	
    	// Monotonic or Fatigue
    	if (load_type == 0)
    	{
      	present_strain += (dt * strain_rate);
      	velocity = strain_rate * size_y;
    	}
    	else if (load_type == 1)
    	{
      	
      	//strain_rate = 0.09* (std::cos(time) +1) ;
      	
      	if (enable_linear && (present_strain < (strain_max + strain_min) /2.))
      	{
      	  present_strain += (dt * strain_rate);
      	  velocity = strain_rate * size_y;
      	}
      	else
      	{ 
      	  enable_linear = false;
      	  if (!stored_time)
      	  {
      	    time_sin_start = time;
      	    stored_time = true;
      	  }
        	velocity = ((strain_max - strain_min) /2.) * std::cos((time - time_sin_start)/0.5) * size_y;
        	present_strain += velocity * dt / size_y;
      	}
    	}
    	else
    	  AssertThrow(load_type < 2, ExcMessage("Unknown load type"));
    	
    	time += dt;
    	timestep_no += 1;
    	
    	//debug stagger
    	//plastic_strain << "\nTimestep: " << timestep_no << std::endl;
    	//#ifdef PFC_DEBUG_STAGGER
      //	const std::string stag = ("debug/stagger-" + Utilities::int_to_string(timestep_no, 4) +
      //	 ".txt");
      //	debug_stagger.open(stag);
      //	if(this_mpi_process == 0)
      //	  debug_stagger_list << stag <<std::endl;
    	//#endif
		  //
		  //#ifdef PFC_DEBUG_NR_U
      //  const std::string newrap_disp_list = ("debug/NR_U-" +
      //    Utilities::int_to_string(timestep_no, 4) + ".txt");
      //  debug_nr_disp_list.open(newrap_disp_list);
      //#endif
      //#ifdef PFC_DEBUG_NR_S
      //  const std::string newrap_dam_list = ("debug/NR_S-" +
      //    Utilities::int_to_string(timestep_no, 4) + ".txt");
      //  debug_nr_dam_list.open(newrap_dam_list);
      //#endif
    	pcout << "\n \n \nTime step: " << timestep_no 
    	      << "\t Time "<< time
    	      << "\t Strain "<< present_strain << std::endl;
      
      #ifdef PFC_DEBUG
        outfile << "\n Time step: " << timestep_no 
                << "\t Time "<< time<<std::endl;
        //matrixfile << "\n Time step: " << timestep_no 
        //           << "\t Time "<< time<<std::endl; 
        //cellmatrixfile << "\n Time step: " << timestep_no 
        //               << "\t Time "<< time<<std::endl;
      #endif      

      //pcout << "\n dt: " << dt << std::endl;
      
      //const std::string st = "results/NRconv" + std::to_string(time) + ".txt";
      //outfile.open(st);
    	// Staggered stepping
      system_rhs_disp_norm0 = 0;
    	system_rhs_dam_norm0 = 0;
      for (unsigned int iter_stagger = 0; iter_stagger <= N_stagger; iter_stagger++)
      {	
      	
      	pcout << "\n Staggered iteration : " << iter_stagger;
        
      	#ifdef PFC_DEBUG
          outfile << "\n Staggered iteration:"<< iter_stagger << std::endl; 
          //cellmatrixfile << "\n Staggered iteration:"<< iter_stagger << std::endl;
        #endif
      	
      	//plastic_strain << "NR: " << iter_stagger << std::endl;
      	// #ifdef PFC_DEBUG_NR_U
       //    const std::string newrap_dispj_list = ("debug/NR_U-" +
       //      Utilities::int_to_string(timestep_no, 4) + "." + 
       //      Utilities::int_to_string(iter_stagger, 4) + ".txt");
       //    debug_nr_disp.open(newrap_dispj_list);
       //    if(this_mpi_process == 0)
      	//     debug_nr_disp_list << newrap_dispj_list <<std::endl;
       //  #endif
      	// NR for u
      	#ifdef PFC_DISP
      	
      	for (unsigned int iter_NR_disp = 0; iter_NR_disp <= N_NR_disp; iter_NR_disp++) // NR solve of displacement
      	{
      	  
         	pcout << "\n Displacement NR iteration : " << iter_NR_disp;
         	setup_system(1);  // reinit system_matrix_u and system_rhs_u
      	  // compute strain  // quadrature computation  // quadrature update
      	        	  
      	  ////assemble_system_u(iter_NR_u);
      	  if ((iter_stagger == 0) && (iter_NR_disp == 0))
      	    flag_BC = 0;
    	    else
    	      flag_BC = 1;

          //pcout<<"\n Assembly - displacement .. started ..."<<std::endl; 
      	  
      	  assemble_system_disp(flag_BC); //

          //pcout<<"\n Assembly - displacement .. finished ..."<<std::endl; 
      	  
      	  if(!(iter_stagger + iter_NR_disp))
            system_rhs_disp_norm0= system_rhs_disp.l2_norm();
      	  
          //pcout<<"\n Solve - displacement .. started ..."<<std::endl; 

      	  solve_NR_disp();
      	
          //pcout<<"\n Solve - displacement NR .. finished ..."<<std::endl; 

          pcout<<"\n Displacement residual norm:"<<system_rhs_disp.l2_norm()<<"\t Displacement residual norm 0:"<<system_rhs_disp_norm0<<std::endl; 

          //pcout << "Blash!" << system_rhs_norm0 << std::endl;
      	  system_rhs_norm1= system_rhs_disp.l2_norm();
      	  #ifdef PFC_DEBUG
      	    //if(this_mpi_process == 0)
    	      /*
    	      outfile << "\n Displacement NR iteration:"<<iter_NR_disp 
    	              << "\t Displacement RHS norm:" << system_rhs_norm1 
    	              <<"\t Displacement residual norm 0:"<<system_rhs_disp_norm0<<std::endl;
            */ 
      	  #endif
      	  
      	  /*if((system_rhs_norm1 < 1.e-1) ||*/ 
      	  if(    ((system_rhs_norm1 / system_rhs_disp_norm0) < rhs_disp_tol) || system_rhs_norm1 < 1e-7)
      	  {
      	    pcout << " Residual norm displacement is :"<< system_rhs_norm1 
      	          <<"\t Steps: " << iter_NR_disp <<"\n"<< std::endl;
      	    // indice transfer on convergence
          	
          	update_Qp_history_disp_NR();
      	    break;
      	  }
      	  //system_rhs_norm0 = system_rhs_u.l2_norm();// STrong
      	  //outfile << iter_NR_u << "\t" << system_rhs_norm0 <<std::endl;
      	  //Assert statement required on not converging
      	  //AssertThrow(iter_NR_disp < N_NR_disp, ExcMessage("NR for displacement did not converge, increase step max. or increase tolerance, or check the rhs total lower limit"));//Only works in debug mode
      	 
      	  
      	}
      	update_Qp_history_disp_NR();// Changed Not tested
      	#endif
      	
      	//#ifdef PFC_DEBUG_NR_U
      	//  debug_nr_disp.close();
      	//#endif
      	
      	
      	// NR for s
/* */   
        //#ifdef PFC_DEBUG_NR_S
        //  const std::string newrap_damj_list = ("debug/NR_S-" +
        //    Utilities::int_to_string(timestep_no, 4) + "." + 
        //    Utilities::int_to_string(iter_stagger, 4) + ".txt");
        //  debug_nr_dam.open(newrap_damj_list);
        //  if(this_mpi_process == 0)
      	//    debug_nr_dam_list << newrap_damj_list <<std::endl;
        //#endif 	
      	
      	#ifdef PFC_DAMAGE
      	
      	for (unsigned int iter_NR_dam = 0; iter_NR_dam <= N_NR_dam; iter_NR_dam++)
      	{
      	  
          pcout << "\n Damage NR iteration : " << iter_NR_dam;

          setup_system(1);  // reinit system_matrix_s and system_rhs_s
      	  // quadrature computation  // quadrature update not needed // all explicit vars
      	  
      	  #ifdef PFC_DEBUG
      	    outfile << "\n Damage NR iteration:\t"<< iter_stagger << std::endl;
      	    //cellmatrixfile << "\n Damage NR iteration:\t"<< iter_stagger << std::endl;
      	  #endif
      	  
      	  assemble_system_dam();
          
          #ifdef PFC_DEBUG
            //matrixfile << "\n System matrix damage:" <<std::endl;
            //system_matrix_dam.print(matrixfile);
            //matrixfile << "\n System rhs damage:" <<std::endl;
            //system_rhs_dam.print(matrixfile, 5);
          #endif
          //pcout<<"\n Assembly - damage .. finished ..."<<std::endl; 
      	  
      	  if(!(iter_stagger + iter_NR_dam))
          {
            system_rhs_dam_norm0= system_rhs_dam.l2_norm();
            //pcout << " Residual norm damage at 0 is :"<< system_rhs_s_norm0 <<"\n"<< std::endl;
          }
            
      	  // del s solve
      	  solve_NR_dam();

          //pcout<<"\n Solve - damage NR .. finished ..."<<std::endl; 

          pcout<<"\n Damage residual norm:"<<system_rhs_dam.l2_norm()<<std::endl; 
      	  
      	  // check Residue s and EXIT
      	  system_rhs_norm1= system_rhs_dam.l2_norm();
      	  #ifdef PFC_DEBUG
      	    //if(this_mpi_process == 0)
      	    outfile << "\n Damage NR iteration:"<<iter_NR_dam 
    	              << "\t Damage RHS norm:" << system_rhs_norm1 
    	              <<"\t Damage residual norm 0:"<<system_rhs_dam_norm0<<std::endl;
      	      //outfile << iter_NR_dam << "\t" << system_rhs_norm1 <<std::endl;
      	  #endif
      	  
      	  if (iter_NR_dam >= N_NR_dam)
      	    enable_dam_debug = true;
      	  
      	  if (enable_dam_debug)
      	  {
        	  
        	  const std::string matrixfilename = ("debug/NR_S_matrix-" +
            Utilities::int_to_string(timestep_no, 4) + "." + 
            Utilities::int_to_string(iter_stagger, 3) + "." + 
            Utilities::int_to_string(iter_NR_dam, 2) + "." + 
            Utilities::int_to_string(this_mpi_process, 3) + ".txt");
        	  
        	  const std::string residual_vectorfilename = ("debug/NR_S_residue-" +
            Utilities::int_to_string(timestep_no, 4) + "." + 
            Utilities::int_to_string(iter_stagger, 3) + "." + 
            Utilities::int_to_string(iter_NR_dam, 2) + "." + 
            Utilities::int_to_string(this_mpi_process, 3) + ".txt");
        	  
        	  const std::string data_vectorfilename = ("debug/NR_S_vector-" +
            Utilities::int_to_string(timestep_no, 4) + "." + 
            Utilities::int_to_string(iter_stagger, 3) + "." + 
            Utilities::int_to_string(iter_NR_dam, 2) + "." + 
            Utilities::int_to_string(this_mpi_process, 3) + ".txt");
            
            matrixfile.open(matrixfilename);
            residual_vectorfile.open(residual_vectorfilename);
            data_vectorfile.open(data_vectorfilename);
            
            
            double stiffness_norm = 0;
            double vector_norm = 0;
            
            stiffness_norm = system_matrix_dam.frobenius_norm();
            
            vector_norm = dam_NR.l2_norm();
            
            //vector_norm = Utilities::MPI::sum(vector_norm, mpi_communicator);
            //vector_norm = std::sqrt(vector_norm);
            
            if(this_mpi_process == 0)
            {
          	  matrixfile << "Timestep no\t" << timestep_no 
          	             << "\tStagger iter:\t" << iter_stagger
          	             << "\tNR iter:\t" << iter_NR_dam;
          	  matrixfile << "\nStiffness matrix:\t" << stiffness_norm << "(norm)\n";
              residual_vectorfile << "Residue vector:\t" << system_rhs_norm1 << "(norm)\n";
              data_vectorfile << "Data vector:\t" << vector_norm << "(norm)\n";
              
              if (iter_NR_dam == 0)
              {
                outfile << "Timestep no:\t" << timestep_no 
        	              << "\tStagger iter:\t" << iter_stagger << "\n";
                outfile << "NR iter\t" 
                        << "Stiffness norm\t"
                        << "Residual norm\t"
                        << "Vector norm\n";
              }
              
              outfile << iter_NR_dam << "\t" 
                      << stiffness_norm << "\t"
                      << system_rhs_norm1 << "\t"
                      << vector_norm << "\n";
              
            }
            
            system_matrix_dam.print(matrixfile);
            system_rhs_dam.print(residual_vectorfile, 5 , true, false);
            dam_NR.print(data_vectorfile, 5 , true, false);
            
            matrixfile.close();
            residual_vectorfile.close();
            data_vectorfile.close();
            
      	  }
      	  
      	  if((system_rhs_norm1 < 1.e-7) || 
      	  /*if(*/    ((system_rhs_norm1 / system_rhs_dam_norm0) < rhs_dam_tol))
      	  {
      	    pcout << "\n Residual norm damage is :"<< system_rhs_norm1 
      	          << "\t Steps: "<<iter_NR_dam <<"\n"<< std::endl;
      	    update_Qp_history_dam_NR();
      	    break;
      	  }
      	  //pcout<<"Not converging"<<std::endl;
      	  /*if (iter_NR_s >= N_NR_s)
      	    throw(ExcMessage("NR for damage did not converge, increase step max. or increase tolerance"));*/
      	    
      	  //AssertThrow(iter_NR_dam < N_NR_dam, ExcMessage("NR for damage did not converge, increase step max. or increase tolerance, or check the rhs total lower limit"));
      	  //Only works in debug mode 
      	  
      	  //#ifdef PFC_DEBUG
      	  
      	  
      	  //#endif
      	}
      	
      	update_Qp_history_dam_NR();// Changed
      	#endif
      	
      	//#ifdef PFC_DEBUG_NR_S
      	//  debug_nr_dam.close();
      	//#endif
      	// indice transfer
      	      	
      	// check energy and EXIT
      	if(!(iter_stagger))
      	  Energy0 = compute_Energy();
      	
      	Energy1 = compute_Energy();
      	//pcout << "\n Energy minimum is :"<< E1 ;//if(!0)// Temporary replacement
      	Energy_err = std::abs((Energy1 - Energy0)/Energy0);
      	
      	#ifdef PFC_DEBUG
        	//if(this_mpi_process == 0)
		        /*outfile << "\n Staggered iteration:"<<iter_stagger 
		                << "\t Energy error:" << Energy_err 
		                << "\t Energy value:" << Energy1 << std::endl;*/
      	#endif
      	// #ifdef PFC_DEBUG_NR_U
      	//   debug_nr_disp_list.close();
      	// #endif
      	// #ifdef PFC_DEBUG_NR_S
      	//   debug_nr_dam_list.close();
      	// #endif
      	if( ((Energy1 <=1e-25) || (Energy_err < Energy_tol)) && (iter_stagger!= 0)) // transfer all variables
    	  {
    	    calc_timestep();
    	    accept_NR_solution();
    	    accept_NR_state();
    	    
    	    pcout << "\n Staggered Steps converged with Energy minimum: "<< Energy1;
    	    break;
    	  }
    	  else // if not converged then reset k variables while keeping g_hat and Psi_plus_max
    	  {
    	    pcout << "\t Current Energy:\t" << Energy1 << std::endl;
    	    //reset_vark();
    	    //transfer_vark();
    	    Energy0 = Energy1;
    	  }//Assert statement required on not converging
    	  AssertThrow(iter_stagger < N_stagger, ExcMessage("Staggered step did not converge, increase step max. or increase Energy tolerance")); //Only works in debug mode
      }
      // Update energy
      
      output_results(timestep_no);
      output_Qp_results(timestep_no);
      output_stress_strain_plot();
      
      // #ifdef PFC_DEBUG_STAGGER
      //   debug_stagger.close();
      // #endif
      
    }//?MArk
    // #ifdef PFC_DEBUG_STAGGER
    //   debug_stagger_list.close();
    // #endif
    
    //plastic_strain.close();
  }


  
}




int main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Problem;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
	  std::cout << "Running" << std::endl;
      PFC<2, 3> fracture_problem;
      fracture_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
