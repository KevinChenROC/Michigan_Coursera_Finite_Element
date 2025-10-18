/*This is a template file for use with 3D finite elements (scalar field).
  The portions of the code you need to fill in are marked with the comment "//EDIT".

  Do not change the name of any existing functions, but feel free
  to create additional functions, variables, and constants.
  It uses the deal.II FEM library.*/

//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace dealii;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM();  // Class constructor 
  ~FEM(); //Class destructor

  //Define your 2D basis functions and derivatives
  double basis_function(unsigned int node, 
			double xi_1,
			double xi_2,
			double xi_3);
  std::vector<double> basis_gradient(unsigned int node, 
				     double xi_1,
				     double xi_2,
				     double xi_3);

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe; 	      //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  //Gaussian quadrature - These will be defined in setup_system()
  unsigned int	      quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double> quad_points; //vector of Gauss quadrature points
  std::vector<double> quad_weight; //vector of the quadrature point weights
    
  //Data structures
  SparsityPattern    	        sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>	        K; 		  //Global stiffness (sparse) matrix
  Vector<double>       	        D, F; 		  //Global vectors - Solution vector (D) and Global force vector (F)
  Table<2,double>	        nodeLocation;     //Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values;  //Map of dirichlet boundary conditions 

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM ()
:
fe (FE_Q<dim>(1), 1), 
  dof_handler (triangulation)
{
  //Nodal Solution names - this is for writing the output file
  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){
  dof_handler.clear ();
}

//EDIT Q1 trilinear basis N_A(ξ,η,ζ) with deal.II node order (see header comment)
template <int dim>
double FEM<dim>::basis_function(unsigned int A, double xi, double eta, double zeta){
  /*"node" specifies which node the basis function corresponds to, 
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the specified basis function and order at the given quadrature pt.
  */
  switch(A){
    case 0: return 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 - zeta); // (-1,-1,-1)
    case 1: return 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 - zeta); // (+1,-1,-1)
    case 2: return 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 - zeta); // (-1,+1,-1)
    case 3: return 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 - zeta); // (+1,+1,-1)
    case 4: return 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 + zeta); // (-1,-1,+1)
    case 5: return 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 + zeta); // (+1,-1,+1)
    case 6: return 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 + zeta); // (-1,+1,+1)
    case 7: return 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 + zeta); // (+1,+1,+1)
    default: return 0.0;
  }
}


/*"node" specifies which node the basis function corresponds to, 
  "xi" is the point (in the bi-unit domain) where the function is being evaluated.
  You need to calculate the value of the derivative of the specified basis function and order at the given quadrature pt.
  Note that this is the derivative with respect to xi (not x)
*/
// ∂N_A/∂(ξ,η,ζ) on the reference element; physical gradient is J^{-T} times this vector
template <int dim>
std::vector<double> FEM<dim>::basis_gradient(unsigned int A, double xi, double eta, double zeta){
  std::vector<double> dN(3,0.0); // [∂N/∂ξ, ∂N/∂η, ∂N/∂ζ]
  // For each A, define signs (sx, sy, sz) ∈ {−1,+1} multiplying (1 + sx*xi), etc.
  double sx=0, sy=0, sz=0;
  switch(A){
    case 0: sx=-1; sy=-1; sz=-1; break;
    case 1: sx=+1; sy=-1; sz=-1; break;
    case 2: sx=-1; sy=+1; sz=-1; break;
    case 3: sx=+1; sy=+1; sz=-1; break;
    case 4: sx=-1; sy=-1; sz=+1; break;
    case 5: sx=+1; sy=-1; sz=+1; break;
    case 6: sx=-1; sy=+1; sz=+1; break;
    case 7: sx=+1; sy=+1; sz=+1; break;
    default: return dN;
  }
  const double fx = (1.0 + sx*xi);
  const double fy = (1.0 + sy*eta);
  const double fz = (1.0 + sz*zeta);
  dN[0] = 0.125 * sx * (      fy) * (      fz); // d/dξ: derivative hits (1+sx*ξ)
  dN[1] = 0.125 * sy * (      fx) * (      fz); // d/dη
  dN[2] = 0.125 * sz * (      fx) * (      fy); // d/dζ
  return dN;
}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //EDIT: Define the limits of your domain
  //Quiz 2 (3D) domain: x ∈ [0,0.04], y ∈ [0,0.08], z ∈ [0,0.02]
   const double x_min = 0.0,   x_max = 0.04;
   const double y_min = 0.0,   y_max = 0.08;
   const double z_min = 0.0,   z_max = 0.02;

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){

  //EDIT - Define the Dirichlet boundary conditions.
	
  /*Note: this will be very similiar to the define_boundary_conds function
    in the HW2 template. You will loop over all nodes and use "nodeLocations"
    to check if the node is on the boundary with a Dirichlet condition. If it is,
    then add the node number and the specified value (temperature in this problem)
    to the boundary values map, something like this:

    boundary_values[globalNodeIndex] = dirichletTemperatureValue

    Note that "nodeLocation" is now a Table instead of just a vector. The row index is
    the global node number; the column index refers to the x, y, or z component (0, 1, or 2 for 3D).
    e.g. nodeLocation[7][2] is the z coordinate of global node 7*/

    boundary_values.clear();

    const double c0 = 1.0/3.0; // K/m
    const double xL = 0.0;
    const double xR = 0.04;
    const double tol = 1e-12;

    const unsigned int totalNodes = dof_handler.n_dofs();
    for (unsigned int a=0; a<totalNodes; ++a){
      const double x = nodeLocation[a][0];
      const double y = nodeLocation[a][1];
      const double z = nodeLocation[a][2];

      if (std::fabs(x - xL) < tol){
        // u(y,z) = 300 * (1 + c0*(y+z)) on x = 0
        const double uD = 300.0*(1.0 + c0*(y + z));
        boundary_values[a] = uD;
      } else if (std::fabs(x - xR) < tol){
        // u(y,z) = 310 * (1 + c0*(y+z)) on x = 0.04
        const double uD = 310.0*(1.0 + c0*(y + z));
        boundary_values[a] = uD;
      }
      // All other faces: natural Neumann (zero flux) — nothing to add here
    }

}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Fill in the Table "nodeLocations" with the x, y, and z coordinates of each node by its global index
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      nodeLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  F.reinit (dof_handler.n_dofs());
  D.reinit (dof_handler.n_dofs());

  //Define quadrature rule - again, you decide what quad rule is needed
  quadRule = 2; //EDIT - Number of quadrature points along one dimension
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = -sqrt(1./3.); //EDIT
  quad_points[1] = sqrt(1./3.); //EDIT

  quad_weight[0] = 1.; //EDIT
  quad_weight[1] = 1.; //EDIT

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  K=0; F=0;

  const unsigned int  	    dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double> 	    Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>      	    Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    /*Retrieve the effective "connectivity matrix" for this element
      "local_dof_indices" relates local dofs to global dofs,
      i.e. local_dof_indices[i] gives the global dof number for local dof i.*/
    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs and quadrature points to populate Flocal and Klocal.
    FullMatrix<double> Jacobian(dim,dim);
    double detJ;

    //Loop over local DOFs and quadrature points to populate Flocal
    Flocal = 0.;
    for(unsigned int q1=0; q1<quadRule; q1++){
      for(unsigned int q2=0; q2<quadRule; q2++){
        for(unsigned int q3=0; q3<quadRule; q3++){
          Jacobian = 0.;
          for(unsigned int i=0;i<dim;i++){
            for(unsigned int j=0;j<dim;j++){
              for(unsigned int A=0; A<dofs_per_elem; A++){
          Jacobian[i][j] += nodeLocation[local_dof_indices[A]][i]
            *basis_gradient(A,quad_points[q1],quad_points[q2],quad_points[q3])[j];
              }
            }
          }
          detJ = Jacobian.determinant();
          for(unsigned int A=0; A<dofs_per_elem; A++){
            //You would define Flocal here if it were nonzero.
          }
        }
      }
    }

    //Loop over local DOFs and quadrature points to populate Klocal
    FullMatrix<double> invJacob(dim,dim), kappa(dim,dim);

    //"kappa" is the conductivity tensor
    kappa = 0.;
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;
    kappa[2][2] = 385.;

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0.;
    // Quadrature over reference cube [-1,1]^3
    for(unsigned int q1=0; q1<quadRule; ++q1){
      for(unsigned int q2=0; q2<quadRule; ++q2){
        for(unsigned int q3=0; q3<quadRule; ++q3){

          // Build Jacobian J_ij = Σ_A x_Ai * ∂N_A/∂ξ_j at this quadrature point
          FullMatrix<double> Jacobian(dim,dim), invJacob(dim,dim);
          Jacobian = 0.0;
          for(unsigned int i=0;i<dim;i++){
            for(unsigned int j=0;j<dim;j++){
              for(unsigned int A=0; A<dofs_per_elem; A++){
                const auto dNA = basis_gradient(A, quad_points[q1], quad_points[q2], quad_points[q3]);
                Jacobian[i][j] += nodeLocation[local_dof_indices[A]][i] * dNA[j];
              }
            }
          }

          const double detJ = Jacobian.determinant();
          invJacob.invert(Jacobian);

          // Precompute weights product
          const double wq = quad_weight[q1]*quad_weight[q2]*quad_weight[q3];

          // Loop over shape functions
          for(unsigned int A=0; A<dofs_per_elem; ++A){
            // ref-grad of N_A
            const auto dNA_ref = basis_gradient(A, quad_points[q1], quad_points[q2], quad_points[q3]);

            // physical grad: grad_x N_A = J^{-T} * grad_ref N_A
            double gradNA[dim] = {0.0,0.0,0.0};
            for(unsigned int I=0; I<dim; ++I)        // physical component
              for(unsigned int j=0; j<dim; ++j)      // ref component
                gradNA[I] += invJacob[j][I] * dNA_ref[j];

            for(unsigned int B=0; B<dofs_per_elem; ++B){
              const auto dNB_ref = basis_gradient(B, quad_points[q1], quad_points[q2], quad_points[q3]);

              double gradNB[dim] = {0.0,0.0,0.0};
              for(unsigned int J=0; J<dim; ++J)
                for(unsigned int k=0; k<dim; ++k)
                  gradNB[J] += invJacob[k][J] * dNB_ref[k];

              // (grad N_A)^T κ (grad N_B)
              double gKg = 0.0;
              for(unsigned int I=0; I<dim; ++I)
                for(unsigned int J=0; J<dim; ++J)
                  gKg += gradNA[I] * kappa[I][J] * gradNB[J];

              Klocal[A][B] += gKg * detJ * wq;
            } //B 
          } //A
        } //q3
      } //q2 
    } // q1
    

    //Assemble local K and F into global K and F
    for(unsigned int A=0; A<dofs_per_elem; A++){
      //You would assemble F here if it were nonzero.
      for(unsigned int B=0; B<dofs_per_elem; B++){
        //EDIT - Assemble K from Klocal (you can look at HW2)
        K.add(local_dof_indices[A], local_dof_indices[B], Klocal[A][B]);
      }
    }

  }

  //Apply Dirichlet boundary conditions
  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}

//Solve for D in KD=F
template <int dim>
void FEM<dim>::solve(){

  //Solve for D
  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}

//Output results
template <int dim>
void FEM<dim>::output_results (){

  //Write results to VTK file
  std::ofstream output1("solution.vtk");
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
			   nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}
