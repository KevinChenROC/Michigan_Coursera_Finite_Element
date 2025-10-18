/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
It uses the deal.II FEM library, dealii.org*/

//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM1.h"
#include "writeSolutions.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);

    unsigned int order = 1;
    unsigned int problem = 1;

    // Get user input for basis function order
    std::cout << "Enter the basis function order (1, 2, or 3): ";
    std::cin >> order;
    if (order < 1 || order > 3) {
      throw std::invalid_argument("Invalid basis function order. Must be 1, 2, or 3.");
    }

    // Get user input for subproblem
    std::cout << "Enter the subproblem (1 or 2): ";
    std::cin >> problem;
    if (problem < 1 || problem > 2) {
      throw std::invalid_argument("Invalid subproblem. Must be 1 or 2.");
    }

    // Get user input for the number of elements in the mesh
    unsigned int num_elements;
    std::cout << "Enter the number of elements in the mesh: ";
    std::cin >> num_elements;
    if (num_elements == 0) {
      throw std::invalid_argument("Number of elements must be greater than 0.");
    }

    FEM<1> problemObject(order,problem);
    
    //Define the number of elements as an input to "generate_mesh"
    problemObject.generate_mesh(num_elements); //e.g. a 10 element mesh
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    std::cout << problemObject.l2norm_of_error() << std::endl;
    
    //write output file in vtk format for visualization
    problemObject.output_results();
    
    //write solutions to h5 file
    char tag[21];
    sprintf(tag, "CA1_Order%d_Problem%d",order,problem);
    writeSolutionsToFileCA1(problemObject.D, problemObject.l2norm_of_error(), tag);
  }
  catch (std::exception &exc){
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
  catch (...){
    std::cerr << std::endl << std::endl
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
