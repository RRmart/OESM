This folder contains the MATLAB program developed for the OESM PROJECT P2.

Files for topology optimization (models -load and boundary condition cases - explained in report file):

"AugLag_model1.m" - Uses augmented Lagrangean method with a Q4 element mesh to obtain results for example model 1.

"OC_Q4_model1.m" - Uses OC method with a Q4 element mesh to obtain results for model 1.

"OC_Q4_model2.m" - Uses OC method with a Q4 element mesh to obtain results for model 2.

"OC_Q9_model1.m" - Uses OC method with a Q9 element mesh to obtain results for model 1.

"OC_Q9_model2.m" - Uses OC method with a Q9 element mesh to obtain results for model 2.

Other Files

"Bisection.m" - Implements the bisection method within a certain interval.

"Deriv.m" - Calculates the derivative of compliance function with regards to the design variables.

"DerivL.m" - Calculates the derivative of the Augmented Lagrangean function with regards to the design variables.

"ElementsQ4.m" - Creates the Q4 mesh

"ElementsQ9.m" - Creates the Q9 mesh

"filtering.m" - Aplies a low pass filter to the sensitivities

"FiniteElement.m" - Runs the finite element method for model 1.

"K_Q4.m" - Calculates the standard element stiffness matrix for Q4 elements.

"K_Q9.m" - Calculates the standard element stiffness matrix for Q9 elements.

"shape_funcQ9.m" - shape functions for the Q9 element (needed for "K_Q9.m"). 