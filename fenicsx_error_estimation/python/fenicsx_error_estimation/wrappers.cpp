#include "caster_petsc.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <fenicsx_error_estimation/projected_local_solver.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(cpp, m)
{
  // Create module for C++ wrappers
  m.doc() = "DOLFINX interface to fenicsx_error_estimation";

  m.def("projected_local_solver_have_fine_space",
        &projected_local_solver<PetscScalar, true>,
        "Local solves on projected finite element space. Computes"
        "Bank-Weiser error solution. Allows imposition of Dirichlet"
        "boundary conditions on Bank-Weiser error solution.");
  m.def("projected_local_solver_no_fine_space",
        &projected_local_solver<PetscScalar, false>,
        "Local solves on projected finite element space. Does not compute "
        "Bank-Weiser error solution. Dirichlet condition on Bank-Weiser "
        "error solution is always zero.");
}
