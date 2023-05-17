# Copyright 2020, Jack S. Hale
# SPDX-License-Identifier: LGPL-3.0-or-later
import numpy as np

from mpi4py import MPI
from petsc4py import PETSc

import dolfinx

import dolfinx
import dolfinx.fem as dfem
import dolfinx.mesh as dmesh

import fenicsx_error_estimation

import ufl

# Structured mesh
msh = dmesh.create_rectangle(MPI.COMM_WORLD, [np.array([0, 0]), np.array([1, 1])], [32, 32], cell_type=dmesh.CellType.triangle)

k = 2
element = ufl.FiniteElement("CG", ufl.triangle, k)
V = dfem.FunctionSpace(msh, element)
dx = ufl.Measure("dx", domain=msh)

x = ufl.SpatialCoordinate(msh)
f = 8.0 * ufl.pi**2 * ufl.sin(2.0 * ufl.pi * x[0]) * ufl.sin(2.0 * ufl.pi * x[1])

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
L = ufl.inner(f, v) * dx

u0 = dfem.Function(V)
with u0.vector.localForm() as u0_local:
    u0_local.set(0.0)

facets = dmesh.locate_entities_boundary(msh, msh.topology.dim - 1, 
                                        lambda x: np.full(x.shape[1], True, dtype=bool))
dofs = dfem.locate_dofs_topological(V, msh.topology.dim - 1, facets)
bcs = [dfem.dirichletbc(u0, dofs)]

problem = dfem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
u_h = problem.solve()

with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "u.xdmf", "w") as of:
    of.write_mesh(msh)
    of.write_function(u_h)

u_exact = ufl.sin(2.0 * ufl.pi * x[0]) * ufl.sin(2.0 * ufl.pi * x[1])
form_error = dfem.form(ufl.inner(ufl.grad(u_h - u_exact), ufl.grad(u_h - u_exact)) * dx(degree=k + 3))
error = MPI.COMM_WORLD.allreduce(dfem.assemble_scalar(form_error), op=MPI.SUM)
print("True error: {}".format(np.sqrt(error)))

# Now we specify the Bank-Weiser error estimation problem.
element_f = ufl.FiniteElement("DG", msh.ufl_cell(), k + 1)
element_g = ufl.FiniteElement("DG", msh.ufl_cell(), k)
element_e = ufl.FiniteElement("DG", msh.ufl_cell(), 0)
N = fenicsx_error_estimation.create_interpolation(element_f, element_g)

# Function spaces
V_f = dfem.FunctionSpace(msh, element_f)
e = ufl.TrialFunction(V_f)
v = ufl.TestFunction(V_f)

# Bilinear form
a_e = ufl.inner(ufl.grad(e), ufl.grad(v)) * dx

# Linear form
n = ufl.FacetNormal(msh)
dS = ufl.Measure("dS", domain=msh)
L_e = ufl.inner(ufl.jump(ufl.grad(u_h), -n), ufl.avg(v)) * dS + ufl.inner(f + ufl.div((ufl.grad(u_h))), v) * dx

# Error form
V_e = dfem.FunctionSpace(msh, element_e)
e_h = dfem.Function(V_f)
v_e = ufl.TestFunction(V_e)
L_eta = ufl.inner(ufl.inner(ufl.grad(e_h), ufl.grad(e_h)), v_e) * dx

# Functions to store results
eta_h = dfem.Function(V_e)

# We must apply homogeneous zero Dirichlet condition on the local problems when
# a (possibly non-zero) Dirichlet condition was applied to the original
# problem. Due to the way boundary conditions are enforced locally it is only
# necessary to compute a sorted list of entities (facets) on which homogeneous
# Dirichlet conditions should be applied.
facets_sorted = np.sort(facets)

# Estimate the error using the Bank-Weiser approach.
fenicsx_error_estimation.estimate(eta_h, a_e, L_e, L_eta, N, facets_sorted)

print("Bank-Weiser error from estimator: {}".format(np.sqrt(eta_h.vector.sum())))

with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "output/eta.xdmf", "w") as of:
    of.write_mesh(msh)
    of.write_function(eta_h)