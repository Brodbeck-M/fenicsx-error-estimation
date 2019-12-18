import numpy as np

from dolfin import *
from dolfin.fem.assembling import _create_dolfin_form
import bank_weiser.cpp as cpp

def estimate(a_e, L_e, N, bcs=[]):
    try:
        len(bcs)
    except TypeError:
        bcs = [bcs]

    a_e_dolfin = _create_dolfin_form(a_e)
    L_e_dolfin = _create_dolfin_form(L_e)

    V_f = FunctionSpace(L_e_dolfin.function_space(0))
    e_V_f = Function(V_f)

    cpp.estimate(e_V_f.cpp_object(), a_e_dolfin, L_e_dolfin, N, bcs)

    return e_V_f

def estimate_python(a_e, L_e, N, bcs=[]):
    L_e_dolfin = _create_dolfin_form(L_e)
    V_f = FunctionSpace(L_e_dolfin.function_space(0))
    mesh = V_f.mesh()

    e_V_f = Function(V_f)

    # Get boundary dofs in case of Dirichlet boundary conditions
    try:
        bc_dofs = bcs[0].get_boundary_values()
    except IndexError:
        bc_dofs = []

    for cell in cells(mesh):
        A_e = assemble_local(a_e, cell)
        b_e = assemble_local(L_e, cell)

        cell_dofs = V_f.dofmap().cell_dofs(cell.index())

        # Apply Dirichlet BC to local system
        for i, cell_dof in enumerate(cell_dofs):
            try:
                value = bc_dofs[cell_dof]
                A_e[i, :] = 0.0
                b_e -= A_e[:, i]*value
                A_e[:, i] = 0.0
                A_e[i, i] = 1.0
                b_e[i] = value
            except KeyError:
                pass
            except IndexError:
                pass

        # Projection of the local residual system onto proper Bank-Weiser space
        A_e_0 = N.T@A_e@N
        b_e_0 = N.T@b_e

        # Solving local residual system in Bank-Weiser space
        e_0 = np.linalg.solve(A_e_0, b_e_0)
        # Sending residual solution back to fine space
        e = N@e_0

        e_V_f.vector()[cell_dofs] = e

    return e_V_f

def weighted_estimate(eta_uh, eta_zh):
    eta_uh_vec = eta_uh.vector()
    eta_zh_vec = eta_zh.vector()

    sum_eta_uh = eta_uh_vec.sum()
    sum_eta_zh = eta_zh_vec.sum()

    eta_wh = Function(eta_uh.function_space(), name="eta")
    eta_wh.vector()[:] = ((sum_eta_zh/(sum_eta_uh + sum_eta_zh))*eta_uh_vec) + \
                         ((sum_eta_uh/(sum_eta_uh + sum_eta_zh))*eta_zh_vec)

    return eta_wh
