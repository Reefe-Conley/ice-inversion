import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, diags, eye
from scipy.special import expit
from scipy.sparse.linalg import spsolve

def D2_mat(grid):
    """
    Construct a second-derivative matrix for a function defined on a regular grid.
    A standard [1 -2 1] central difference stencil is used in the interior of the matrix.
    
    A [2 -5 4 -1] one-sided stencil is used on the edges of the matrix where central differences
    cannot be calculated. This stencil gives O(Δx^2) accuracy despite being one-sided, at the cost of
    using more function evaluations.

    Parameters
    ----------
    grid : ndarray
        The regular grid on which to construct the D2 matrix.

    Returns
    -------
    D2 : Second-order accurate second-derivative matrix.

    Raises
    ------
    ValueError
        If the grid is not regular.
    """

    # Store grid properties and ensure it is regular
    grid = grid.flatten()
    N = len(grid)
    dx = grid[1] - grid[0]
    if not np.allclose(np.diff(grid), dx):  
            raise ValueError(f'Spatial grid must be regular!') 

    # Construct D2 matrix entries
    lower = 1 * np.ones(N - 1)
    diagonal = -2 * np.ones(N)
    upper = 1 * np.ones(N - 1)
    one_sided = np.array([2, -5, 4, -1])

    # Construct matrix interior with central differencing
    D2 = diags([lower, diagonal, upper], offsets=[-1, 0, 1]).tocsr()
    # Apply one-sided differencing at domain edges
    D2[0,:4] = one_sided
    D2[-1,-4:] = np.flip(one_sided)

    D2 = D2 / dx**2

    return D2

def solve_newton(x, h, b0, zm0, kappa, alpha, eps, D2, max_iter=50, tol=1e-8, print_flag=False):
    """
    Solve for the nondimensional midline profile z_m(x) using Newton's method.

    The nonlinear equation uses logistic approximations to the grounded and
    floating indicator functions. Four Dirichlet boundary conditions are imposed:
    the first two grid points use the grounded condition, and the last two use
    the floating condition.

    Vertical lengths are nondimensionalised by the characteristic vertical scale
    H. The horizontal grid `x` must be nondimensionalised by the chosen horizontal
    scale L (typically the bending–buoyancy length), consistent with the
    nondimensional governing equation and the supplied second-derivative matrix
    `D2`.

    Parameters
    ----------
    x : ndarray, shape (N,)
        One-dimensional nondimensional horizontal grid.
    h : ndarray, shape (N,)
        Nondimensional ice thickness profile on `x`.
    b0 : ndarray, shape (N,)
        Nondimensional bed profile on `x`, using the sign convention of the
        governing equation.
    zm0 : ndarray, shape (N,)
        Initial guess for the nondimensional midline profile.
    kappa : float
        Nondimensional bed-stiffness parameter.
    alpha : float
        Density ratio parameter, typically rho_i / rho_w.
    eps : float
        Positive logistic transition-width parameter, in nondimensional vertical
        units.
    D2 : scipy sparse matrix, shape (N, N)
        Discrete second-derivative operator on `x`, including the appropriate
        grid-spacing scaling.
    max_iter : int, optional
        Maximum number of Newton iterations. Default is 50.
    tol : float, optional
        Convergence tolerance applied to the infinity norm of the Newton update.
        Default is 1e-8.
    print_flag : bool, optional
        If True, print the iteration number and residual norm at each iteration.
        Default is False.

    Returns
    -------
    z : ndarray, shape (N,)
        Converged nondimensional midline profile.
    z_lower : ndarray, shape (N,)
        Lower surface profile, calculated as `z - h/2`.
    z_upper : ndarray, shape (N,)
        Upper surface profile, calculated as `z + h/2`.

    Raises
    ------
    ValueError
        If `x`, `h`, `b0`, and `zm0` do not have matching shapes.
    RuntimeError
        If Newton's method does not converge within `max_iter` iterations.

    Notes
    -----
    The discretized elastic operator is constructed as
    `D2 @ diag(h**3) @ D2`. Boundary residual rows are replaced by the four
    Dirichlet conditions. The convergence check uses the update norm; callers
    should also verify the final equation residual and boundary-condition errors.
    """

    # Check input dimensions match
    if not (x.shape == h.shape == b0.shape == zm0.shape):
         raise ValueError(
              "x, h, b0 and zm0 must all have the same shape;"
              f"got {x.shape}, {h.shape}, {b0.shape}, {zm0.shape}."
         )

    h = np.asarray(h)
    b0 = np.asarray(b0)
    z = np.array(zm0, dtype=float) # Make a new copy of zm0
    N = len(x)
    h3 = h**3
    H3 = diags(h3)
    # Construct elastic operator
    E = D2 @ H3 @ D2

    # Boundary conditions
    grounded_bc = (0.5 - alpha / (1 + kappa))*h - (kappa/(1+kappa))*b0
    floating_bc = (0.5 - alpha)*h
    bc_idx = np.array([0, 1, N-2, N-1])
    bc_val = np.array([grounded_bc[0], grounded_bc[1], floating_bc[-2], floating_bc[-1]])

    # Iteratively solve linearised system using Newton-Raphson
    iterno = 0
    for iteration in range(max_iter):
        q1 = h/2 - z - b0
        q2 = h/2 - z
        S1 = expit(q1/eps)
        S2 = expit(q2/eps)

        # Full nonlinear residual
        F = E @ z + alpha*h - kappa*q1*S1 - q2*S2
        # Apply boundary conditions
        F[bc_idx] = z[bc_idx] - bc_val

        # Calculate Jacobian by differentiating the restoring forces
        drestore = kappa*(S1 + q1/eps * S1 * (1 - S1)) + S2 + q2/eps * S2 * (1 - S2)

        J = E + diags(drestore) # The drestore matrix consists of pointwise forces, so its Jacobian contribution is diagonal.

        # Replace corresponding Jacobian rows with identity rows to satisfy BCs
        J = J.tolil()
        for i in bc_idx:
            J.rows[i] = [i]
            J.data[i] = [1.0]
        J = J.tocsr()

        delta = spsolve(J, -F)
        z += delta
        iterno += 1

        if print_flag == True:
            print(iterno)
            print(np.linalg.norm(F))

        if np.linalg.norm(delta, ord=np.inf) < tol:
             return z, z - h/2, z + h/2

    raise RuntimeError("Newton iteration did not converge")