import numpy as np

def analytical_profile(x_start, x_end, Nx, h, height_bb, alpha, kappa, angle_deg):
    """
    Returns the analytical centreline and upper/lower surfaces of a constant-thickness elastic beam
    resting on a Winkler foundation submerged in fluid at z=0 (Sayag and Worster, 2011).

    Parameters
    ----------
    x_start : float
        The leftmost (grounded) domain boundary.
    x_end : float
        The rightmost (floating) domain boundary.
    Nx : int
        The number of grid cells in the domain.
    h : ndarray
        The thickness profile h(x) along the domain.
    height_bb: float
        The ratio of characteristic vertical and horizontal length scales, H/L. L is the bending-buoyancy length scale.
    alpha: float
        The density ratio rho_i / rho_w.
    kappa: float
        The elasticity-buoyancy ratio k / (rho_w g).
    angle_deg: float
        The bed slope angle in degrees.
    """

    x = np.linspace(x_start, x_end, Nx)
    s = np.deg2rad(angle_deg)
    b0 = np.tan(s)* x / height_bb

    # Constants used later in the solution
    gamma0 = height_bb / np.sqrt(2) * kappa**(1/4)
    gamma1 = height_bb / np.sqrt(2)

    # Calculate grounding line position
    xg = height_bb * alpha / np.tan(s) * (1 + kappa**(-1/2)) - np.sqrt(2)/(1 + kappa**(-1/4))
    xg_ind = np.where(x < xg)
    xf_ind = np.where(x >= xg)

    gamma0 = height_bb / np.sqrt(2) * kappa**(1/4)
    gamma1 = height_bb / np.sqrt(2)
    y_float = 1/2 - alpha + np.exp(-gamma1 * (x - xg)) * (np.cos(gamma1*(x - xg))*(alpha - np.tan(s)*xg / height_bb) - np.sin(gamma1*(x-xg))/(1+gamma1/gamma0) * ((gamma0/gamma1)**2 * alpha/kappa + gamma1/gamma0*(alpha - np.tan(s)*xg / height_bb)))
    y_ground = 1/2 - alpha/kappa - np.tan(s)*x/height_bb + np.exp(gamma0 * (x - xg)) * (np.cos(gamma0*(x - xg))*alpha/kappa + np.sin(gamma0*(x - xg))/(1 + gamma1/gamma0)*(alpha/kappa + (gamma1/gamma0)**3*(alpha - np.tan(s)*xg/height_bb)))

    y_sln = np.zeros(Nx)
    y_sln[xg_ind] = y_ground[xg_ind]
    y_sln[xf_ind] = y_float[xf_ind]

    lower = y_sln - h/2
    upper = y_sln + h/2

    return y_sln, lower, upper, xg