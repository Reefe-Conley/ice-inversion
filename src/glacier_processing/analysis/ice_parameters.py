from dataclasses import dataclass

@dataclass(frozen=True)
class DimensionalParameters:
    """Dimensional inputs (SI units) for the ice beam, bedrock and ocean."""

    E: float        # [Pa]      Young's modulus.  
    nu: float       # [-]       Poisson's ratio.
    H: float        # [m]       Characteristic (reference) ice beam thickness.
    rho_i: float    # [kg/m3]   Ice density.
    rho_w: float    # [kg/m3]   Water density.
    k: float        # [Pa/m]    Bed stiffness.
    g: float        # [m/s2]    Acceleration due to gravity.


@dataclass(frozen=True)
class NondimensionalParameters:
    """Nondimensional inputs for the ice beam, bedrock and ocean."""

    alpha: float    # [-]       Density ratio rho_i/rho_w
    kappa: float    # [-]       Force ratio k/rho_w*g
    H: float        # [m]       Characteristic (reference) ice beam thickness.
    L: float        # [m]       Bending-buoyancy length scale

def nondimensionalise(material: DimensionalParameters) -> NondimensionalParameters:
    """
    Return the nondimensional parameters alpha and kappa along with the vertical
    and horizontal length scales H and L.
    """
    alpha = material.rho_i / material.rho_w
    kappa = material.k / (material.rho_w * material.g)
    H = material.H
    L = (material.E * material.H**3 / (12 * (1 - material.nu**2) * material.rho_w * material.g))**(1/4) # Bending-buoyancy length scale.

    return NondimensionalParameters(
        alpha=alpha,
        kappa=kappa,
        H=H,
        L=L
    )