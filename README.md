# Ice inversion
An inversion scheme to analyse ice sheets using a modified elastic beam model. 

## Calculations

We model the ice sheet in two dimensions. We asume that the static shape and thickness of the ice is determined by the balance between its weight, any elastic bending stresses within the ice, and any normal stresses applied to the ice. These normal forces are taken to be elastic compression and buoyancy where the ice is underlain by bedrock and ocean water, respectively.

The full force balance used in the inversion scheme is therefore

```math
\frac{d^2}{dx^2}\left[B(x) \frac{d^2z_m}{dx^2} \right]+\rho_i g h = k\left(\frac{h}{2} - z_m - b_0 \right)\mathcal{H}(h/2 - z_m - b_0) + \rho_w g \left(\frac{h}{2} - z_m \right)\mathcal{H}(h/2 - z_m),
```
where $x$ is the horizontal downstream coordinate, $h$ is the ice thickness, $z_m$ is the ice midline position, $k$ is the stiffness of the bedrock (modelled here as a Winkler foundation), $\rho_{i,w}$ are the densities of ice and water, $b_0$ is the uncompressed bedrock position, $g$ is the acceleration due to gravity, and $B(x)$ is the bending stiffness, given by
```math
B(x) = \frac{E[h(x)]^3}{12(1-\nu^2)},
```
where $E$ and $\nu$ are the Young's modulus and Poisson's ratio of ice (here assumed to be uniform constants throughout the beam). The Heaviside step function
```math
\mathcal{H}(x) =
\begin{cases}
1 & \text{if } x \geq 0, \\
0 & \text{if } x < 0.
\end{cases}
```
applies the appropriate normal reaction force from the bedrock and ocean along the beam.

