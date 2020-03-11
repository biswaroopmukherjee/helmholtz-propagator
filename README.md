Helmholtz Propagator
============
Propagate a light field by solving the Helmholtz equation. 


## How does it work?
This is numerically achieved by applying the following transfer function:


![H_{12} = \exp \left(-i kz_{12} \sqrt{1-\lambda^2\left(k_x^2+k_y^2\right)}\right)](https://render.githubusercontent.com/render/math?math=H_%7B12%7D%20%3D%20%5Cexp%20%5Cleft(-i%20kz_%7B12%7D%20%5Csqrt%7B1-%5Clambda%5E2%5Cleft(k_x%5E2%2Bk_y%5E2%5Cright)%7D%5Cright))

in the Fourier domain, so that if the initial E-field is ![U_0(x,y)](https://render.githubusercontent.com/render/math?math=U_0(x%2Cy)) then


![U(z) = \mathfrak{F}^{-1}\left\{H_{12}\cdot \mathfrak{F}\{U_0\}\right\}](https://render.githubusercontent.com/render/math?math=U(z)%20%3D%20%5Cmathfrak%7BF%7D%5E%7B-1%7D%5Cleft%5C%7BH_%7B12%7D%5Ccdot%20%5Cmathfrak%7BF%7D%5C%7BU_0%5C%7D%5Cright%5C%7D)

Run the jupyter notebook to see in action.

## Sample propagations

We can plot the propagation of something simple like a Gaussian beam:

![gaussian](figures/gaussian.png)

or something more sophisticated like an axicon with a lens:

![axilens](figures/lens.png)