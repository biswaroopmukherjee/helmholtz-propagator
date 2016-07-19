x=linspace(-10,10,300);
[X,Y] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
z=exp(-5*(R-2).^2);
% z_p=z.*exp(1i*R*0.01*2*pi/532e-9);
z_p=exp(-1i*R).*z;
imagesc(z)
axis image
% figure;
% imagesc(abs(fresnel_advance(z_p,1e-6,1e-6,-0.00005,532e-9)));
% axis image
