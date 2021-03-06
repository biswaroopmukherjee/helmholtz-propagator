function zslice
 
%% Prepare parameters
lambda= .532;   % wavelength [�m]
N = 1024; % Resolution
type = 'axicon'
img = [];

for z=0:1:50
    
    %% Prepare transfer function
    freqx=-10:20/N:10-1/N;    % setup frequency axis [1/�m]
    [freqx,freqy] = meshgrid(freqx,freqx);
    H = exp(1i*2*pi*(z/lambda)*sqrt(1-(lambda*freqx).^2-(lambda*freqy).^2));

    %% Calculate Initial Field
    x=linspace(-10,10,N);
    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2 + Y.^2);
    theta = atan(Y./X);
    switch type
        case 'axicon'
            Field=exp(-1*(R).^2).*exp(-18*1i.*R);
        case 'spiralaxicon'
            Field=exp(-0.5*(R).^2).*exp(-8*1i.*R+6*1i.*theta);
        case 'pseudo-LG'
            Field=exp(-0.5*(R).^2).*exp(6*1i.*theta);
    end

    %% Fresnel Propagation
    u0=Field;    % setup aperture
    a0=(fftshift(fft2(u0)));    % fourier transform
    az=a0.*H;    % multiply with transfer function
    uz=ifft2(fftshift(az));    % inverse fourier transform
    p=uz.*conj(uz);
    img = [img ,p(:,N/2)];
end

%% Plot
figure(1)
imagesc(img)
