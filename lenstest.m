function lenstest
%% Prepare parameters
lambda= .532;   % wavelength [µm]
N = 1024; % Resolution
img = [];

f = 200e3; 

x = linspace(-10,10,N);
[X, Y] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);
        
        
img = fftshift(fft2(exp(((1j*2*pi/lambda)*(R.^2))/(2*f))));

imagesc(abs(ifft2(ifftshift(img))))

end

function H=transfer_function(z)
H = exp(1j*2*pi*(z/wavelength) * sqrt(1-(wavelength*freqx).^2 - ((wavelength*freqy).^2)+0j));
end