ii=sqrt(-1); 
lambda= .365;   % wavelenght [µm]
z=1;    % distance [µm]
w=1;    % width of slit is 2*w [µm]

x=-12.75:0.05:12.8;    % setup spatial axis [µm]
freqx=-10:20/512:10-1/512;    % setup frequency axis [1/µm]
freqy=freqx;


u0=zeros(512);    % field at z=0
a0=zeros(512);    % angular spectrum at z=0
H=zeros(512);    % transfer function
az=zeros(512);    % angular spectrum at z=z
uz=zeros(512);    % field at z=z

for nx=1:512    % setup transfer function
   for ny=1:512
     H(nx,ny)=exp(ii*2*pi*(z/lambda)*...
     sqrt(1-(lambda*freqx(nx))^2-(lambda*freqy(ny))^2));
   end
end

u0(257-w*20:256+w*20,257-w*20:256+w*20)=1;    % setup aperture
a0=(fftshift(fft2(u0)));    % fourier transform
az=a0.*H;    % multiply with transfer function
uz=ifft2(fftshift(az));    % inverse fourier transform
p=uz.*conj(uz);

figure(1)
plot(x, p(:,256));    % plot of cross-section of intensity at z
xlabel('x'); ylabel('I');

figure(2)
imagesc(x, x, p);    %diffraction pattern at z
xlabel('x'); ylabel('y'); colormap(gray);
colorbar;