%imagen de partida 
A       =double(imread('cameraman.tif'));
%stripes en diferentes direcciones
lambda  =7.15;
x       =-127:128;
coef1   =[1 1 0 0];
coef2   =[1 0 0 0];
coef3   =[1 0 0 0];
coef4   =[1 0 0 0];
frec    =2*pi/lambda;
k       =(0:length(coef2)-1)';

%patron cartesiano
stripe1x = coef1* cos(frec*k*x);
stripe1y = (coef2* cos(frec*k*x))';
patron1  = repmat(stripe1x,length(x),1).*repmat(stripe1y,1,length(x));

%patron +/- 45
buffer45 = repmat(2*pi/lambda*(repmat(x*1/sqrt(2),length(x),1)+repmat(x'*1/sqrt(2),1,length(x))),[1 1 length(k)]);
bufferm45= repmat(2*pi/lambda*(repmat(x*1/sqrt(2),length(x),1)-repmat(x'*1/sqrt(2),1,length(x))),[1 1 length(k)]);
buffer2  = repmat(permute(k,[2 3 1]),[length(x) length(x) 1]);
coefvol3  = repmat(permute(coef3,[1 3 2]),[length(x) length(x) 1]);
coefvol4  = repmat(permute(coef4,[1 3 2]),[length(x) length(x) 1]);
patron2 = sum(cos(buffer45.*buffer2).*coefvol3,3).*sum(cos(bufferm45.*buffer2).*coefvol4,3);


figure(1),imshow(patron1,[])
figure(2),imshow(patron2,[])
figure(3),imshow(fftshift(fftshift(log10(abs(fft2(A.*(patron1.*patron2)))),2),1),[]);
