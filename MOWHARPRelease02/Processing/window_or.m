function ret = window_or(N, sigma, theta)
%only squared windows 

ret  = zeros(N);
theta   = (theta/180)*pi; %angulo

rbeg  = -round(N / 2); %valor medio

if mod(N,2)==1
    center=0; %Ventana impar maximo en medio
else
    center=1/2; %Ventana par
end

% [xm, ym]=meshgrid(1:N,1:N);
% ret2=exp(-((((rbegin+xm-center)*cos(theta) + (cbegin+ym-center)*sin(theta))/sigma).^2)/2);

%gaussiana girada y extendida (customgauss)
%[c,r]=meshgrid(1:N);
for r=1:N
    for c=1:N
        ret(r,c)=exp(-((((rbeg+r-center)*cos(theta) - (rbeg+c-center)*sin(theta))/sigma)^2)/2);
    end
end

%equivalente a gausswin(N,(N-1)/(2*sigma))
%if theta==0, plot(1:N,gausswin(N,(N-1)/(2*sigma)),'x',1:N,ret(:,1),'r')
%elseif theta==pi/2, plot(1:N,gausswin(N,(N-1)/(2*sigma)),'x',1:N,ret(1,:),'r'), end

% xm      = (x-xc)*cos(theta) - (y-yc)*sin(theta);
% ym      = (x-xc)*sin(theta) + (y-yc)*cos(theta);
