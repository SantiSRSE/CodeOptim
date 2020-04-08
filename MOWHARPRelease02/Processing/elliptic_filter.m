function [ Ell,Circ ] = elliptic_filter( NW, center, or, rad_min, r_cent, ecc )
%radios al cuadrado
NT=size(rad_min); cc=(NT(1:2)+1)/2; %rad_min=rad_min^2;
rad_max=rad_min/(1-ecc^2); 
[Y0,X0]=meshgrid(1:NW,1:NW);
X=repmat(X0,[1 1 NT(3:4)]); Y=repmat(Y0,[1 1 NT(3:4)]);
Xc=repmat(center{1},[NW NW 1 1]); Yc=repmat(center{2},[NW NW 1 1]);

Ell=((( ((X-Xc)*cos(or) + (Y-Yc)*sin(or)) .^2)./rad_max)+ ...
    (( ((Y-Yc)*cos(or) - (X-Xc)*sin(or)) .^2)./rad_min))<1;
Circ_BW=0.5;
Circ=((X-cc(1)).^2)+ ((Y-cc(2)) .^2)<(r_cent*(Circ_BW^2)); 
Ell=(Ell-Circ)>0.5;
end