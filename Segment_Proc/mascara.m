function [X areaX caja]=mascara(radio,tam)
% MASCARA dibuja un c�rculo de cierto radio con unos en el centro de una matriz 
%
%Entradas
%radio: indica el radio del c�rculo de unos 
%tam: tama�o de la matriz donde se dibuja el c�rculo
%
%Salidas
%X: matriz devuelta donde los unos indican el interior del c�rculo y los
%ceros el exterior
%areaX: Area del circulo correspondiente a la suma de todos los elementos de la matriz 
%caja: bounding box sobre el c�rculo que le circunscribe 

radio2=radio^2;%prec�mputo del radio al cuadrado
%coordenadas del centro
centro=[(tam(2)+1)/2 (tam(1)+1)/2];

%trazo la malla de distancias al centro para cada pixel
[C D]=meshgrid((1-centro(1):tam(2)-centro(1)).^2,(1-centro(2):tam(1)-centro(2)).^2);

%valores interiores a la ROI -> 1, exteriores ->0
X=double(C+D<=radio2);

[A(:,2), A(:,1)]=find (X==1);
caja=zeros(2,2);
%Caja: valores l�mite de la ROI
for l=1:2
    caja(l,1)=min(A(:,l));
    caja(l,2)=max(A(:,l));
end
areaX=sum(sum(X));%area de la ROI
        