function xr=regularizeBSpline(x,tS,Field)
% REGULARIZAR actualiza los t�rminos de suavidad espacial y temporal
% envueltos en el c�lculo de la transformaci�n 
%
%Entradas
%ts: estructura de almacenamiento de las variables bspline
%TAux: matriz de un tama�o conveniente que almacena los valores de la
%matriz de la transformaci�n
%caja: determina la zona de la imagen donde se opera
%term: estructura de almacenamiento de las variables de suavidad temporal y
%espacial
%
%Salidas
%term: estructura de almacenamiento de las variables de suavidad temporal y
%espacial, con los nuevos valores necesarios para la transformaci�n

%NOTE NOW DERIVATIVES ARE ALONG 4TH AND 9TH DIMENSIONS

IL=eye(L,'like',x);
IL=permute(IL,[3 4 5 1 6 7 8 9 2]);

%Derivada espacial dada en (2)
if Field.lambda(1)>0
    xr.dr1=extractROI(bsxfun(@plus,multDimSum(bsxfun(@times,x,extractROI(ts.BB1,ts.ROI,1)),6:8),IL),ts.ROI,0);%NOT SURE IF IL SHOULD BE ADDED BEFORE OR AFTER FILLING THE ROI. DOES NOT SEEM CRITICAL AT THE MOMENT
end
if Field.lambda(2)>0
    xr.dr2=extractROI(multDimSum(bsxfun(@times,x,extractROI(tS.BB2,ts.ROI,1)),6:8),tS.ROI,0);%d2taul/dx2
    xr.dr11=extractROI(2*multDimSum(bsxfun(@times,x,extractROI(tS.BB11,ts.ROI,1)),6:8),tS.ROI,0);%d2taul/dxidxj
end


