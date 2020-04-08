function dtheta  = gradient_transformation( ts,x,T )
% GRADIENT_TRANSFORMATION calculate the gradient of the transformation
% traslational, rigid  or bspline
%
%	Copyright (C) 2016 Santiago Sanz Estébanez <ssanest@lpi.tel.uva.es>
%	Laboratorio de Procesado de Imagen, Universidad de Valladolid
%	www.lpi.tel.uva.es
%
%   When indicated, some parts are based on code from Lucilio Cordero Grande:
%       https://kclpure.kcl.ac.uk/portal/en/persons/lucilio-corderogrande.html
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%input
%ts transformation parameters (matrix bspline)
%x: intial mesh
%T: transformation matrix with displacements
%
%output
%dtheta: gradient of the transformation

%nº of images according to la trasnformation
if  length(size(T))==2   
    N=size(T,2);
else
    N=size(T,3);
end
tama=size(x);
L=2;
switch ts.type
    % derivatiive in each pixel for each image with respect to each parameter 
    case 't'
        %traslational
        dtheta= zeros(tama(1),tama(2),N,L,ts.nt);                 
        dtheta(:,:,:,1,1)=1; %dx1/dt1; dx1/dt2=0; dx2/dt1=0; 
        dtheta(:,:,:,2,2)=1; %dx2/dt2
        
    case 'tyr'
        %rgid
        dtheta= zeros(tama(1),tama(2),N,L,ts.nt);
        %Precómputing the difference
        resta(:,:,1)=x(:,:,1)-ts.c(1);
        resta(:,:,2)=x(:,:,2)-ts.c(2);
        
        dtheta(:,:,:,1,1)=1; %dx1/dt1; dx1/dt2=0; dx2/dt1=0
        dtheta(:,:,:,2,2)=1;%dx2/dt2
        for k=1:N
          % precómputing sines, cosines.
          seno=-sin(T(3,k));
          coseno=-cos(T(3,k));
          
          dtheta(:,:,k,1,3)= seno*resta(:,:,1)-coseno*resta(:,:,2);%dx1/dalfa
          dtheta(:,:,k,2,3)= coseno*resta(:,:,1)+seno*resta(:,:,2);%dx2/dtalfa
        end
                
    case 'b'
        %bsplines
        dtheta=permute(repmat(ts.BBg,[1 1 1 1 N]),[1 2 5 3 4]);
               
    otherwise 
       error('no implementada')
 
end

end