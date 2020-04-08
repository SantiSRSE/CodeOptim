function dV = gradient_joint( IT,I,ts, data,x,T,term,bound_box )
% GRADIENTE_JOINT assemble each part of the gradiente calculation: image
% gradient, metric gradient and transformation gradiente with smoothness
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
%Input
%x: initial mesh grid
%IT: images 
%ts: bspline variables structure
%data: type of transformation and metric employed
%T: transformation matrix with displacements
%term: smoothness structure
%
%Ouput
%dV: gradient obtained using gradient descent method

tam=size(IT);
L=2;
%gradient independently form interpolation
[dx(:,:,:,1), dx(:,:,:,2)] = gradient_image(IT,I,x,ts,T,bound_box,1);
dy = gradient_metric( IT,data,term,bound_box );
dtheta = gradient_transformation( ts,x,T);

if ndims(T)==2 %rigid
    if ndims(IT)==2, dV=zeros([tam 1 ts.nt]); end
    
    for k=1:ts.nt
       %gradient at each pixel for each parameter: (t1,t2,alfa)    
       dV(:,:,:,k)=dy.*(dx(:,:,:,1).*dtheta(:,:,:,1,k)+dx(:,:,:,2).*dtheta(:,:,:,2,k));
       %dV/Dyn(dyn/dx1*1+dyn/dx2*0)
       %dV/Dyn(dyn/dx1*0+dyn/dx2*1)
       %dV/Dyn(dyn/dx1* -sin(alfa)*(x-c1)+cos(alfa)(y-c2) + dyn/dx2*-sin(alfa)*(x-c1)-cos(alfa)(y-c2) )
    end 

elseif ndims(T)==4 %bsplines
    lim=size(ts.BBg);
    term  = gradient_regularization( ts,term );% smoothness
    %Precomputing the regularization term
    Aux=repmat(sum(sum(term.landa(1)*sum(term.dthetax,7)+ term.landa(2)*sum(term.dthetax2,7) + term.landa(2)*2*term.dthetaxy + ...
        term.landa(3)*term.dthetat + term.landa(4)*term.dthetat2 ,6),3),[1 1 tam(3) 1 1]);
    clear term.dthetax term.dthetat term.dthetat2 term.dthetax2 term.dthetaxy
    
    %initialisation
    dV=zeros(lim(1),lim(2),tam(3),lim(3),lim(4),ts.nt);
    dim=zeros(L,1);        
    tic
    for k=1:lim(3)
        for l=1:lim(4)
            %not null range            
            range{1}=ts.coefg{1}(k,1):ts.coefg{1}(k,2);
            range{2}=ts.coefg{2}(l,1):ts.coefg{2}(l,2);

            for m=1:ts.nt
              dim(m)=size(range{m},2);
            end
                        
            %Joint Implementation for the initial gradient calculus with
            %the smoothness terms
            for m=1:ts.nt
               dV(1:dim(2),1:dim(1),:,k,l,m)=dy(range{2},range{1},:).*dx(range{2},range{1},:,m).*dtheta(1:dim(2),1:dim(1),:,k,l)+Aux(1:dim(2),1:dim(1),:,k,l);
            end            
        end
    end
    toc
    clear dx dy dtheta Aux
end
