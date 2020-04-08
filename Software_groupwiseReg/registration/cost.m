function cost_total = cost( V,mask,ts )
% COST determines the total cost of a given metric array 
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
%V: metric array
%mask: determines the subset to take into account
%ts: bsplines variables  struct, only in 6 %dimensional case
%
%Output
%cost_total:total cost for each frame with 1,2 o 4 dimensiones 

%Avoid NaN 
V(isnan(V))=0;

if ndims(V)==2
    %Sum of the metric values only in the ROI
    cost_total=sum(sum(V.*mask));
elseif ndims(V)==4   
    %dH (cost) traslational and rigid
    %Sum of the metric values only in the ROI for each image and dimension
    Aux=permute(sum(sum(bsxfun(@times,V,mask),2),1),[3 4 1 2]);
    cost_total=Aux(:,:,1,1);        
elseif ndims(V)==6
    tam=size(V);
    cost_total=zeros(tam(5),tam(4),tam(3),tam(6)); %init. of gradient costs
    V=permute(V,[1 2 4 5 3 6]);% dimensional reshaping; dH is [py px N L] and dV is [y x N px py L]
        
    %Sum of the metric values only in the ROI for each image and point of the mesh
    for k=1:tam(5) %py
        for l=1:tam(4) %px
            %auxiliar Matrix
            Aux=repmat(mask(ts.coefg{2}(k,1):ts.coefg{2}(k,2),ts.coefg{1}(l,1):ts.coefg{1}(l,2)),[1 1 1 1 tam(3) tam(6)]);
            
            cost_total(k,l,:,:)=sum(sum(V(1:ts.coefg{2}(k,2)-ts.coefg{2}(k,1)+1,1:ts.coefg{1}(l,2)-ts.coefg{1}(l,1)+1,l,k,:,:).*Aux,2),1);
        end
    end    
end