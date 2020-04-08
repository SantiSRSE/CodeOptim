function V = metric( IT,data,term,bound_box )
% METRICA metric calculation for a given image set
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
%Entradas
%IT: images for metric calculation
%data: struct that stores the metric 
%term: matriz that stores the regularization terms
%bound_box: ROI
%
%Salidas
%V: metric at each pixel
tic
tam=size(IT);
V=zeros(tam(1),tam(2));%initialisation

switch data.metric
    case 'variance'
        %variance of the images        
        media=mean(IT,3);       
        for k=1:tam(3)
             V=V+(IT(:,:,k)-media).^2;
        end
        V=V/tam(3);                
       
    case 'reference'
        V=mean((IT-term.ref).^2,3);
                           
    otherwise
        sprintf('error')

end

%smooth terms added
if sum(term.landa)~=0
    r2=bound_box(2,1):bound_box(2,2); r1=bound_box(1,1):bound_box(1,2);
    V(r2,r1)=V(r2,r1)+sum(sum( (term.landa(1)*sum((term.dtaux(r2,r1,:,:,:).^2),5))+(term.landa(2)*sum((term.dtaux2(r2,r1,:,:,:).^2),5))+(term.landa(2)*2*(term.dtauxy(r2,r1,:,:).^2))+ ...
        ((term.dtaut2(r2,r1,:,:).^2)*term.landa(4)) + ((term.dtaut(r2,r1,:,:).^2)*term.landa(3)) ,4),3);
end

toc