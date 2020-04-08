function IT = interpolator( x,xn,I,interpolation,bound_box)
% INTERPOLATOR interpolate a mesh over the images ( original mesh) to 
% generate new imáges corresponding to the new mesh
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
%x: initial mesh
%xn:final mesh
%I: imáges 
%interpolatión: method of interpolation; lineal, cubic,...
%bound_box: determine the ROI
%
%Output
%IT: transformed images accordign to xn

tam=size(I);
if ndims(I)==2
    tam=[tam 1];
end
IT=I;% transformed images
margin=ceil(tam/40);% margin
r{1}=bound_box(2,1)-margin(1):bound_box(2,2)+margin(1);%margin accorging to image size
r{2}=bound_box(1,1)-margin(2):bound_box(1,2)+margin(2);
for k=1:tam(3)
    IT(r{1},r{2},k)=interp2(x(r{1},r{2},1),x(r{1},r{2},2),I(r{1},r{2},k),xn(r{1},r{2},k,1),xn(r{1},r{2},k,2),interpolation);
    %lineal interpolation  (default), careful with NaN     
    
end


