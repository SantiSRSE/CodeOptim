function [X, areaX, bound_box]=set_mask(radius,Im_size)
% set_mask draws a circle of a given radius over the ROI in the center of the image 
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
%radius: radiius of the circle
%Im_size: image size
%
%Output
%X: matrix of the circular mask
%areaX: total area of the ROI 
%bound_box: bounding box of X 

radius2=radius^2;%precomputing the squared radius
%center coordinates
center=[(Im_size(2)+1)/2 (Im_size(1)+1)/2];

%mesh of the distances to the center of each pixel
[C D]=meshgrid((1-center(1):Im_size(2)-center(1)).^2,(1-center(2):Im_size(1)-center(2)).^2);

%inner values of ROI -> 1, outter ->0
X=double(C+D<=radius2);

[A(:,2), A(:,1)]=find (X==1);
bound_box=zeros(2,2);
%bound_box: límits of the ROI
for l=1:2
    bound_box(l,1)=min(A(:,l));
    bound_box(l,2)=max(A(:,l));
end
areaX=sum(sum(X));%area of the ROI
        