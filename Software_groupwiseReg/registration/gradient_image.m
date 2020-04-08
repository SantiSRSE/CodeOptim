function [ dx1, dx2 ] = gradient_image( IT,I,x,ts,T,caja,version2 )
% GRADIENT_IMAGE gradient (horizontal and vertical) of the images
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
%IT: images
%
%output
%dx1: gradient in horizontal direction 
%dx2: gradient in vertical direction 

L=ndims(I)-1;
%tic
if version2==0
    [dx1, dx2]=gradient(IT); % gradient of transformed images
else
    [d1, d2]=gradient(I); % transformed gradient of the images
    [xn, ~] = transformation( x,ts,T,caja,1 );
    dx1 = interpolator( x,xn,d1,'linear',caja);
    dx2 = interpolator( x,xn,d2,'linear',caja);
        
end


end



