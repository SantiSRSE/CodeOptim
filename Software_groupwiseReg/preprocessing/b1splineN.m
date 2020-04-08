function y=b1splineN(x,n)
% B1SPLINEN first derivative of order n B-spline 
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
%
%input
%x: values
%n: bsplineorder 
%
%output
%y: valoues of bspline 

%selection
if(n==1)
    y = b1spline1(x);
elseif(n==2)
    y = b1spline2(x);
elseif(n==3)
    y = b1spline3(x);
end

end