function y = b1spline2( x )
% B1SPLINE2 first derivative of second order B-spline%

%	Copyright (C) 2016 Santiago Sanz Est�banez <ssanest@lpi.tel.uva.es>
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
y = zeros(size(x));
%prec�mpute abs
val=abs(x);
% not null range
x0range = val < 0.5;
x1range = val >=0.5 & val<1.5;

val=x(x1range);%precompute
% asigne value to range 
y(x0range) = -2*x(x0range);
y(x1range) = val-1.5*sign(val);

end