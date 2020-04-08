function y = b1spline3( x )
% B1SPLINE3 first derivative of third order B-spline%

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
y = zeros(size(x));

% not null range
x0range = abs(x) <= 1;
x1range = x >1 & x<2;
x2range = x <-1 & x>-2;

val=x(x0range);%precompute
% asigne values to range
y(x0range) = 1.5*sign(val).*val.^2 -2*val;
y(x1range) = -(x(x1range)-2).^2/2;
y(x2range) = (x(x2range)+2).^2/2;

end
