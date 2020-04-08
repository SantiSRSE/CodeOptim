function y = bspline3( x )
% BSPLINE3 B-spline order 3

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

%precompute abs
val=abs(x);

y = zeros(size(x));
% not null range
x0range = val < 1;
x1range = val >=1 & val<2;

% asigne values to range
y(x0range) = 0.6666 - val(x0range).^2 + 0.5*val(x0range).^3;
y(x1range) = (2-val(x1range)).*(2-val(x1range)).*(2-val(x1range))/6;
end