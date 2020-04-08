function y = bspline2( x )
% BSPLINE2 B-spline order 2

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
%precompute abs

val=abs(x);
%not null range
x0range = val < 0.5;
x1range = val >=0.5 & val<1.5;

% asigne values to range
y(x0range) = -x(x0range).^2 + 0.75;
y(x1range) = 0.5*(val(x1range)-1.5).^2;
end