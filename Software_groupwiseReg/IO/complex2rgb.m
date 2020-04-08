function [rgb] = complex2rgb(A, scale)

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

if nargin == 1
    reduce = abs(A);
    ndim = length(size(reduce));
    for i=1:ndim
        reduce = max(reduce);
    end
    scale = reduce;
end

H = angle(A)+pi/2;
H (H<0)=H(H<0)+2*pi;
H = H/(pi/3);
S = 1;
V = abs(A)/scale;

%C = (1 - abs(2 * Y - 1)) * S/2.481;
C = V * S;
X = C .* (1 - abs(rem(H,2) - 1));
Z = zeros(size(V));

R = C;
G = X;
B = Z;

R(H>1) = X(H>1);
G(H>1) = C(H>1);
B(H>1) = Z(H>1);

R(H>2) = Z(H>2);
G(H>2) = C(H>2);
B(H>2) = X(H>2);

R(H>3) = Z(H>3);
G(H>3) = X(H>3);
B(H>3) = C(H>3);

R(H>4) = X(H>4);
G(H>4) = Z(H>4);
B(H>4) = C(H>4);

R(H>5) = C(H>5);
G(H>5) = Z(H>5);
B(H>5) = X(H>5);

%m = Y - (0.30*R+0.59*G+0.11*B);
m = V - C;
R = R + m;
G = G + m;
B = B + m;
rgb = reshape([R(:), G(:), B(:)], [size(R), 3]);
rgb(rgb<0)=0;
rgb(rgb>1)=1;
end
