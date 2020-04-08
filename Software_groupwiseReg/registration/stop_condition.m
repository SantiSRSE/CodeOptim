function  fin  = stop_condition( p1,p3,p2,p4,p5,p6 )
% stop_condition condition for the exit of the loop 
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
%p1,p3,p2,p4,p5,p6: values to compare
%
%output
%fin: boleean indicating the exit of the loop


%Retrun 1 if all conditions are true
fin=0;
%compare to achieve convergence
if (p5==p6)|| ((p1<p2)&&(p3<p4))
    fin=1;
end

end