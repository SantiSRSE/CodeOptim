function dy = gradient_metric( IT,data,term,bound_box )
% GRADIENTE_METRICA calculates the gradient of the metric over the images

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
%data: structure where the metric is defined
%
%output
%dy: gradient of the metric 
tic
tam=size(IT); dy=zeros(tam);
switch data.metric
    case 'variance'
        %calculation of the derivative of the image intensity variance dV/dy
        dy=(2/tam(3)).*bsxfun(@minus,IT,mean(IT,3));
    case 'reference'
        dy=2*(IT-term.ref);
    otherwise
        sprintf('error')

end

toc
