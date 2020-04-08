function [ Tdif, Hdif, Wn] = evolution( T,H,iter,Wn,flagW )
% EVOLUTION calculate the differences between inputs (norm and y metric)
% and recast them to adequeate them to output conditions
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
%T: transformation matrix (displacement)
%H:values of the metric for each iteration
%iter: current iteration 
%Wn: weight matrix
%flagW: flag for adative W algorithmn
%
%Output
%Tdif: Norm difference of the transformation matrix
%Hdif: Metric variation
%Wn: new weight matrix

%calculus of the differences
%norm
if  ndims(T)==2   %traslational
    if size(T,1)==2
        Tdif=norm(T);
    elseif  size(T,1)==3  %rígid
        %ponderation 
        Tdif= norm(1./(Wn*(sum(1./(Wn(:,1))))).*T);
    end
elseif ndims(T)==4 % b-splines
    Tdif=norm(T(:));
end

%metric
Hdif=H(1,iter)-H(1,iter+1);
% adaptation of Wn, flagW
if flagW
    if Hdif>0
        Wn=Wn*1.2;
    else
        Wn=Wn/2;
    end
end
