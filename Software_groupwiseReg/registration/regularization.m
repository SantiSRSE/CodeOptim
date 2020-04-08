function term = regularization( ts,TAux,bound_box,term )
% REGULARIZATION set the spatio/temporal regularization terms for the
% transformation calculation
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
%ts: bspline variables struct 
%TAux: matrix that stores the transformation
%bound_box: ROI
%term: struct of the spatio/temporal regularization terms
%
%Output
%term: Actualizating struct of the spatio/temporal regularization terms

%Dimensions of the submatrix       
r{1}=bound_box(1,1):bound_box(1,2);
r{2}=bound_box(2,1):bound_box(2,2);
tam=size(ts.BB); 
L=tam(4);
N=tam(3);
TAux=permute(TAux,[2 1 3 4 5 6 7 8]);
tic
%initialiseo dtau/dx y dtau/dt
term.dtaux=zeros(tam(1),tam(2),N,L,2,'single');
term.dtaut=zeros(tam(1),tam(2),N,L,'single');
term.dtaux2=zeros(tam(1),tam(2),N,L,2,'single');
term.dtaut2=zeros(tam(1),tam(2),N,L,'single');
term.dtauxy=zeros(tam(1),tam(2),N,L,'single');


%temporal Derivativeda given by (2) ; 
term.dtaut(r{1},r{2},:,:)=sum(sum((TAux - circshift(TAux,[0 0 1 0 0 0])).*ts.BB(r{1},r{2},:,:,:,:),5),6);
term.dtaut2(r{1},r{2},:,:)=sum(sum((TAux - 2*circshift(TAux,[0 0 1 0 0 0]) + circshift(TAux,[0 0 2 0 0 0])).*ts.BB(r{1},r{2},:,:,:,:),5),6);


%espatial Derivative given by (2)
for l=1:L
    term.dtaux(r{1},r{2},:,:,l)=sum(sum(TAux.*ts.BB1(r{1},r{2},:,:,:,:,l),5),6);%dtaul/dx 
    term.dtaux(r{1},r{2},:,l,l)=term.dtaux(r{1},r{2},:,l,l)+1;%dtau1/dx1 delta added
    term.dtaux2(r{1},r{2},:,:,l)=sum(sum(TAux.*ts.BB2(r{1},r{2},:,:,:,:,l),5),6);%dtau2/dx2
end
term.dtauxy(r{1},r{2},:,:)=sum(sum(TAux.*ts.BB11(r{1},r{2},:,:,:,:),5),6);%dtau2/dxdy
toc
end