function term  = gradient_regularization( ts,term )
% GRADIENT_REGULARIZATION recast the spatio-temporal smooth terms for the 
% calculation of the gradient 
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
%term: initial struct where the spatio-tempora terms are store
%
%Output
%term: struct where the spatio-temporal terms are stored (needed for the
%gradient)

%Dimensions
N=size(term.dtaux,3);
L=2;
lim=size(ts.BBg);
%initialisation
term.dthetax=zeros([lim(1:2) N lim(3:4) 2 2],'single');
term.dthetax2=zeros([lim(1:2) N lim(3:4) 2 2],'single');
term.dthetaxy=zeros([lim(1:2) N lim(3:4) 2],'single');
term.dthetat=zeros([lim(1:2) N lim(3:4) 2],'single');
term.dthetat2=zeros([lim(1:2) N lim(3:4) 2],'single');
tic
dtaux(:,:,:,1,1,:,:)=term.dtaux;%dimensional reshaping
dtaux2(:,:,:,1,1,:,:)=term.dtaux2;
dtauxy(:,:,:,1,1,:)=term.dtauxy;
clear term.dtaux

dim=zeros(2,1);
for k=1:lim(3)
   for l=1:lim(4)
     %not null ranges          
     range{1}=ts.coefg{1}(k,1):ts.coefg{1}(k,2);
     range{2}=ts.coefg{2}(l,1):ts.coefg{2}(l,2);

     for m=1:ts.nt
       dim(m)=size(range{m},2);
     end
          
     %Spatial Derivative according to (8)
     for m=1:ts.nt
       term.dthetax(1:dim(2),1:dim(1),:,k,l,:,m)=dtaux(range{2},range{1},:,1,1,:,m).*ts.BB1g(1:dim(2),1:dim(1),:,k,l,:,m);%d/dx1
       term.dthetax2(1:dim(2),1:dim(1),:,k,l,:,m)=dtaux2(range{2},range{1},:,1,1,:,m).*ts.BB2g(1:dim(2),1:dim(1),:,k,l,:,m);%d/dx2
     end      
     term.dthetaxy(1:dim(2),1:dim(1),:,k,l,:)=dtauxy(range{2},range{1},:,1,1,:).*ts.BB11g(1:dim(2),1:dim(1),:,k,l,:);%d/dxdy
     
     %Precomputing ts.BBg out of the loop
     
     Aux=repmat(ts.BBg(1:dim(2),1:dim(1),k,l),[1 1 N 1]);
     %temporal derivative according to (10)
     for q=1:L
        term.dthetat(1:dim(2),1:dim(1),:,k,l,q)=2*(term.dtaut(range{2},range{1},:,q)-circshift(term.dtaut(range{2},range{1},:,q),[0 0 -1 0])).*Aux;
        term.dthetat2(1:dim(2),1:dim(1),:,k,l,q)=2*(term.dtaut2(range{2},range{1},:,q)-2*circshift(term.dtaut2(range{2},range{1},:,q),[0 0 -1 0]) + ...
            circshift(term.dtaut2(range{2},range{1},:,q),[0 0 -2 0])).*Aux;
     end
     
   end
end
clear term.dtaut term.dtaut2 dtaux dtaux2 dtauxy Aux

toc
