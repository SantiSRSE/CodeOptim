function [xn, TAux] = transformation( x,ts,T,bound_box,flagopt )
% TRANSFORMATION generate a new control point mesh from the original
% according to their displacements
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
%x: initial mesh 
%ts: bspline variables struct 
%T: transformation matrix  with  displacements
%bound_box: determine the region of interest
%flagopt: flag activating low memory method (active in zero)
%
%output
%xn: new mesh generated from transformation
%TAux: matrix with new values reshaped to fit 

%Value of N, according to transformation 
if  ndims(T)==2   
    N=size(T,2);
else
    N=size(T,3);
end
%inicialisation of xn
tama=size(x);
L=tama(3);
xn=zeros(tama(1),tama(2),N,L);
TAux=0;

switch ts.type
    case 't' %traslational
        %tic
        for k=1:N
            for l=1:2
                xn(:,:,k,l)=x(:,:,l)+T(l,k);%displacement of pixels 
            end
        end    
     case 'tyr' %rigid
        %c and x (coord xy) for precomputation
        resta(:,:,1)=x(:,:,1)-ts.c(1);
        resta(:,:,2)=x(:,:,2)-ts.c(2);
        for k=1:N
            %precomputing sinus cosines 
            coseno=cos(T(3,k));
            seno=-sin(T(3,k));
            xn(:,:,k,1)=coseno*resta(:,:,1)-seno*resta(:,:,2)+ts.c(1)+T(1,k); %desplazamiento de pixel en horizontal (l=1)
            xn(:,:,k,2)=seno*resta(:,:,1)+coseno*resta(:,:,2)+ts.c(2)+T(2,k);%y vertical (l=2)        
        end  
    case 'b'
        
        if (flagopt)% normal procedure
            r{1}=bound_box(1,1):bound_box(1,2);
            r{2}=bound_box(2,1):bound_box(2,2);
           
            tic
            %TAux: auxiliar variable of size ts.BB with the values of T
            TAux=zeros(size(r{2},2),size(r{1},2),N,2,ts.E+1,ts.E+1,'single');
            for k=0:ts.E
               for l=0:ts.E
                   TAux(:,:,:,:,l+1,k+1)=T(k+ts.coef{2}(r{2},1),l+ts.coef{1}(r{1},1),:,:);
               end
            end
            %obtaining xn
            xn(r{2},r{1},:,:)=sum(sum(TAux.*ts.BB(r{2},r{1},:,:,:,:),5),6);
            toc
        else %low memory procedure
            for m=1:L
                ts.coef{m}(ts.coef{m}(:,1)<0,1)=0; %Se modifican los coeficientes para eliminar efectos de bordes
                ts.coef{m}(ts.coef{m}(:,1)>ts.coef{m}(bound_box(m,2),1)+1,1)=ts.coef{m}(bound_box(m,2),1)+1;
                ts.coef{m}=ts.coef{m}+1;
            end
            tic
            ts.BB=permute(ts.BB,[6 5 3 4 2 1]); %resahping ts.BB          
            %out of rangeo,T is 0 for smooth transition
            TAux=zeros(size(T)+[2 2 0 0]);
            TAux(2:end-1,2:end-1,:,:)=T;
            Dp=ceil(ts.Dp);
            %alternative procedure for low memory
            for m=bound_box(1,1)-Dp(1):bound_box(1,2)+Dp(1) %margin of Dp
                for n=bound_box(2,1)-Dp(2):bound_box(2,2)+Dp(2)
                    xn(n,m,:,:)=sum(sum(ts.BB(:,:,:,:,m,n).*TAux(ts.coef{2}(n,1):ts.coef{2}(n,1)+ts.E,ts.coef{1}(m,1):ts.coef{1}(m,1)+ts.E,:,:),2),1);
                end
            end
            toc
       
        end
        xn=permute(repmat(x,[1 1 1 N]),[1 2 4 3])+xn;

    otherwise
        sprintf('error')
        
end
