%OPTIMIZATOR  principal loop of the registration algorithmno
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
%%%%%%%%%%%%%%%%%%Entradas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gradient parameters
% nmax=60;
% et=0.01;%0.01 pixels.
% eh=0.005;% a 0.5% of initial

%data struct for optimisation
% interpolation='linear';
% metric='varianza';

% transformation; struct ts

%x original meshgrid 
%T transformation matrix 
%I original images 
%X  mask that defines the ROI
%ts transformation parameters (matrix bspline)
%term smooth terms; weights and temporal and espatialderivatives 
%Wn weight matrix
%flagW flag for adaptive Wn

%%%%%%%%%%%%%%%%%%%%%%Salidas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H metric value at each iteration
%timing total time at each iteration
%T final transformation matrix
%IT transformed images 
%xn transformed mesh grid

function [ H, timing, T, IT, xn] = optimizator( x,ts,term,T,I,data,X,parameters,Wn,bound_box,flagW )
tam=size(I);%initial size
L=2;
if ndims(I)==3 
    N=tam(3);
    
    %initialisation time , cost
    if strcmp(ts.type,'b')
        term.dtaux=zeros(size(ts.BB,2),size(ts.BB,1),N,L,2,'single');%term. smooth
        term.dtaut=zeros(size(ts.BB,2),size(ts.BB,1),N,L,'single');
        term.dtaux2=zeros(tam(1),tam(2),N,L,2,'single');
        term.dtaut2=zeros(tam(1),tam(2),N,L,'single');
        term.dtauxy=zeros(tam(1),tam(2),N,L,'single');
    end
else
    N=1;
    I(:,:,1)=I;
end

H(1)=cost(metric(I,data,term,bound_box),X);
timing(1)=0;
IT=I;%inital images

p2=parameters.et*N*ts.nt;% parameters of stop (conditions)
p4= H(1)*parameters.eh;
iter=1;

while (1)
    tic 
    %Gradients
    clear xn
    dV= gradient_joint(IT,I,ts,data,x,T,term,bound_box);
    clear IT    
    dH=cost(dV,X,ts);
    clear dV
    %Proyection of gradient given by (2) en [2]:
    if ndims(I)>2, proydH = proyection( dH ); else proydH=dH'; end
    clear dH
    Dif=Wn.*proydH;
    T=T-Dif;% New T matrix
    clear proydH
       
    %transformtion of images
    [xn, TAux] = transformation( x,ts,T,bound_box,1 );
    %new regularization terms 
    if ndims(I)>2 
        if strcmp(ts.type,'b')
            term = regularization( ts,TAux,bound_box,term );
            clear TAux;
        end
    end
    
    %interpolation of original images 
    IT = interpolator( x,xn,I,data.interpolation,bound_box);
        
    %metric cost (1) with snmooth terms
    H(iter+1)=cost(metric(IT,data,term,bound_box),X);
            
    timing(iter+1)=toc;% time of current iteration 
    timing(iter+1)=timing(iter+1)+timing(iter);%time total
    
    %out condition
    [Tdif(iter), Hdif(iter), Wn]=evolution(Dif,H,iter,Wn,flagW);
    clear Dif   
    %Info 
    sprintf('Norm Variation in last iteration: %0.5f,threshold of %0.5f', Tdif(iter),p2)
    sprintf('Metric Variation in last iteration: %d,threshold of %d', Hdif(iter), p4)
    sprintf('Current value of the metric: %d', H(iter+1))
    sprintf('Iteration nº %d de %d', iter, parameters.nmax)
    sprintf('New weighting matrix: %0.5f', sum(sum(X))*Wn(1))
        
    %Compare: Tdif -> p2 , Hdif -> p4 , nmax -> iter
    if stop_condition( Tdif(iter),Hdif(iter),p2,p4,parameters.nmax,iter ) 
        break;
    end
          
    iter=iter+1;
end
term.landa=[0 0 0 0];
sprintf('original final Metric (landa=0): %d', cost(metric(IT,data,term,bound_box),X)) %original final metric 