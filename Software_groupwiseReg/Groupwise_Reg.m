% Sample script to perform Groupwise registration on 2D MR-CINE sequences.
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


%For notation details see reference:
%L. Cordero-Grande, S. Merino-Caviedes, S. Aja-Fernandez, and C. Alberola-Lopez,
%"Groupwise elastic registration by a new sparsity-promoting metric: application 
%to the alignment of cardiac magnetic resonance perfusion images,"
%IEEE Trans. Pattern Anal. Mach. Intell., vol. in press, 2013.

addpath ./registration;
addpath ./IO;
addpath ./preprocessing;

route=[pwd '\images\']; 
%read_images
[I,Im_size,Mask,Tdias]=Im_read(route);%N, number of images

N=Im_size(3); %temporal phases
margin=ceil(Im_size/40);% 2.5%-margin for boundaries
[x(:,:,1), x(:,:,2)] =meshgrid(1:Im_size(2),1:Im_size(1));
%bounding box of the ROI
radius=75; %set the radius of the circular ROI
[X, areaX, bound_box]=set_mask(radius,Im_size);

%E spline order, Dp point density
Dp=[4 4]; E=3; % set bspline variables, stored in ts
tic, ts=creabspline(x,bound_box,Dp,E,N); toc
term.landa=[0 0.005 0 0]; % set spatio-temporal weigths [1st spat 2nd spat 1st temp 2nd temp]
data= struct('interpolation','linear','metric','variance');
%weight matrix for each point displacement
W=25*ones(sum(ts.C(2,:))+1,sum(ts.C(1,:))+1,N,ts.nt); % set by Armijo
Wn=W/areaX;%normalisation
flagW=1; %Flag for adaptative set of W matrix

%Defining parameteos for gradient descent optimisation
% nmax=60;
% et=0.01; 0.01 pixels 
% eh=0.005; a 0.5% of initial
parameters= struct('nmax',60,'et',0.01,'eh',0.005);

%initialisation for transformation
T=zeros(sum(ts.C(2,:))+1,sum(ts.C(1,:))+1,N,ts.nt);

%Groupwise Regsitration
[ H, tiempo, T, IT, xn] = optimizator( x,ts,term,T,I,data,X,parameters,Wn,bound_box,flagW );

%mask propagation
difxn=(xn-repmat(xn(:,:,Tdias,:),[1 1 N 1]));
Mfin = interpolator( x,permute(repmat(x,[1 1 1 N]),[1 2 4 3])-difxn,abs(repmat(Mask,[1 1 N])) ...
    ,'linear',[1+margin(2) Im_size(2)-margin(2); margin(1)+1 Im_size(1)-margin(1)]);

for i=1:N %plot results
    figure(1)
    [Lines,Vertices]=isocontour(double(Mfin(:,:,i)>0.5),0.5);
    imshow(I(:,:,i),[]), hold on; %original image
    V1=Vertices(Lines(:,1),:); V2=Vertices(Lines(:,2),:); 
    plot([V1(:,2) V2(:,2)]',[V1(:,1) V2(:,1)]','b'); %propagated masks
    % gif generation
    F=getframe(gcf);
    [X,map] = frame2im(F);
    gif(:,:,:,i)=X;
end

sec2gifRGB( double(gif(75:225,125:275,:,:)), 'Group_Reg')
