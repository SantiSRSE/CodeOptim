function ts = creabspline(x,bound_box,Dp,E,N)
% CREABSPLINE store the neccesary elements for the gradient and bspline
% transformation calculus, B-spline product matrixes with the correspondent
% coefficients
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
%
%Input
%x: initial mesh
%bound_box: Region of interest
%Dp: Density of control points
%E: bspline order
%N: number of images
%
%Output
%ts: struct where the spline variables are stored

%Definition of parameters (coord. xy)
tamxy=[size(x,2) size(x,1) size(x,2) size(x,1)];%size

L=size(x,3);%Dimension number
ts= struct('type','b','nt',2,'Dp',Dp, 'E',E,'c',sum(bound_box(1:2,1:2),2)/2);%(1+tamxy(1:2))/2);
ts.rp=(ts.E+1)*ts.Dp/2;%radius of influence

for l=1:L%dimensions
    for b=1:2%min-max
        ts.C(l,b)=floor((abs(ts.c(l)-bound_box(l,b))+ts.rp(l))/ts.Dp(l));
    end
end    
% reordering each control point
[ts.u(:,:,1), ts.u(:,:,2)] =meshgrid(-ts.C(1,1):ts.C(1,2),-ts.C(2,1):ts.C(2,2));
for l=1:L
   ts.pu(:,:,l)=ts.c(l)+ts.Dp(l).*ts.u(:,:,l); %location of each control point in the mesh
end


%coefficients and product matrix for transformation
for l=1:L
   coef{l}=zeros(tamxy(l),2);%original coefficients  without offset
end

for b=1:2%min-max
   if rem(ts.E,2)~=0%order 2 
      coef{1}(:,b)=floor((x(1,:,1)-ts.c(1))/ts.Dp(1))+((ts.E+(-1)^b)/2)*(-1)^b;
      coef{2}(:,b)=floor((x(:,1,2)-ts.c(2))/ts.Dp(2))+((ts.E+(-1)^b)/2)*(-1)^b;
   elseif rem(ts.E,2)==0%order 1 y 3
      coef{1}(:,b)=round((x(1,:,1)-ts.c(1))/ts.Dp(1))+(ts.E/2)*(-1)^b;
      coef{2}(:,b)=round((x(:,1,2)-ts.c(2))/ts.Dp(2))+(ts.E/2)*(-1)^b;
   end
end

%initialisation of auxiliar matrix (2->derivative) 
for l=1:L
   BBAux{l}=zeros(tamxy(l), ts.E+1,3,'single');
end    
   
%Auxiliar matrix for the coefficients
for l=1:L
    xbis=1:tamxy(l);% image size 
    for n=1:tamxy(l) %margin for the image for the transformation
       for k=coef{l}(n,1):coef{l}(n,2)
           %Matrix deriv.1 in 3y4
           BBAux{l}(n,k-coef{l}(n,1)+1,1)=bsplineN(((xbis(n)-ts.c(l))/ts.Dp(l))-k,ts.E);
           BBAux{l}(n,k-coef{l}(n,1)+1,2)=b1splineN(((xbis(n)-ts.c(l))/ts.Dp(l))-k,ts.E);
           BBAux{l}(n,k-coef{l}(n,1)+1,3)=b2splineN(((xbis(n)-ts.c(l))/ts.Dp(l))-k,ts.E);
        end
    end
end 

BB=zeros(ts.E+1,ts.E+1,tamxy(1),tamxy(2),'single');%matrix of bspline products
BB1=zeros(ts.E+1,ts.E+1,tamxy(1),tamxy(2),2,'single');%matrix of bspline products Bprimax*By y Bprimay*Bx
BB2=zeros(ts.E+1,ts.E+1,tamxy(1),tamxy(2),2,'single');%matrix of bspline products B2primax*By y B2primay*B
BB11=zeros(ts.E+1,ts.E+1,tamxy(1),tamxy(2),'single');%matrix of bspline products Bprimax*Bprimay
for m=1:tamxy(1)%margin for the image 
     for n=1:tamxy(2)
         BB(:,:,m,n)=BBAux{2}(n,:,1)'*BBAux{1}(m,:,1);
         %matricial product 
         BB1(:,:,m,n,1)=BBAux{2}(n,:,1)'*BBAux{1}(m,:,2);%border effects.
         BB1(:,:,m,n,2)=BBAux{2}(n,:,2)'*BBAux{1}(m,:,1);
         BB2(:,:,m,n,1)=BBAux{2}(n,:,1)'*BBAux{1}(m,:,3);
         BB2(:,:,m,n,2)=BBAux{2}(n,:,3)'*BBAux{1}(m,:,1);
         BB11(:,:,m,n)=BBAux{2}(n,:,2)'*BBAux{1}(m,:,2);
     end
end

%offset added to coefficients, used in transformation 
for l=1:L
   ts.coef{l}=coef{l}+1-coef{l}(bound_box(l,1),1);% coefficients with offset
end

%implementation in 6D, great memory costs 
for l=1:L
    ts.BB1(:,:,:,:,:,:,l)=ts.Dp(l)*repmat(permute(BB1(:,:,:,:,l),[4 3 5 6 2 1]),[1 1 N L 1 1 1]);
    ts.BB2(:,:,:,:,:,:,l)=(ts.Dp(l)^2)*repmat(permute(BB2(:,:,:,:,l),[4 3 5 6 2 1]),[1 1 N L 1 1 1]);
end

ts.BB11=ts.Dp(l)*ts.Dp(2)*repmat(permute(BB11,[4 3 5 6 2 1]),[1 1 N L 1 1 1]);
ts.BB=repmat(permute(BB,[4 3 5 6 2 1]),[1 1 N L 1 1]);
clear BB BB1 BB2 BB11 BBAux





%coefficients and products matrix  for the gradient calculus
limxy=[size(ts.pu,2) size(ts.pu,1)] ;

%separate coefficients method
for l=1:L %dimensions
    ts.coefg{l}=zeros(limxy(l),2); 
end

for b=0:1%min-max
    ts.coefg{1}(:,b+1)=ceil(ts.pu(1,:,1)-ceil(ts.rp(1)))+b*2*ceil(ts.rp(1)); 
    ts.coefg{2}(:,b+1)=ceil(ts.pu(:,1,2)-ceil(ts.rp(2)))+b*2*ceil(ts.rp(2));     
end

%same ts.pu for BBg, BB1g
xbis2{1}=ts.pu(1,:,1)';   
xbis2{2}=ts.pu(:,1,2);    

for l=1:L
    BBAux{l}=zeros(limxy(l),2*ceil(ts.rp(l))+1,2,'single');
end

%auxiliar matrix for coefficients in each dimension ( xbis2)
for l=1:L    
    for n=1:size(xbis2{l}) %  l=1 corresponds withlim(2) and l=2 with lim(1)
        for k=ts.coefg{l}(n,1):ts.coefg{l}(n,2)
            BBAux{l}(n,k-ts.coefg{l}(n,1)+1,1)=bsplineN(((k-xbis2{l}(n))/ts.Dp(l)),ts.E);
            BBAux{l}(n,k-ts.coefg{l}(n,1)+1,2)=b1splineN(((k-xbis2{l}(n))/ts.Dp(l)),ts.E);
            BBAux{l}(n,k-ts.coefg{l}(n,1)+1,3)=b2splineN(((k-xbis2{l}(n))/ts.Dp(l)),ts.E);
        end
    end
end  
% bspline products  matriz 
ts.BBg=zeros(2*ceil(ts.rp(2))+1,2*ceil(ts.rp(1))+1,limxy(1),limxy(2),'single');
ts.BB1g=zeros(2*ceil(ts.rp(2))+1,2*ceil(ts.rp(1))+1,limxy(1),limxy(2),2,'single');
ts.BB11g=zeros(2*ceil(ts.rp(2))+1,2*ceil(ts.rp(1))+1,limxy(1),limxy(2),'single');
ts.BB2g=zeros(2*ceil(ts.rp(2))+1,2*ceil(ts.rp(1))+1,limxy(1),limxy(2),2,'single');
for m=1:limxy(2)
    for n=1:limxy(1)
        ts.BBg(:,:,n,m)=BBAux{2}(m,:,1)'*BBAux{1}(n,:,1); 
        ts.BB1g(:,:,n,m,1)=BBAux{2}(m,:,1)'*BBAux{1}(n,:,2); 
        ts.BB1g(:,:,n,m,2)=BBAux{2}(m,:,2)'*BBAux{1}(n,:,1);
        ts.BB2g(:,:,n,m,1)=BBAux{2}(m,:,1)'*BBAux{1}(n,:,3); 
        ts.BB2g(:,:,n,m,2)=BBAux{2}(m,:,3)'*BBAux{1}(n,:,1);
        ts.BB11g(:,:,n,m)=BBAux{2}(m,:,2)'*BBAux{1}(n,:,2); 
    end
end

% 2*Dp as stated in (8)
for l=1:L
    ts.BB1g(:,:,:,:,l)=2*ts.Dp(l)*ts.BB1g(:,:,:,:,l);
    ts.BB2g(:,:,:,:,l)=2*(ts.Dp(l)^2)*ts.BB2g(:,:,:,:,l);
end
ts.BB11g=2*ts.Dp(1)*ts.Dp(2)*ts.BB11g;
% bspline  products matriz  
ts.BB1g=permute(repmat(ts.BB1g,[1 1 1 1 1 N L]),[1 2 6 3 4 7 5]);
ts.BB2g=permute(repmat(ts.BB2g,[1 1 1 1 1 N L]),[1 2 6 3 4 7 5]);
ts.BB11g=permute(repmat(ts.BB11g,[1 1 1 1 N L]),[1 2 5 3 4 6]);

%Modification of coefficients to avoid borders
for l=1:ts.nt
   ts.coefg{l}(ts.coefg{l}<1)=1;
   ts.coefg{l}(ts.coefg{l}>tamxy(l))=tamxy(l);
end
clear BBAux