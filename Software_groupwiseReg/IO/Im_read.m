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

function [I, Im_size, Mask, Tdias]=Im_read(route)
% Im_read read DICOM images from the folder indicated
%
%Input
%route: Image folder  
%
%Output
%I: Image set (double) 
%Im_size: images size
%Mask: miocardical mask in diastole 
%Tdias: correspondence to diastolic phase
dicom=0; % read from dicom or not
if dicom
    A=dir(strcat(route,'IM*'));% IM image identifier
    N=length(A);%number of images
    
    for k=1:N
        im=strcat(route,['IM' num2str(k)]);
        filename=sprintf('%s',im); %build direction
        Image= dicomread(filename);%load
        I(:,:,k)=double(Image);%3D array without displacement
        
    end
    Im_size=size(I);%tamaño final de las imágenes 
    load([route 'Mask']); Mask=M; %diastolic mask
else
    %subsampling factor x2
    load([route 'Images']);
    load([route 'Mask']); Mask=imresize(M,0.5); %diastolic mask
end
    Im_size=size(I);%tamaño final de las imágenes 

Tdias=1;