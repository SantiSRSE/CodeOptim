function sec2gif( im, outfile, DelayTime)
%sec2gif Saves stack of images as a animated gif file  
%   sec2gif( im, outfile, DelayTime) saves 3D array im in outfile as a gif
%   animated image. DelayTime is the elapsed time between frames.
%
% 2014 Javier Royuela del Val <jroyval@lpi.tel.uva.es>
% Laboratorio de Procesado de Imagen, Universidad de Valladolid
% www.lpi.tel.uva.es
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

im = uint8(abs(im)*(255/max(abs(im(:)))));

if nargin < 3
    DelayTime  = 0.05;
end

map = gray(256);

im = squeeze(im);
imwrite(im(:,:,1),map,outfile, ...
    'gif','LoopCount',Inf,'DelayTime', DelayTime);

for i=2:size(im,3)
    imwrite(im(:,:,i),map,outfile, ...
        'gif','WriteMode','append','DelayTime',DelayTime);
end


end
