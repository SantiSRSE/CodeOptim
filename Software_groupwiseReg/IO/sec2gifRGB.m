
function sec2gifRGB( im, outfile)
%sec2gifRGB Saves stack of images as a animated gif file  
%   sec2gif( im, outfile, DelayTime) saves 3D RGB image im in outfile as a gif
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
rgb = complex2rgb(im);
DelayTime  = 0.05;

[ind, map] = rgb2ind(rgb(:,:,:,1), 256);
imwrite( ind,map,outfile, ...
    'gif','LoopCount',Inf,'DelayTime', DelayTime);

for i=2:size(im,4);
    [ind, map] = rgb2ind(rgb(:,:,:,i), 256);
    imwrite(ind,map,outfile, ...
        'gif','WriteMode','append','DelayTime',DelayTime);
end


end
