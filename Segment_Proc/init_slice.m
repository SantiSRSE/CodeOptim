function varargout = init_slice(varargin)
% VISVOL3DRT M-file for visvol3drt.fig
%      VISVOL3DRT, by itself, creates a new VISVOL3DRT or raises the existing
%      singleton*.
%
%      H = VISVOL3DRT returns the handle to a new VISVOL3DRT or the handle to
%      the existing singleton*.
%
%      VISVOL3DRT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISVOL3DRT.M with the given input arguments.
%
%      VISVOL3DRT('Property','Value',...) creates a new VISVOL3DRT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visvol3drt_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visvol3drt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visvol3drt

% Last Modified by GUIDE v2.5 04-May-2016 13:02:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       'init_slice', ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @init_slice_OpeningFcn, ...
                   'gui_OutputFcn',  @init_slice_OutputFcn, ... % 
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before demoMascaraSA is made visible.
function init_slice_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demoMascaraSA (see VARARGIN)

% Choose default command line output for demoMascaraSA
handles.output = hObject;

%% Enlazamos los datos a la interfaz
SA_C18=varargin{1};
handles.vol3drt = single(SA_C18-min(SA_C18(:)))./single(max(SA_C18(:))--min(SA_C18(:)));
handles.pos = [1 1 1 1];
handles.oldpos = [0 0 0 0];
handles.range = double([min(handles.vol3drt(:)) max(handles.vol3drt(:))]);
    
%Modificamos los sliders para que se ajusten al tamaño del volumen
%... y que sus pasos sean enteros. POS y OLDPOS siguen la convencion XYZT
SIZE = size(SA_C18);
handles.pos([2 1 3]) = max(round(0.5*SIZE(1:3)),1);
set(handles.sliderX,'Min',1,'Max',SIZE(2),'Value',handles.pos(1),'SliderStep',[1/(SIZE(2)-1) 5/(SIZE(2)-1)]);
set(handles.sliderY,'Min',1,'Max',SIZE(1),'Value',handles.pos(2),'SliderStep',[1/(SIZE(1)-1) 5/(SIZE(1)-1)]);
set(handles.sliderZ,'Min',1,'Max',SIZE(3),'Value',handles.pos(3),'SliderStep',[1/(SIZE(3)-1) 5/(SIZE(3)-1)]);
if length(SIZE)>=4
    set(handles.sliderT,'Min',1,'Max',SIZE(4),'Value',handles.pos(4),'SliderStep',[1/(SIZE(4)-1) 5/(SIZE(4)-1)]);
else
    set(handles.sliderT,'Min',0,'Max',1,'Value',1,'Enable','off');
end

%Mostramos el valor del slider en los cuadros de texto
set(handles.txtEjeX,'String',num2str(handles.pos(1)));
set(handles.txtEjeY,'String',num2str(handles.pos(2)));
set(handles.txtEjeZ,'String',num2str(handles.pos(3)));
set(handles.txtEjeT,'String',num2str(handles.pos(4)));

%Mostramos los cortes iniciales
showVisualization(hObject,handles)

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = init_slice_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.pos;


% --- Executes on slider movement.
function sliderX_Callback(hObject, eventdata, handles)
% hObject    handle to sliderX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pos = handles.pos;
handles.oldpos = pos;
pos(1) = round(get(hObject,'Value'));
handles.pos = pos;
guidata(hObject, handles);

set(handles.txtEjeX,'String',num2str(handles.pos(1)));

showVisualization(hObject,handles);





% --- Executes during object creation, after setting all properties.
function sliderX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderY_Callback(hObject, eventdata, handles)
% hObject    handle to sliderY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pos = handles.pos;
handles.oldpos = pos;
pos(2) = round(get(hObject,'Value'));
handles.pos = pos;
guidata(hObject, handles);

set(handles.txtEjeY,'String',num2str(handles.pos(2)));

showVisualization(hObject,handles);



% --- Executes during object creation, after setting all properties.
function sliderY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pos = handles.pos;
handles.oldpos = pos;
pos(3) = round(get(hObject,'Value'));
handles.pos = pos;
guidata(hObject, handles);

set(handles.txtEjeZ,'String',num2str(handles.pos(3)));

showVisualization(hObject,handles);



% --- Executes during object creation, after setting all properties.
function sliderZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderT_Callback(hObject, eventdata, handles)
% hObject    handle to sliderT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pos = handles.pos;
handles.oldpos = pos;
pos(4) = round(get(hObject,'Value'));
handles.pos = pos;
guidata(hObject, handles);

set(handles.txtEjeT,'String',num2str(handles.pos(4)));

showVisualization(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnPointsYZ.
function btnPointsYZ_Callback(hObject, eventdata, handles)
% hObject    handle to btnPointsYZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axesX)
%Tomamos la posicion de los puntos de la imagen
[x2,x3] = ginput(1);
x2 = round(x2); x3 = round(x3);
x1 = round(get(handles.sliderX,'Value'));
x4 = round(get(handles.sliderT,'Value')); 
%Sacamos los puntos por pantalla
disp([x1 x2 x3 x4])


% --- Executes on button press in btnPointsXZ.
function btnPointsXZ_Callback(hObject, eventdata, handles)
% hObject    handle to btnPointsXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesY)
%Tomamos la posicion de los puntos de la imagen
[x1,x3] = ginput(1);
x1 = round(x1); x3 = round(x3);
x2 = round(get(handles.sliderY,'Value'));
x4 = round(get(handles.sliderT,'Value')); 
%Sacamos los puntos por pantalla
disp([x1 x2 x3 x4])


% --- Executes on button press in btnPointsXY.
function btnPointsXY_Callback(hObject, eventdata, handles)
% hObject    handle to btnPointsXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesZ)
%Tomamos la posicion de los puntos de la imagen
[x1,x2] = ginput(1);
x1 = round(x1); x2 = round(x2); 
x3 = round(get(handles.sliderZ,'Value'));
x4 = round(get(handles.sliderT,'Value')); 
%Sacamos los puntos por pantalla
disp([x1 x2 x3 x4])


function showVisualization(hObject,handles)
%Se encarga de llamar a unas funciones u otras segun que visualizacion
%queramos mostrar
showslices(hObject,handles)


function showslices(hObject,handles)
% Esta funcion se encarga de mostrar por pantalla los cortes que han
% cambiado respecto a la nueva posicion de los slices.

%Mostramos los cortes iniciales
axes(handles.axesX); %PLANO YZ
I =  im2double(squeeze(handles.vol3drt(:,handles.pos(1),:,handles.pos(4)))');
imshow( cat(3,I,I,I),handles.range);
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(3)*[1 1],'ZData',[1 1],'Color','b')
line('XData',handles.pos(2)*[1 1],'YData',[1 size(handles.vol3drt,3)],'ZData',[1 1],'Color','g')
hold off

axes(handles.axesY); %PLANO XZ
I = im2double(squeeze(handles.vol3drt(handles.pos(2),:,:,handles.pos(4)))');
imshow( cat(3,I,I,I),handles.range);
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(3)*[1 1],'ZData',[1 1],'Color','b')
line('XData',handles.pos(1)*[1 1],'YData',[1 size(handles.vol3drt,3)],'ZData',[1 1],'Color','r')
hold off

axes(handles.axesZ); %PLANO XY
I =  im2double(squeeze(handles.vol3drt(:,:,handles.pos(3),handles.pos(4))));
imshow( cat(3,I,I,I),handles.range);
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(2)*[1 1],'ZData',[1 1],'Color','g')
line('XData',handles.pos(1)*[1 1],'YData',[1 size(handles.vol3drt,1)],'ZData',[1 1],'Color','r')
hold off

handles.oldpos = handles.pos;
guidata(hObject,handles);


function showlabels(hObject,handles)
% Esta funcion se encarga de mostrar por pantalla los cortes que han
% cambiado respecto a la nueva posicion de los slices.

%Mostramos los cortes iniciales
axes(handles.axesX); %PLANO YZ
I =  im2double(squeeze(handles.vol3drt(:,handles.pos(1),:,handles.pos(4)))');
L =  squeeze(handles.labels3drt(:,handles.pos(1),:,handles.pos(4)))';
HSV = cat(3,double(handles.hues(L+1)),0.5*ones(size(L)),I);
imshow(hsv2rgb(HSV))
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(3)*[1 1],'ZData',[1 1],'Color','b')
line('XData',handles.pos(2)*[1 1],'YData',[1 size(handles.vol3drt,3)],'ZData',[1 1],'Color','g')
hold off

axes(handles.axesY); %PLANO XZ
I = im2double(squeeze(handles.vol3drt(handles.pos(2),:,:,handles.pos(4)))');
L =  squeeze(handles.labels3drt(handles.pos(2),:,:,handles.pos(4)))';
HSV = cat(3,double(handles.hues(L+1)),0.5*ones(size(L)),I);
imshow(hsv2rgb(HSV))
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(3)*[1 1],'ZData',[1 1],'Color','b')
line('XData',handles.pos(1)*[1 1],'YData',[1 size(handles.vol3drt,3)],'ZData',[1 1],'Color','r')
hold off

axes(handles.axesZ); %PLANO XY
I =  im2double(squeeze(handles.vol3drt(:,:,handles.pos(3),handles.pos(4))));
L =  squeeze(handles.labels3drt(:,:,handles.pos(3),handles.pos(4)));
HSV = cat(3,double(handles.hues(L+1)),0.5*ones(size(L)),I);
imshow(hsv2rgb(HSV))
hold on
line('XData',[1 size(handles.vol3drt,2)],'YData',handles.pos(2)*[1 1],'ZData',[1 1],'Color','g')
line('XData',handles.pos(1)*[1 1],'YData',[1 size(handles.vol3drt,1)],'ZData',[1 1],'Color','r')
hold off

handles.oldpos = handles.pos;
guidata(hObject,handles);


% --- Executes on button press in Button.
function Button_Callback(hObject, eventdata, handles)
% hObject    handle to Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout = init_slice_OutputFcn(hObject, eventdata, handles); 
corte=varargout(1:3); save corte corte
close all;
