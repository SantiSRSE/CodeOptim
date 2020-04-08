%Para detalles de la notacion consultar la referencia:
%L. Cordero-Grande, S. Merino-Caviedes, S. Aja-Fernandez, and C. Alberola-Lopez,
%"Groupwise elastic registration by a new sparsity-promoting metric: application 
%to the alignment of cardiac magnetic resonance perfusion images,"
%IEEE Trans. Pattern Anal. Mach. Intell., vol. in press, 2013.
%clear all

addpath ./registrado;
addpath ./IO;
addpath ./preprocesado;

%ruta=('/nfs/export/pfc/ssanest/');% carpeta de imagenes, que es la de trabajo
%En este caso se utiliza la actual, pero se puede cambiar con rutas absolutas
ruta= [pwd '\']; 
%leer_imagenes
%[N I tam]=leer(ruta);%N, numero de imagenes
%[I_ref,meta_ref]=nrrdread([ruta 'SA-LE-label.nrrd']);
% I=dicomread('');
[I,meta]=nrrdread([ruta 'SA-C09.nrrd']);
%seleccion de corte inicial
init_slice(I)
uiwait(gcf)

load corte
inicio=corte(3) % para SA
delete corte.mat
tam=size(I);
%Mask=zeros([tam(1:2) tam(3) 3]);%LVen LVep LVmi
Mask=zeros(tam);

%sirve para hacer zoom
radio=(tam(1)-1)/2.1;%radio para la máscara
[~, ~, caja]=mascara(radio,tam);
%caja=[1 tam(1);1 tam(2)];

im=squeeze(I(caja(1,1):caja(1,2),caja(2,1):caja(2,2),:));
primer_corte=1;
secciones_completas=0;

while(1)
   if primer_corte, i=inicio; else i=numero; end
   figure(i)
   imshow(im(:,:,i),[]);%mostrar imagen del area seleccionada por caja
   xlabel(['LV endocardio corte nº' num2str(i) ])
   set(gcf,'Position',get(0,'screensize'))%maximizar
   while(1)
        h=imfreehand(gca);%LVen
        Mask(caja(1,1):caja(1,2),caja(2,1):caja(2,2),i)=h.createMask();%recoloco la máscara
        if(waitforbuttonpress)
            break;
        end  
        clear h;
        imshow(im(:,:,i),[]);
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);        
   end
%    xlabel(['LV epicardio corte nº' num2str(i) ])
%    while(1)
%         h=imfreehand(gca);%LVep
%         Mask(caja(1,1):caja(1,2),caja(2,1):caja(2,2),i,2)=h.createMask();
%         if(waitforbuttonpress)
%             break;
%         end
%         clear h;
%         imshow(im(:,:,i),[]);
%         set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%    end
%   Mask(:,:,i,3)=Mask(:,:,i,2)-Mask(:,:,i,1); Segm=squeeze(Mask(:,:,i,:));
%   save(['Corte#' num2str(i)],'Segm') 
   
%    while(1)
%         h=imfreehand(gca);%RVen
%         V(caja(1,1):caja(1,2),caja(2,1):caja(2,2),i,3)=h.createMask();
%         if(waitforbuttonpress)
%             break;
%         end
%         clear h;
%         imshow(im(:,:,i),[]);
%         set(gcf,'Position',get(0,'screensize'))
%    end
   close;
   
   if primer_corte
       choice = questdlg('Direccion?', 'Seleccion de corte','apice', 'base','base');
       switch choice
           case 'apice'
               numero = inicio+1; direccion=1;
           case 'base'
               numero = inicio-1; direccion=-1;
       end
       primer_corte=0;
   else
       choice = questdlg('Finalizada seccion?', 'Seleccion de corte','Si', 'No','No');
       switch choice
           case 'Si'
               numero = inicio-direccion; direccion=-direccion;
               secciones_completas=secciones_completas+1
           case 'No'
               numero = numero+direccion; 
           
       end
   end
   if secciones_completas==2, break, end
end


%guardo mascaras en un solo fichero
save('MascarasMiocardio','Mask')

%visualizo los resultados
for i=1:tam(3)
    if exist(['Corte#' num2str(i) '.mat'],'file')>0
        load(['Corte#' num2str(i)])
        Imagen = I(:,:,i);
        [Lines,Vertices]=isocontour(Segm(:,:,3),0.5);
        figure(i), imshow(Imagen,[]), hold on;
        V1=Vertices(Lines(:,1),:); V2=Vertices(Lines(:,2),:);
        plot([V1(:,2) V2(:,2)]',[V1(:,1) V2(:,1)]','b');
    end
end

