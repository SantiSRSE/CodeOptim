function [xn,TAux] = transformBSpline(x,tS,T)
% TRANSFORMAR genera una nueva malla de puntos a partir de una inicial
% teniendo en cuenta los desplazamientos definidos anteriormente
%
%Entradas
%x: malla de pixeles inicial 
%ts: estructura de almacenamiento de las variables bspline y del tipo de
%transformaci�n a implementar
%T: matriz de transformaci�n con los desplazamientos
%caja: determina la zona de las im�genes donde se operar�
%flagopt: flag de activaci�n del m�todo de baja memoria (activo en cero)
%
%Salidas
%xn: nueva malla generada a partir de la transformaci�n
%TAux: matriz con los valores de la transformaci�n recolocada para una
%comunicaci�n posterior con otras funciones
TAux=0;
switch tS.type
    case 't'%Translations
        xn=bsxfun(@plus,x,T);%desplazamiento de pixel       
    case 'r'%Rotations
        error('Not implemented')
    case 'tyr'%Translations and rotations   
        error('Not implemented')
    case 'b'
        %coeficientes separados
        %tic
        if tS.flagRlast==0
        r=cell(1,tS.L);
        NR=zeros(1,tS.L);
        tSE3=(tS.E+1)*ones(1,tS.L);
        for l=1:tS.L;r{l}=tS.ROI(l,1):tS.ROI(l,2);NR(l)=length(r{l});end
        
        %TAux=zeros([NR tS.L tS.N tSE3],'like',T);%TAux:variable auxiliar de tama�o ts.BB con los valores de T
        TAux=zeros([NR tS.L size(T,5) tSE3],'like',T);%para el caso en el que no esten todos los shots
                
        %transformation is assigned by dimension to avoid asign big matrixes
        for k=0:tS.E; coef{1}=k+tS.coef{1}(r{1},1); 
            T1=T(coef{1},:,:,:,:);
            for l=0:tS.E;coef{2}=l+tS.coef{2}(r{2},1);
                T2=T1(:,coef{2},:,:,:);
                for m=0:tS.E;coef{3}=m+tS.coef{3}(r{3},1);
                    T3=T2(:,:,coef{3},:,:);
                    TAux(:,:,:,:,:,k+1,l+1,m+1)=T3;
                end
            end
        end       
        %update fields only over the ROI
        xn=extractROI(multDimSum(bsxfun(@times,TAux,extractROI(tS.BB,tS.ROI,1)),6:8),tS.ROI,0);          
        xn=bsxfun(@plus,x,real(xn)); 
        %toc  
        
        %xn2=xn;
        %tic
        elseif tS.flagRlast==-1
            
        r=cell(1,tS.L);
        NR=zeros(1,tS.L);
        tSE3=(tS.E+1)*ones(1,tS.L);
        for l=1:tS.L;r{l}=tS.ROI(l,1):tS.ROI(l,2);NR(l)=length(r{l});end
        
        %TAux=zeros([NR tS.L tS.N tSE3],'like',T);%TAux:variable auxiliar de tama�o ts.BB con los valores de T
        TAux=zeros([NR 1 1 tSE3],'like',T);%para el caso en el que no esten todos los shots
        xn=zeros([NR tS.L size(T,5)],'like',T);        
        
        for p=1:tS.L
        for q=1:size(T,5)
        %transformation is assigned by dimension to avoid asign big matrixes
        for k=0:tS.E; coef{1}=k+tS.coef{1}(r{1},1); 
            T1=T(coef{1},:,:,p,q);
            for l=0:tS.E;coef{2}=l+tS.coef{2}(r{2},1);
                T2=T1(:,coef{2},:);
                for m=0:tS.E;coef{3}=m+tS.coef{3}(r{3},1);
                    T3=T2(:,:,coef{3});
                    TAux(:,:,:,1,1,k+1,l+1,m+1)=T3;
                end
            end
        end       
        %update fields only over the ROI
        xn(:,:,:,p,q)=extractROI(multDimSum(bsxfun(@times,TAux,extractROI(tS.BB,tS.ROI,1)),6:8),tS.ROI,0);          
        end
        end
        xn=bsxfun(@plus,x,real(xn)); 

        elseif tS.flagRlast==1
        %procedimiento alternativo de baja memoria para la iteracion final
        BB=permute(extractROI(tS.BB,tS.ROI,1),[6 7 8 4 5 1 2 3]); %recoloco ts.BB
        xn=zeros([size(x) tS.N],'like',T);
        for m=tS.ROI(1,1):tS.ROI(1,2) %margen de Dp
            for n=tS.ROI(2,1):tS.ROI(2,2)
                for o=tS.ROI(3,1):tS.ROI(3,2)
         xn(m,n,o,:,:)=multDimSum(bsxfun(@times,T(tS.coef{1}(m,1):tS.coef{1}(m,1)+tS.E, ...
         tS.coef{2}(n,1):tS.coef{2}(n,1)+tS.E,tS.coef{3}(o,1):tS.coef{3}(o,1)+tS.E,:,:),...
         BB(:,:,:,1,1,m,n,o)),1:3);

         %xn(m,n,o,ti,nd)=sum(sum(sum(BB(:,:,:,1,1,m,n,o).*TAux(ts.coef{1}(m,1):ts.coef{1}(m,1)+ts.E,ts.coef{2}(n,1):ts.coef{2}(n,1)+ts.E,ts.coef{3}(o,1):ts.coef{3}(o,1)+ts.E,ti,nd),3),2),1);
                end
            end
        end
        xn=bsxfun(@plus,x,extractROI(real(xn),tS.ROI,0));
        
        end
        %toc
        %sum(abs(xn(:)-xn2(:)))/numel(xn), pause
          
        
    otherwise
        error('Invalid argument: t for translation, r for rotation, tyr for both and b for bspline');
end
% IMPLEMENTAR CASO DE BAJA MEMORIA, ESTA EN TRANSFORMARND
