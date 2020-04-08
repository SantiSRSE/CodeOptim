function tS=buildBSpline(x,Field,N,permT)
% CREABSPLINE almacena los elementos necesarios para la transformaci�n
% bspline, tanto las matrices de productos B-spline (hasta la derivada
% primera) como los coeficientes asociados, as� como la ubicaci�n de los puntos
% de control como su densidad, el orden del bspline, de manera que
% faciliten un posterior procesado
%
%Entradas
%x: malla inicial de pixeles
%caja: determina la zona donde se operar� posteriormente
%Dp: Separaci�n entre puntos de control
%E: orden del bspline
%N: n�mero de im�genes
%
%Salidas
%tS: estructura donde se almacenan las variables relacionadas con el
%bspline

gpu=isa(x{1},'gpuArray');

L=length(x);%numero de dimensiones
NX=zeros(1,L);
for l=1:L;NX(l)=length(x{l});end

%Definici�n de los parametros requeridos (coord. imagen)
NT=[NX(1:3) NX(1:3)];%tama�o

tS=struct('type','b','dist','unif','L',L,'Dp',Field.Dp, 'E',Field.order,'c',(1+NT(1:3))/2,'N',N,'ROI',Field.ROI,'rp',(Field.order+1)*Field.Dp/2);%rp: radio de influencia
for l=1:L%dimensiones de la imagen
    for b=1:2%min-max
        tS.C(l,b)=single(floor((abs(tS.c(l)-tS.ROI(l,b))+tS.rp(l))/tS.Dp(l)));
    end
end    
% ordenaci�n de cada punto de control
[tS.u(:,:,:,1),tS.u(:,:,:,2),tS.u(:,:,:,3)]=ndgrid(-tS.C(1,1):tS.C(1,2),-tS.C(2,1):tS.C(2,2),-tS.C(3,1):tS.C(3,2));

permD=[1 3 4 2];
tS.pu=bsxfun(@plus,permute(tS.c,permD),bsxfun(@times,permute(tS.Dp,permD),tS.u));%ubicaci�n en la malla de cada punto de control


%coeficientes y matriz de productos para la transformacion
coef=cell(1,L);
for l=1:L;coef{l}=zeros([NX(l) 2],'like',x{l});end%coeficientes originales sin offset
for b=1:2%min-max
   if rem(tS.E,2)~=0%orden 2 
       for l=1:L;coef{l}(:,b)=floor((x{l}(:)-tS.c(l))/tS.Dp(l))+((tS.E+(-1)^b)/2)*(-1)^b;end
   elseif rem(tS.E,2)==0%orden 1 y 3
       for l=1:L;coef{l}(:,b)=round((x{l}(:)-tS.c(l))/tS.Dp(l))+(tS.E/2)*(-1)^b;end
   end
end

%inicializacion de matriz auxiliar para el c�lculo (2->derivada) 
tS.BBs=cell(1,L);
for l=1:L;tS.BBs{l}=zeros([NX(l) tS.E+1 3],'like',coef{l});end    
   
%matriz auxiliar para los coeficientes de cada dimension
for l=1:L
    xaux=(x{l}-tS.c(l))/tS.Dp(l);
    for n=1:NX(l) %margen de toda la imagen para la transformacion
        vk=coef{l}(n,1):coef{l}(n,2);
        vki=vk-coef{l}(n,1)+1;
        for d=1:3;tS.BBs{l}(n,vki,d)=BSplineFun(xaux(n)-vk,tS.E,d-1);end
    end
end

%ts.BB mesh
tSE3=(tS.E+1)*ones(1,3);
tS.BBs{1}=permute(tS.BBs{1},[2 4 5 1 6 7 3]);
tS.BBs{2}=permute(tS.BBs{2},[4 2 5 6 1 7 3]);
tS.BBs{3}=permute(tS.BBs{3},[4 5 2 6 7 1 3]);
tS.BB=bsxfun(@times,bsxfun(@times,dynInd(tS.BBs{1},1,7),dynInd(tS.BBs{2},1,7)),dynInd(tS.BBs{3},1,7));

pv=[1 2 3];
if Field.lambda(1)>0
    tS.BB1=zeros([tSE3 NX L],'like',tS.BB);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Primera derivada
    for l=1:3
        tS.BB1=dynInd(tS.BB1,l,7,tS.Dp(pv(1))*bsxfun(@times,bsxfun(@times,dynInd(ts.BBs{pv(1)},2,7),dynInd(ts.BBs{pv(2)},1,7)),dynInd(ts.BBs{pv(3)},1,7)));
        pv=circshift(pv,[0 -1]);
    end
end
if Field.lambda(2)>0        
    tS.BB11=zeros([tSE3 NX L],'like',tS.BB);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Derivada cruzada
    tS.BB2=zeros([tSE3 NX L],'like',tS.BB);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Segunda derivada
    for l=1:3
        tS.BB11=dynInd(tS.BB11,l,7,(tS.Dp(pv(2))*tS.Dp(pv(3)))*bsxfun(@times,bsxfun(@times,dynInd(ts.BBs{pv(1)},1,7),dynInd(ts.BBs{pv(2)},2,7)),dynInd(ts.BBs{pv(3)},2,7)));
        tS.BB2=dynInd(tS.BB2,l,7,(tS.Dp(pv(1))^2)*bsxfun(@times,bsxfun(@times,dynInd(ts.BBs{pv(1)},3,7),dynInd(ts.BBs{pv(2)},1,7)),dynInd(ts.BBs{pv(3)},1,7)));
        pv=circshift(pv,[0 -1]);
    end
end
for l=1:L;tS.BBs{l}=dynInd(tS.BBs{l},1,7);end%So far only interested in the 0-order term

%Se a�ade un offset a los coeficientes de caja, se utilizan en transformacion 
for l=1:L;tS.coef{l}=coef{l}+1-coef{l}(Field.ROI(l,1),1);end% coeficientes con offset


if Field.Full==0
    if ~isempty(permT)
    [tS.BB,tS.BBs{1},tS.BBs{2},tS.BBs{3}]=parUnaFun({tS.BB,tS.BBs{1},tS.BBs{2},tS.BBs{3}},@permute,permT);
    if Field.lambda(1)>0;[tS.BB1]=parUnaFun({tS.BB1},@permute,permT);end
    if Field.lambda(2)>0;[tS.BB11,tS.BB2]=parUnaFun({tS.BB11,tS.BB2},@permute,permT);end
    end
    tS=rmfield(tS,{'pu','u'});
    return,
end


%coeficientes y matriz de productos para el c�lculo del gradiente
NP=size(tS.pu);

%metodo de coeficientes separados
for l=1:L;tS.coefg{l}=zeros([NP(l) 2],'like',tS.pu);end

xbis{1}=tS.pu(:,1,1,1);   
xbis{2}=tS.pu(1,:,1,2); 
xbis{3}=tS.pu(1,1,:,3);


for b=0:1
    for l=1:L;tS.coefg{l}(:,b+1)=ceil(xbis{l}(:)-ceil(tS.rp(l)))+b*2*ceil(tS.rp(l));end
end%min-max

ts.BBgs=cell(1,L);
for l=1:L;tS.BBgs{l}=zeros([NP(l) 2*ceil(tS.rp(l))+1 3],'like',xbis{l});end

%matriz auxiliar para los coeficientes de cada dimension (dada por xbis2)
for l=1:L    
    for n=1:length(xbis{l})
        vk=tS.coefg{l}(n,1):tS.coefg{l}(n,2);
        vki=vk-tS.coefg{l}(n,1)+1;
        for d=1:3;tS.BBgs{l}(n,vki,d)=BSplineFun(((vk-xbis{l}(n))/tS.Dp(l)),tS.E,d-1);end
    end
    if gpu;tS.BBgs{l}=gpuArray(tS.BBgs{l});end
end

% matriz de productos bspline 
NG=2*ceil(tS.rp(1:3))+1;
tS.BBgs{1}=permute(tS.BBgs{1},[1 6 7 2 4 5 3]);
tS.BBgs{2}=permute(tS.BBgs{2},[6 1 7 4 2 5 3]);
tS.BBgs{3}=permute(tS.BBgs{3},[6 7 1 4 5 2 3]);

tS.BBg=bsxfun(@times,bsxfun(@times,dynInd(tS.BBgs{1},1,7),dynInd(tS.BBgs{2},1,7)),dynInd(tS.BBgs{3},1,7));
if Field.lambda(1)>0
    tS.BB1g=zeros([NP(1:3) NG L],'like',tS.BBg);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Primera derivada
    for l=1:3
        tS.BB1g=dynInd(tS.BB1g,l,7,tS.Dp(pv(1))*bsxfun(@times,bsxfun(@times,dynInd(tS.BBgs{pv(1)},2,7),dynInd(tS.BBgs{pv(2)},1,7)),dynInd(tS.BBgs{pv(3)},1,7)));
        pv=circshift(pv,[0 -1]);
    end
end
if Field.lambda(2)>0        
    tS.BB11g=zeros([NP(1:3) NG L],'like',tS.BBg);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Derivada cruzada
    tS.BB2g=zeros([NP(1:3) NG L],'like',tS.BBg);%matriz de productos bspline Bprimax*By y Bprimay*Bx-Segunda derivada
    for l=1:3
        tS.BB11g=dynInd(tS.BB11g,l,7,(tS.Dp(pv(2))*tS.Dp(pv(3)))*bsxfun(@times,bsxfun(@times,dynInd(tS.BBgs{pv(1)},1,7),dynInd(tS.BBgs{pv(2)},2,7)),dynInd(tS.BBgs{pv(3)},2,7)));
        tS.BB2g=dynInd(tS.BB2g,l,7,(tS.Dp(pv(1))^2)*bsxfun(@times,bsxfun(@times,dynInd(tS.BBgs{pv(1)},3,7),dynInd(tS.BBgs{pv(2)},1,7)),dynInd(tS.BBgs{pv(3)},1,7)));
        pv=circshift(pv,[0 -1]);
    end
end
for l=1:L;tS.BBgs{l}=dynInd(tS.BBgs{l},1,7);end%So far only interested in the 0-order term

%Modificaci�n de los coeficientes para evitar desbordes
for l=1:tS.L
    tS.coefgr{l}=min(max(tS.coefg{l},1),NX(l));
    tS.dim{l}=1+tS.coefgr{l}(:,1)-tS.coefg{l}(:,1);
    tS.dim{l}(:,2)=tS.coefgr{l}(:,2)-tS.coefgr{l}(:,1)+tS.dim{l};
end
tS=rmfield(tS,{'pu','u'});

%take care of dimensions
if ~isempty(permT)
    [tS.BB,tS.BBg,tS.BBs{1},tS.BBs{2},tS.BBs{3},tS.BBgs{1},tS.BBgs{2},tS.BBgs{3}]=parUnaFun({tS.BB,tS.BBg,tS.BBs{1},tS.BBs{2},tS.BBs{3},tS.BBgs{1},tS.BBgs{2},tS.BBgs{3}},@permute,permT);
    if Field.lambda(1)>0;[tS.BB1,tS.BB1g]=parUnaFun({tS.BB1,tS.BB1g},@permute,permT);end
    if Field.lambda(2)>0;[tS.BB11,tS.BB2,tS.BB11g,tS.BB2g]=parUnaFun({tS.BB11,tS.BB2,tS.BB11g,tS.BB2g},@permute,permT);end
end