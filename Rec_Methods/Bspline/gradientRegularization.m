function gR=gradientRegularization(rec,tR)
% GRADIENTE_REGULARIZAR actualiza los t�rminos de suavidad espacial 
% envueltos en el c�lculo del gradiente 
%
%Entradas
%ts: estructura de almacenamiento de las variables bspline
%term: estructura de almacenamiento de las variables de suavidad temporal y
%espacial
%
%Salidas
%term: estructura de almacenamiento de las variables de suavidad temporal y
%espacial, con los nuevos valores necesarios en el gradiente

NG=size(tS.BBg);

if any(rec.Field.lambda~=0)
    %Derivada espacial dada en (2)
    if rec.Field.lambda(1)>0;gR.dr1=zeros([NG(1:3) 1 rec.tS.N NG(6:8)],'like',tS.BBg);end
    if rec.Field.lambda(2)>0;gR.dr2=zeros([NG(1:3) 1 rec.tS.N NG(6:8)],'like',tS.BBg);end

    %support of the bspline
    coefgr=cell(1,3);dim=cell(1,3);
    for k=1:NG(6);coefgr{1}=rec.tS.coefgr{1}(k,1):rec.tS.coefgr{1}(k,2);dim{1}=rec.tS.dim{1}(k,1):rec.tS.dim{1}(k,2);
        for l=1:NG(7);coefgr{2}=rec.tS.coefgr{2}(l,1):rec.tS.coefgr{2}(l,2);dim{2}=rec.tS.dim{2}(l,1):rec.tS.dim{2}(l,2);
            for n=1:NG(8);coefgr{3}=rec.tS.coefgr{3}(n,1):rec.tS.coefgr{3}(n,2);dim{3}=rec.tS.dim{3}(n,1):rec.tS.dim{3}(n,2);       
                %transformation * precomputed derivative
                if rec.Field.Lambda(1)>0
    gR.dr1(dim{1},dim{2},dim{3},1,:,k,l,n)=gR.dr1(dim{1},dim{2},dim{3},1,:,k,l,n)+multDimSum(bsxfun(@times,tR.dr1(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:).*ts.BB1g(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:)),[4 9]);
                end
                if rec.Field.Lambda(2)>0
    gR.dr2(dim{1},dim{2},dim{3},1,:,k,l,n)=gR.dr2(dim{1},dim{2},dim{3},1,:,k,l,n)+multDimSum(bsxfun(@times,tR.dr2(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:).*ts.BB2g(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:)),[4 9]);       
    gR.dr2(dim{1},dim{2},dim{3},1,:,k,l,n)=gR.dr2(dim{1},dim{2},dim{3},1,:,k,l,n)+multDimSum(bsxfun(@times,tR.dr11(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:).*ts.BB11g(coefgr{1},coefgr{2},coefgr{3},:,:,k,l,n,:)),[4 9]);
                end
            end
        end
    end
else
    gR=[];
end
