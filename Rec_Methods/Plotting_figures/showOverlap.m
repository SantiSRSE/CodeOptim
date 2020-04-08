function showOverlap(y1,y2,typ)
%show initial and final shots in different colors to notice overlap

N=size(y1);N(end+1:5)=1;
if ~strcmp(typ,'sequential')
    Nss=sqrt(N(5));
    y1=reshape(y1,[N(1) 1 N(2) 1 Nss Nss]);
    y1=permute(y1,[1 5 3 6 2 4]);
    IMS(:,:,1)=reshape(y1,Nss*N(1:2));
    if size(y2,5)==1
        IMS(:,:,2)=repmat(y2,[Nss Nss]);
    else
        y2=reshape(y2,[N(1) 1 N(2) 1 Nss Nss]);
        y2=permute(y2,[1 5 3 6 2 4]);
        IMS(:,:,2)=reshape(y2,Nss*N(1:2));
    end
else
    IMS(:,:,1)=reshape(y1,[N(1) N(2)*N(5)]);
    if size(y2,5)==1
        IMS(:,:,2)=repmat(y2,[1 N(5)]);
    else
        IMS(:,:,2)=reshape(y2,[N(1) N(2)*N(5)]);
    end
end
IMS(:,:,3)=0;
IMS=sqrt(abs(IMS));%SQRT TO IMPROVE CONTRAST A BIT
IMS=IMS/max(IMS(:));

figure
imshow(IMS)
set(gcf,'Position',get(0,'ScreenSize'))

