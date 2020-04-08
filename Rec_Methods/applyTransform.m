function x=applyTransform(x,T,d)
%ellastic transformation to a 2D or 3D image

NST=size(T,5);
NSX=size(x,5);
if NSX==1;x=repmat(x,[1 1 1 1 NST]);end %check if applying to an image or shots
for s=1:NST
    Ts=T(:,:,:,:,s); %extract mesh
    if d==2
       x(:,:,:,:,s)=interpn(x(:,:,:,:,s),Ts(:,:,:,1),Ts(:,:,:,2),'linear');
    else
       x(:,:,:,:,s)=interpn(x(:,:,:,:,s),Ts(:,:,:,1),Ts(:,:,:,2),Ts(:,:,:,3),'linear'); 
    end
    %linear interpolation, only one available in gpu
end
x(isnan(x))=0;
