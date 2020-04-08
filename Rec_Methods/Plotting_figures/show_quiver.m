function  show_quiver( x,xn,x0,Dimen)
%quiver over images
%revise permute
xn=permute(xn,[6 7 8 5 10 1 2 3 4 9]); x0=permute(x0,[6 7 8 5 10 1 2 3 4 9]);
x=abs(x);
if Dimen==2
    for i=1:size(xn,4)
        if size(x,4)>1
            figure, imshow(x(:,:,1,i-1),[])
        else
            %figure, 
            imshow(x(:,:,1),[])
        end
        hold on
        int1=1:10:size(xn,1); int2=1:10:size(xn,2);
        quiver(x0(int1,int2,1,1,2),x0(int1,int2,1,1,1),xn(int1,int2,1,i,2)-x0(int1,int2,1,1,2),xn(int1,int2,1,i,1)-x0(int1,int2,1,1,1),0), pause(0.01)
        hold off
        return
    end
end

end