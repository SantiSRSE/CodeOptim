function [xn,xnBplus] = ResampleField( C, sm,Mult_field,sm2,Dimen,gpu)
%From XCAT structure to rec.T structure
% the applied displacement is modulated by Mult


    %sm=[x y time]
    d=Dimen;
   
    if gpu>0
        sm=gpuArray(sm); sm2=gpuArray(sm2);
    end
    
    % BUILD XN and interpollating fields
    [xn(:,:,:,1,1), xn(:,:,:,1,2), xn(:,:,:,1,3)] = ndgrid(1:sm(1), 1:sm(2), 1:sm(3));
    [xnLin(:,:,:,1,1), xnLin(:,:,:,1,2), xnLin(:,:,:,1,3)] = ndgrid(linspace(1,sm(1),sm2(1)), linspace(1,sm(2),sm2(2)),linspace(1,sm(3),sm2(3)));
    [xnBplus(:,:,:,1,1), xnBplus(:,:,:,1,2), xnBplus(:,:,:,1,3)] = ndgrid(1:sm2(1), 1:sm2(2), 1:sm2(3));
    xn = repmat(xn,[1 1 1 sm(5) 1]);
    
    
    for nn = 2:2%sm(5)
        if gpu>0
            Cn = gpuArray(C{nn});
        else
            Cn = C{nn};
        end
        np = size(Cn,1);
        
        %idxn1 = sub2ind(size(xn),Cn(:,1), Cn(:,2), nn*ones(np,1), 1*ones(np,1));
        %idxn2 = sub2ind(size(xn),Cn(:,1), Cn(:,2), nn*ones(np,1), 2*ones(np,1));
        %xn(idxn1) = Cn(:,5);
        %xn(idxn2) = Cn(:,4);
        
        idxn1 = sub2ind(size(xn),Cn(:,1)+1, Cn(:,2)+1, Cn(:,3)+1, nn*ones(np,1), 1*ones(np,1));
        idxn2 = sub2ind(size(xn),Cn(:,1)+1, Cn(:,2)+1, Cn(:,3)+1, nn*ones(np,1), 2*ones(np,1));
        idxn3 = sub2ind(size(xn),Cn(:,1)+1, Cn(:,2)+1, Cn(:,3)+1, nn*ones(np,1), 3*ones(np,1));
        %idxn1 = sub2ind(size(xn),Cn(:,1), Cn(:,2), max(Cn(:,3),1), nn*ones(np,1), 1*ones(np,1));
        %idxn2 = sub2ind(size(xn),Cn(:,1), Cn(:,2), max(Cn(:,3),1), nn*ones(np,1), 2*ones(np,1));
        %idxn3 = sub2ind(size(xn),Cn(:,1), Cn(:,2), max(Cn(:,3),1), nn*ones(np,1), 3*ones(np,1));
        xn(idxn1) = Cn(:,4);
        xn(idxn2) = Cn(:,5);
        xn(idxn3) = Cn(:,6);
        
        disp=(xn(:,:,:,nn,:)-xn(:,:,:,1,:))*Mult_field; %sm2(1)/sm(1)
        
        if gpu>0, disp=gpuArray(disp); end
                
        if d==2
            %gaussian smoothing kernel
            H = fspecial('gaussian',60,50);
            disp2(:,:,1,1,1) = imfilter(disp(:,:,1,1,2),H,'replicate'); %ojo dimensiones
            disp2(:,:,1,1,2) = imfilter(disp(:,:,1,1,1),H,'replicate');
            disp2(:,:,1,1,3) = disp(:,:,1,1,3);
            xnBplus(:,:,1,nn,1) = interpn(xn(:,:,1,1,1),xn(:,:,1,1,2),disp2(:,:,1,1,1),xnLin(:,:,1,1,1),xnLin(:,:,1,1,2),'linear')+xnBplus(:,:,1,1,1);
            xnBplus(:,:,1,nn,2) = interpn(xn(:,:,1,1,1),xn(:,:,1,1,2),disp2(:,:,1,1,2),xnLin(:,:,1,1,1),xnLin(:,:,1,1,2),'linear')+xnBplus(:,:,1,1,2);
            xnBplus(:,:,1,nn,3) = xnBplus(:,:,1,1,3);
            
        elseif d==3
            %gaussian smoothing kernel
            H = fspecial('gaussian',60,50);
            disp2(:,:,:,1,1) = imfilter(disp(:,:,:,1,2),H,'replicate');
            disp2(:,:,:,1,2) = imfilter(disp(:,:,:,1,1),H,'replicate');
            disp2(:,:,:,1,3) = imfilter(disp(:,:,:,1,3),H,'replicate');
            xnBplus(:,:,:,nn,1) = interpn(xn(:,:,:,1,1),xn(:,:,:,1,2),xn(:,:,:,1,3),disp2(:,:,:,1,1),xnLin(:,:,:,1,1),xnLin(:,:,:,1,2),xnLin(:,:,:,1,3),'linear')+xnBplus(:,:,:,1,1);
            xnBplus(:,:,:,nn,2) = interpn(xn(:,:,:,1,1),xn(:,:,:,1,2),xn(:,:,:,1,3),disp2(:,:,:,1,2),xnLin(:,:,:,1,1),xnLin(:,:,:,1,2),xnLin(:,:,:,1,3),'linear')+xnBplus(:,:,:,1,2);
            xnBplus(:,:,:,nn,3) = interpn(xn(:,:,:,1,1),xn(:,:,:,1,2),xn(:,:,:,1,3),disp2(:,:,:,1,3),xnLin(:,:,:,1,1),xnLin(:,:,:,1,2),xnLin(:,:,:,1,3),'linear')+xnBplus(:,:,:,1,3);
        end
    end
    clear disp2 disp xnLin
    %revise permute
    xn=permute(xn,[6 7 8 9 4 1 2 3 10 5]);
    xnBplus=permute(xnBplus,[6 7 8 9 4 1 2 3 10 5]);
    
    
end

