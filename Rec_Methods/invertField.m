 function xn_inv = invertField(xn,lambda,thres,max_iter,gpu)
%   Inv motion field numericacly by min ||T^-1(T(x)) - x||^2
%   as in Metz et al 2011, and then obtain back
%
%   Javier Royuela del Val <jroyval@lpi.tel.uva.es>
 
    rule='combined';
    %revise permute
    xn=permute(single(xn),[1 2 3 5 4]);
        
    s = size(xn);

    ROW = s(1);    COL = s(2);     SLI = s(3);     N   = s(4);
    
    [x(:,:,:,1), x(:,:,:,2), x(:,:,:,3)] = ndgrid(1:ROW,1:COL,1:SLI);
    x=single(x);
    if gpu >0, x=gpuArray(x); xn=gpuArray(xn);end 
    
    n = ROW*COL*SLI*N;
    stop_criteria = false;

    xn_inv=permute(repmat(x,[1 1 1 1 N]),[1 2 3 5 4]);
    
    [~, g] = fun(xn_inv,x,xn,SLI);
    
    %conjugate gradient
    Dif = -g;
    iter=0;
    while ~stop_criteria && iter < max_iter
        iter = iter + 1;

        % Derivada direccional: proyeccion gradiente sobre Dif
        %g0 = Dif(:).'*g(:);
        if (Dif(:).'*g(:))>=0 %No direccion de descenso
            fprintf(1,'%d: warn: No direccion de descenso\n',iter);
            break;
        end

        gold = g;
        [f, g] = fun(xn_inv,x,xn,SLI); %likelihood term
        
        % Actualizar punto:
        xn_inv = xn_inv + lambda*Dif;

        %gnorm = sum(g(:).^2);
        
        if sqrt(abs(sum(g(:).^2))) < thres * n
            stop_criteria = true;
        end
        
      
        % Gradiente conjugado
        y = g - gold; 

        if strcmp(rule,'PR')
            % Polak-Ribiere-Polyak variant
            bk = g(:)'*y(:)/((gold(:)'*gold(:)) + eps);
        elseif strcmp(rule,'PR')
            % Fletcher-Reeves variant
            bk = g(:)'*g(:)/((gold(:)'*gold(:)) + eps);
        elseif strcmp(rule,'HS')
            % Hestenes-Stiefel variant
            bk=(g(:)'*y(:))/((Difold(:)'*y(:))+eps);
        elseif strcmp(rule,'DY')
            % Day-Yuan variant
            bk=(g(:)'*g(:))/((Difold(:)'*y(:))+eps);
        elseif strcmp(rule,'combined')
            % Combined variant
            bpr = g(:)'*y(:)/((gold(:)'*gold(:)) + eps);
            bfr = g(:)'*g(:)/((gold(:)'*gold(:)) + eps);
            if bpr < -bfr
                bk = -bfr;
            elseif bpr> bfr
                bk = bfr;
            else
                bk = bpr;
            end
        end
        
        Dif = bk*Dif - g; 
        
        Difold=Dif;
        
        %if (mod(iter,40)==0)
        %    fprintf('Iter %d, f = %f, g = %f\n',iter,f,sqrt(gnorm)); 
        %end
    end
    
    xn_inv=permute(xn_inv,[1 2 3 5 4]);
 end

%cost function for the optimization   
function [f, g] = fun(xn_inv, x, xn,SLI)
    % xn  = varargin{1};
    s = size(xn);
    g = zeros(s,'like',xn);
    %s_ = [s(1) s(2) s(3) 1 3];
    x = single(x);
        
    for i=1:s(4)
        %g_=zeros(s_,'like',xn);        
        %xn_ = xn(:,:,:,i,:);
        %xn_inv_ = xn_inv(:,:,:,i,:);
        %interpolation of field
        if SLI==1
        g(:,:,1,i,1) = interpn(x(:,:,1,1),x(:,:,1,2),xn(:,:,1,i,1),xn_inv(:,:,1,i,1),xn_inv(:,:,1,i,2),'linear')-x(:,:,1,1);
        g(:,:,1,i,2) = interpn(x(:,:,1,1),x(:,:,1,2),xn(:,:,1,i,2),xn_inv(:,:,1,i,1),xn_inv(:,:,1,i,2),'linear')-x(:,:,1,2);
        g(:,:,1,i,3) = 0;
        else
        g(:,:,:,i,1) = interpn(x(:,:,:,1),x(:,:,:,2),x(:,:,:,3),xn(:,:,:,i,1),xn_inv(:,:,:,i,1),xn_inv(:,:,:,i,2),xn_inv(:,:,:,i,3),'linear')-x(:,:,:,1);
        g(:,:,:,i,2) = interpn(x(:,:,:,1),x(:,:,:,2),x(:,:,:,3),xn(:,:,:,i,2),xn_inv(:,:,:,i,1),xn_inv(:,:,:,i,2),xn_inv(:,:,:,i,3),'linear')-x(:,:,:,2);
        g(:,:,:,i,3) = interpn(x(:,:,:,1),x(:,:,:,2),x(:,:,:,3),xn(:,:,:,i,3),xn_inv(:,:,:,i,1),xn_inv(:,:,:,i,2),xn_inv(:,:,:,i,3),'linear')-x(:,:,:,3);
        end
        %g(:,:,:,i,:)=g_;
    end
    g(isnan(g))=0;
    f = 0.5*sum(g(:).^2);
end