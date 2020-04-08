function y=BSplineFun(x,n,d)

if nargin<3;d=0;end

% BSPLINEN funci�n B-spline de orden n and its derivatives
%
%Entradas
%x: valores de entrada a modificar
%n: orden del bspline
%
%Salidas
%y: valores bspline correspondientes

xa=abs(x);
y=x;
y(:)=0;
%Non-null values range
li=0.5+0.5*mod(n,2);
x0r=xa<li;
if n>=2;x1r=xa>=li & xa<li+1;end
    
if d==0
    if n==0        
        y(x0r)=1;
    elseif n==1        
        y(x0r)=1-xa(x0r);
    elseif n==2        
        y(x0r)=-x(x0r).^2+0.75;
        y(x1r)=0.5*(xa(x1r)-1.5).^2;
    elseif n==3
        y(x0r)=2/3-xa(x0r).^2+0.5*xa(x0r).^3;
        y(x1r)=(2-xa(x1r)).*(2-xa(x1r)).*(2-xa(x1r))/6;        
    end
elseif d==1
    if n==0
        error('Derivative of BSpline not defined for order 0\n');
    elseif n==1
        y(x0r)=-sign(x(x0r));
    elseif n==2
        y(x0r)=-2*x(x0r);
        y(x1r)=x(x1r)-1.5*sign(x(x1r));
    elseif n==3
        y(x0r)=1.5*sign(x(x0r)).*x(x0r).^2 -2*x(x0r);
        y(x1range)=-sign(x(x1r)).*(x(x1r)-sign(x(x1r))*2).^2/2;
    end
elseif d==2
    if n<2
        error('Derivative of BSpline not defined for order %d\n',n);
    elseif n==2
        y(x0r)=-2;
        y(x1r)=1;
    elseif n==3
        y(x0r)=3*xa(x0r)-2;
        y(x1r)=-xa(x1r)+2;
    end
end
