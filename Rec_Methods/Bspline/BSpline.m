function y=BSpline(x,n)
% BSPLINEN funci�n B-spline de orden n
%
%Entradas
%x: valores de entrada a modificar
%n: orden del bspline
%
%Salidas
%y: valores bspline correspondientes

%seleccion del orden por parametro
if(n==1)
    y = bspline1(x);
elseif(n==2)
    y = bspline2(x);
elseif(n==3)
    y = bspline3(x);
elseif(n==0)
    y = bspline0(x);
end

end

