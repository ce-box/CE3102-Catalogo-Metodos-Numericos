function archivo_descenso_coordinado
  %Ejemplo Numérico
  %P1: Definir parámetros de entrada
  %f='((x-2)^2)+((y+3)^2)+x*y';
  %valores=[1,1];
  %variables=['x','y'];
  f='(x^3)+(y^3)+(z^3)-2*x*y-2*x*z-2*y*z'
  valores=[1,1,1];
  variables=['x','y','z'];
  tol=10^-8;
  %iterMax=9;
  iterMax=10;
  %P2: Llamar a la función
  [xk,error]=coordinado(f,variables,valores,iterMax,tol)

end


function [xk,error]=coordinado(f,variables,valores,iterMax,tol)
  pkg load symbolic %Invoca a la libreria symbolic
  nTotal=length(variables);
  valores_iniciales=valores;
  nuevos_valores=[];
  symbolicas=[];
  n=1;
 
  e=[];
  while n<nTotal+1
    symbolicas = [symbolicas, sym(variables(n))]
    n=n+1;
  end
  f1=sym(f); %Convierte el texto a simbolico
  f=matlabFunction(f1);%Función f en formato del lenguaje M
  i=0;
  while i<iterMax+1
    p=1;
    while p<nTotal+1
      f2=f;
      n2=1;
      while n2<nTotal+1
        clc;
        if n2!=p
          f2=subs(f2,symbolicas(n2),valores(n2));
          n2=n2+1;
        else
          n2=n2+1;
        end
      end
      fcn=matlabFunction(f2);
      [xmin, fval] = fminsearch (fcn, 0);
      valores(p)=xmin;
      nuevos_valores=[nuevos_valores, xmin];
      p=p+1;
    end
    
    error=norm(nuevos_valores-valores_iniciales)
    if error<tol

      break
    else
      valores=nuevos_valores;
      valores_iniciales=valores;
      nuevos_valores=[];
      e=[e error];
      i=i+1;
    end  
  end
  xk=f;
  t=1;
  while t<nTotal+1
    xk=subs(xk,symbolicas(t),valores(t));
    t=t+1;
  end
  xk=double(xk);
  clc;
  display(["Iteracion final #:" num2str(i-1)]);
  title ("Descenso Coordinado error vs iteraciones");
  xlabel("Iteraciones");
  ylabel("Error");
  plot(1:i,e)
  
end

