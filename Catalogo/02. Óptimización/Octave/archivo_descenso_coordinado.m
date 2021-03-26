function archivo_descenso_coordinado
  %Ejemplo Numérico
  %P1: Definir parámetros de entrada
  f='((x-2)^2)+((y+3)^2)+x*y';
  valores=[1,1];
  variables=['x','y'];
  tol=10^-8;
  iterMax=10;
  
  %Otro ejemplo que puede ser utilizado simplemente se quitan los % y se ponen en la otra funcion
  %f='(x^3)+(y^3)+(z^3)-2*x*y-2*x*z-2*y*z'
  %valores=[1,1,1];
  %variables=['x','y','z'];

  %P2: Llamar a la función
  [xk,error,e]=coordinado(f,variables,valores,iterMax,tol)%Los argumentos para la funcion biseccion

end


%Funcion que realiza el metodo de descenso coordinado
%f-> funcion que puede estar cierto numero de variables 
%variables->  Un array o lista con las variables como char que tiene la funcion ejemplo ['x','y','z']
%valores->  Un array o lista con los valores de las variables en el mismo orden ejemplo [1,2,3]
%iterMax = el numero total de iteraciones
%tol = el valor de tolerancia
%Retorna el valor de aproximacion como "xk" y el valor del error final
%Muestra la grafica de error vs iteraciones
function [xk,error,e]=coordinado(f,variables,valores,iterMax,tol)
  pkg load symbolic %Invoca a la libreria symbolic
  nTotal=length(variables); # Se obtiene el largo de la lista
  valores_iniciales=valores; #Copiamos los valores ingresados y se guardan como iniciales
  nuevos_valores=[]; #Lista donde se almacenaran los valores de las nuevas variables
  symbolicas=[];#Lista donde se almacenaran las variables simbolico
  e=[];#Lista donde se almacenaran el valor de los errores
  
  #While que almacenara las variables en symbolicas
  n=1;
  while n<nTotal+1
    symbolicas = [symbolicas, sym(variables(n))]
    n=n+1;
  end
  
  f1=sym(f); %Convierte el texto a simbolico
  f=matlabFunction(f1);%Función f en formato del lenguaje M

  
  #While que llevara la cantidad de iteraciones que tendra el metodo
  #Se dentendra al llegar a la cantidad de iteraciones maximas
  i=0;
  while i<=iterMax   
    nuevos_valores=[];%Se reinicia la lista de nuevos valores para la siguiente iteraciones
    #While encargado de llevar el orden de la variable que quedara como incognita ejemplo p=1 seria x en ['x','y','z']
    p=1;
    while p<nTotal+1
      f2=f;
      n2=1;
      
      #While donde se cambiaran las variables symbolicas a su valor numerico
      while n2<nTotal+1
        clc;
        if n2!=p %Compruebacion para evitar cambiar la incognita en la ecuacion 
          f2=subs(f2,symbolicas(n2),valores(n2));
          n2=n2+1;
        else
          n2=n2+1;
        end
      end
      
      fcn=matlabFunction(f2); %Función f en formato del lenguaje M
      [xmin, fval] = fminsearch (fcn, 0); %Se obtiene el minimo de la funcion 
      valores(p)=xmin;%Se actualiza el valor de la variable que se dejo como incognita
      nuevos_valores=[nuevos_valores, xmin]; %Se agrega dicho valor a los nuevos valores
      p=p+1;
    end
    
    error=norm(nuevos_valores-valores_iniciales)%Se obtiene el valor de error de la funcion
    if error<tol%Se comprueba si error es menor a la tolerancia
      e=[e error];
      break;
    else   
      valores=nuevos_valores;%Se transforman los valores obtenidos a los actuales
      valores_iniciales=valores;
      e=[e error]; %se almacena el error a la lista de errores
      
      if i==iterMax %Evita que se sume un valor extra i 
        break;
      end
      
      i=i+1;
    end  
  end
  
  
  #Se realiza el calculo de la aproximacion
  xk=f;
  t=1;
  #While donde se cambiaran las variables symbolicas a su valor numerico
  while t<nTotal+1
    xk=subs(xk,symbolicas(t),valores(t));
    t=t+1;
  end
  xk=double(xk);%Transforma el valor a numero decimal
  
  clc;
  display(["Iteracion final #:" num2str(i)]);
  title ("Descenso Coordinado error vs iteraciones");
  xlabel("Iteraciones");
  ylabel("Error");
  plot(0:i,e)%Grafica los errores con la cantidad de iteraciones
  
end

