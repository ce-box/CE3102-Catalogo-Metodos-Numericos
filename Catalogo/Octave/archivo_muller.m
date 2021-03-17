function archivo_muller
  %Ejemplo Numérico
  %P1: Definir parámetros de entrada
  f='sin(x)-(x/2)';
  x0=2;
  x1=2.2;
  x2=1.8;;
  iterMax=5;
  tol=10^-8;
  muller(f,x0,x1,x2,iterMax,-0.1);%Argumentos de la funcion muller

end

%Funcion que realiza el metodo de biseccion
%f=funcion
%x0 y x0 = los valores iniciales
%iterMax = el numero total de iteraciones
%tol = el valor de tolerancia
%Muestra la grafica de error vs iteraciones
function muller(f,x0,x1,x2,iterMax,tol)
  pkg load symbolic%Invoca a la libreria symbolic
  syms x%Define a x como simbolica
  f1=sym(f);%convierte el Texto en  simbolico
  f=matlabFunction(f1);%Función f en formato del lenguaje M
  e=[];%La lista donde se almacenaran los valores de error
  con=0;%Contador que ademas llevara el numero de iteraciones
  while(con<iterMax)
    lista1=[];#Matriz 3x3 con que contendra los valores del sistema de ecuaciones
    val=[];#Matriz 3x1 que contendra los valores "resultado" de las ecuaciones
    a=[];#Lista que contendra los valores de la ecuacion 1
    b=[];#Lista que contendra los valores de la ecuacion 2
    c=[0,0,1];#Lista que contendra los valores de la ecuacion 3
    
    a=[a ((x0-x2)^2)];%Calcula y almacena en la lista correspondiente, el valor de las constantes en la ecuacion 1
    a=[a ((x0-x2))];
    a=[a 1];
  
    b=[b ((x1-x2)^2)];%Calcula y almacena en la lista correspondiente, el valor de las constantes en la ecuacion 2
    b=[b ((x1-x2))];
    b=[b 1];
  
    lista1=[lista1;a];%Se agregan a la Matriz principal los valores de las 3 ecuaciones
    lista1=[lista1;b];
    lista1=[lista1;c];
  
    val=[val;f(x0)];%Se calculan y se agregan a la matriz de resultados
    val=[val;f(x1)];
    val=[val;f(x2)];
    
    r=resolver(lista1,val,x0,x1,x2);%Llama a la funcion resolver, la cual le retornara el valor de R
    [t0,t1,t2]=masCercanos(x0,x1,x2,r);%Llama a la funcion masCercanos, la cual le retornara el valor de t0,t1,t2
    x0=t0;
    x1=t1;
    x2=t2;
    
    if(abs(f(r))<tol)%Comprueba si el valor de error es menor a la tolerancia
      break;%Corta el ciclo del while si se cumple la condicion
    else
      e=[e abs(f(r))];%En caso de no cumplir la condicion se agrega el error a la lista e
      con=con+1;%Se aumenta el contador para continuar con el ciclo
    end
  end
    clc;
    display(['-------------'])
    display(['Iteracion #' num2str(con)])
    display(['la aproximacion R=' num2str(r)])
    display(['El error es=' num2str(e)])
    plot(1:con,e);%Se muestra la grafica de ERROR VS Iteraciones de la funcion
    title ("Muller error vs iteraciones");
    xlabel("Iteraciones");
    ylabel("Error");
    
end



%Funcion que obtendra el valor de las variables del sistema de ecuaciones
%matr = array matriz A 3x3
%val = array matriz B 3x1
%x0, x1 y x2 = valores iniciales
%retorna el valor de "r"
function [r]=resolver(matr,val,x0,x1,x2)
  resultado=[];
  resultado=inv(matr)*val;%Operacion de matrices que dara como resultado el valor de las variables del sistema de ecuaciones
  r=obtenerR(resultado(1),resultado(2),resultado(3),x2);%Invoca la funcion "" con la finalidad de obtener el valor de R
end

function [r]=obtenerR(a1,b1,c1,x2)%Funcion encargada de calcular el valor de r para el metodo
  syms x a b c s %Se definen las variables simbolicas
  signo=obtenerSig(b1);%Se llama a la funcion para obtener el simbolo de la variable b
  r='x-((2*c)/(b+s*sqrt((b^2)-4*a*c)))';%Se define como texto la ecuacion para r
  r1=sym(r);%convierte el Texto en  simbolico
  f2=matlabFunction(r1);%Función f en formato del lenguaje M
  r=f2(a1,b1,c1,signo,x2);%(a,b,c,s,x)Orden para ingresar las variables y se obtiene el valor de t=r
end

%Funcion para obtener el simbolo de la variable b
%b = el valor de la variable "b"
%retorna el valor del signo de b, -1 o 1
function [signo]=obtenerSig(b)
  signo=0;
  if b<0%Condicion para determinar si b es negativa
    signo=-1;%Si se cumple la condicion se retornara un -1, que al estar multiplicado similara el signo negativo
  else %Si b es positiva entonces retornara un 1
    signo=1;
  end
end

%Funcion encargada de determinar que valores estan mas cercanos a R
%%x0, x1 y x2 = variables iniciales/iteraciones siguientes
%retorna t0, t1 y t2 = los nuevos valores para las variables a evaluar
function [t0,t1,t2]=masCercanos(x0,x1,x2,r)
  m0=abs(r-x0);%se realiza una resta entre los valores para determinar quien esta mas cerca
  m1=abs(r-x1);
  m2=abs(r-x2);
  if(m0<m2 && m1<m2)%Se evalua que m0 y m1 sean mas pequeños que m2
    t0=x0;%Las variables tomaran el valor de x0,x1 y r
    t1=x1;
    t2=r;
  elseif(m0<m1 && m2<m1)%Se evalua que m0 y m2 sean mas pequeños que m1
    t0=x0;%Las variables tomaran el valor de x0,x2 y r
    t1=x2;
    t2=r;
  elseif(m1<m0 && m2<m0)%Se evalua que m1 y m2 sean mas pequeños que m0
    t0=x1;%Las variables tomaran el valor de x1,x2 y r
    t1=x2;
    t2=r;
  end
end

