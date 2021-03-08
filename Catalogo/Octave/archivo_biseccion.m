
function archivo_biseccion
  %Ejemplo Numérico
  %P1: Definir parámetros de entrada
  f='exp(x)-x-2';
  %f='exp(-x)-x'; %a=0 b=1 tol=1*10^-3;
  a=0;
  b=2;
  tol=10^-4;
  iterMax=100;
  %P2: Llamar a la función
  [xk,error]=biseccion(f,a,b,iterMax,tol)
end

function[xk,error]=biseccion(f,a,b,iterMax,tol)
  pkg load symbolic %Invoca a la libreria symbolic
  syms x
  f1=sym(f); %Convierte el texto a simbolico
  f=matlabFunction(f1);%Función f en formato del lenguaje M
  k=1;%Se define iteracion inicial
  e=[];%Lista que contendra los valores de error
  for i=1:iterMax
    xk=(a+b)/2; %Calculo y actualizacion de la aproximacion
    fx=f(xk);%Valor de la funcion con respecto a la aproximacion
    error=(b-a)/2^k;%Se calcula el error
    display(['---------------']);
    display(['iteracion #' num2str(k)]);
    display(['a:' num2str(a)]);
    display(['b:' num2str(b)]);
    e=[e error];%Se agrega a la lista el valor de error actual
    k=k+1;%Actualiza iteraciones
    if f(a)*fx<0 %Se verifica la condicion de f(a)*f(x) sea menor a 0  
       b=xk;%Se actualiza el valor de b
    elseif fx*f(b)<0%Sino verifica la condicion de f(b)*f(x) sea menor a 0  
       a=xk;%Se actualiza el valor de b
    end
    if abs(error)<tol%Condicion de parada
      break;
    end
  endfor
  plot(2:k,e);%Se muestra la grafica de ERROR VS Iteraciones de la funcion
  title ("Biseccion error vs iteraciones");
  xlabel("Iteraciones")
  ylabel("Error")
end
