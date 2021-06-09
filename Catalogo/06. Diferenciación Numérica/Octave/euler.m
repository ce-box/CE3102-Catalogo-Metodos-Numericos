function archivo_euler
  clc; clear;
  pkg load symbolic
  warning('off', 'all');
  display("Ejemplo para metodo euler:")
  f='y-(x^2)+1' 
  intervalo=[0,2]
  yinicial=0.5
  pasosh=11
  display("\nSe obtiene:\n ")
  [pares,poli]=euler(f,intervalo,yinicial,pasosh)#Ejemplo de ingreso de datos

end


%Funcion que realiza el metodo de Euler
%Parametros de entrada
%f -> funcion a evaluar, intervalo ->Valor del rango donde se evaluara el metodo
%yinicial -> Valor inicial de 'y', pasosh -> cantidad de puntos
%Parametros de salida
%pares -> Pares x y, %poli -> polinomio de interpolacion es tipo symbolic
function [pares,poli]=euler(f,intervalo,yinicial,pasosh)
  a=intervalo(1);#Separto los extremos
  b=intervalo(2);
  f1=matlabFunction(sym(f));#Se transforma la funcion a tipo matlab
  h=(b-a)/(pasosh-1);#Se calcula el valor de h
  xk=[a];#Inicializa la lista donde se almacenara los valores de x
  yk=[yinicial];#Inicializa la lista donde se almacenara los valores de xy
  x=a;
  for i=1:pasosh-1  
    y=yk(i)+h*f1(x,yk(i));#Se calcula la y siguiente
    x=x+h;#Se calcula la x siguiente
    xk=[xk x];#Se guarda el valor de x_siguiente en la lista xk
    yk=[yk y];#Se guarda el valor de y_siguiente en la lista yk 
  end
  
  pares=[xk' yk'];#Se crea la lista de pares ordenados
  poli=dd_newton(pares);#Se crea el polinomio de interpolacion por el metodo de Diferencias Divididas de Newton
  poli_g=matlabFunction(poli);#Se vuelve tipo matlab
  ezplot(poli_g,intervalo);#Se hace el plot del polinomio de interpolacion en el invertavalo 
  hold on
  stem(xk,yk)#Se colocan los puntos del polinomio en la grafica
end

%Funcion que realiza el metodo de Diferencias Divididas de Newton
%Parametros de entrada
%puntos -> una matriz (m x 2), sera un matriz que en la columna 1 tenga los valores x, en columna 2 las y 
%Parametros de salida
%poli_inter polinomio de interpolacion es tipo symbolic
function poli_inter = dd_newton(puntos)
  [n, m] = size(puntos);
  if m ~= 2#Comprueba que la lista tenga solo pares ordenados
    error('No son pares ordenados');
  end 
  x = sym ('x');#Define variable symbolic
  poli_inter = puntos(1,2);#Se almacena la primera diferencia dividida
  #resultado = dd_newton2(puntos)
  lista_y = puntos(1:n,2)';
  lista_x = puntos(1:n,1)';
  variable=1; #
  con=n-1;
  for i=2:n
    variable=variable*(x -puntos(i-1,1));#Se calcula la variable que sera multiplicado por las diferencias
    nuevos=[];#Se almacenara los multiplicadores dividendo/divisor
    for j=1:con
      dividendo=lista_y(j+1)-lista_y(j);#Se calcula el dividendo ejemplo: F[x1,x2]-[F[x0,x1]
      divisor=lista_x(j + i - 1)-lista_x(j);#Se calcula el divisor ejemplo x2-x0
      nuevos = [nuevos (dividendo / divisor)];#Lista de nuevos resultados
    endfor
    con=con-1;#Se disminuye el numero de diferencias divididas a calcular
    poli_inter=poli_inter+nuevos(1)*variable;#Se agrega el siguiente valor de polinomo de interpolacion
    lista_y=nuevos;#Se actualizan la lista de valores de y
  endfor
  poli_inter=expand(poli_inter);#Combina las variables simbolicas ejemplo x*x -> x^2
end


