function runge_kutta_4 ()
  pkg load symbolic
  clc;
  warning('off','all');
  display("------------- MÉTODO RUNGE-KUTTA ORDEN 4 -------------");
  
  syms x y
  f = -x*y + 4*x/y
  y0 = 1;
  a = 0;
  b = 1;
  N = 11;
  
  % Calcula los puntos y genera los pares ordenados
  [x,y] = RK2(f,a,b,y0,N);
  puntos = [];
  for k=1:N
    puntos = [puntos; x(k) y(k)];
  endfor
  
  % Muestra los puntos y el polinomio de interpolación
  display("Puntos:");
  for k=1:N
    display(["\t(",num2str(x(k)) ,",",num2str(y(k)),")"]);
  endfor
  display("Polinomio de interpolación: ");
  polinomio = dd_newton(puntos)
  
  % Grafica los puntos, el polinomio y la solución exacta
  stem(x,y,'blue')
  hold on
  ezplot(polinomio,[a,b],'red')
  
  
end


% Metodo de Runge-Kutta para aproximar la solucion de un problema de valor inicial
    % :param str_funcion: string con la funcion que se debe evaluar
    % :param a: Extremo inferior del intervalo
    % :param b: Extremo superior del intervalo
    % :param n: Cantidad de puntos
    % :param y: Valor de y para el x inicial
    % :return [xv,yv]: Set de puntos (x,y) de la solución dy/dx = y(x)
function [xv,yv] = RK2(f,a,b,y0,N)
  
  h = (b-a)/(N-1);
  xv = a:h:b;
  yv = zeros(1,N);
  
  f1 = matlabFunction(f);
  
  yv(1) = y0;
  
  for k=1:N-1
    k1 = f1(xv(k), yv(k));
    k2 = f1(xv(k)+h/2, yv(k)+h*k1/2);
    k3 = f1(xv(k)+h/2, yv(k)+h*k2/2);
    k4 = f1(xv(k)+h, yv(k)+h*k3);
    yv(k+1) = yv(k) + h/6 * (k1+2*k2+2*k3+k4);
  endfor
  
end


 % Diferencias Divididas de Newton para resolver un problema de interpolacion
    % :param puntos: Los puntos deben ser ingresados en una lista, donde cada punto
    % es una lista, donde la primera posicion es el x, y la segunda posicion es el y. [x, y]
    % :return: Polinomio de interpolacion obtenido
function poli_inter = dd_newton(puntos)
  [n, m] = size(puntos);
  if m ~= 2
    error('No son pares ordenados');
  end 
  x = sym ('x');
  poli_inter = puntos(1,2);
  
  lista_y = puntos(1:n,2)';
  lista_x = puntos(1:n,1)';
  variable=1; 
  con=n-1;
  for i=2:n
    variable=variable*(x -puntos(i-1,1));
    nuevos=[];
    for j=1:con
      dividendo=lista_y(j+1)-lista_y(j);
      divisor=lista_x(j + i - 1)-lista_x(j);
      nuevos = [nuevos (dividendo / divisor)];
    endfor
    con=con-1;
    poli_inter=poli_inter+nuevos(1)*variable;
    lista_y=nuevos;
  endfor
  poli_inter=expand(poli_inter);
end