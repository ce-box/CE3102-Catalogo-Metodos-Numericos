function archivo_dd_newton
  pkg load symbolic
  warning('off', 'all');
  puntos = [-2 0; 0 1; 1 -1];#Ejemplo de ingreso de pares ordenados
  puntos2 = [0 1; 2 0; 3 4;4 0; 5 5];
  display("Difencias DIvididas de Newton")
  display("Ejemplo 1:")
  display(["Puntos: (",num2str(puntos(1,1)),",",num2str(puntos(1,2)),")",",(",num2str(puntos(2,1)),",",num2str(puntos(2,2)),")",",(",num2str(puntos(3,1)),",",num2str(puntos(3,2)),")"])
  display("Polinomio de interpolación:")
  poli_inter = dd_newton(puntos)
  
  display("Ejemplo 2:")
  display(["Puntos: (",num2str(puntos2(1,1)),",",num2str(puntos2(1,2)),")",",(",num2str(puntos2(2,1)),",",num2str(puntos2(2,2)),")",",(",num2str(puntos2(3,1)),",",num2str(puntos2(3,2)),")",",(",num2str(puntos2(4,1)),",",num2str(puntos2(4,2)),")",",(",num2str(puntos2(5,1)),",",num2str(puntos2(5,2)),")"])
  display("Polinomio de interpolación:")
  poli_inter2 = dd_newton(puntos2)
  
end



%Funcion que realiza el metodo de Diferencias Divididas de Newton
%Parametros de entrada
%puntos -> una matriz (m x 2), sera un matriz que en la columna 1 tenga los valores x, en columna 2 las y 
%Parametros de salida
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

