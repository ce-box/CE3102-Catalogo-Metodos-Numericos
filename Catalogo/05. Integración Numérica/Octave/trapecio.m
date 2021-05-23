function archivo_trapecio
  clc;
  pkg load symbolic
  warning('off', 'all');
  f='ln(x)';
  intervalo=[2,5];
  [error,aprox]=trapecio(f,intervalo)
  
end
function [error,aprox]=trapecio(f,intervalo)
  f=sym(f);
  f1=matlabFunction(f);
  h=intervalo(2)-intervalo(1); 
  aprox=(h/2)*(f1(intervalo(2))+f1(intervalo(1)));
  error=cota_error_trapecio(f,intervalo);
end
function [error]=cota_error_trapecio(f,intervalo)
  a=intervalo(1);#Separto los extremos
  b=intervalo(2);
  f=sym(f);
  g=diff(diff(f,'x'));  
  fd = diff(g,'x')==0;#Calculo la primera derivada con respecto a x
  puntos_criticos = double(cell2mat(solve(fd,'x')));#Obtenego los puntos criticos  
  puntos_a_evaluar=[a b puntos_criticos];#Obtengo los puntos a evaluar en la funcion 
  f1=matlabFunction(g);#Obtengo la funcion en tipo matlab
  valores_evaluados= [f1(puntos_a_evaluar)];#Obtengo array con los valores evaluados en la funcion
  valores_evaluados=abs(valores_evaluados);#Aplica el valor absoluto
  [fmax]=max(valores_evaluados);#Se obtiene el valor maximo de F
  error=(((b-a)^3)/12)*fmax;
end