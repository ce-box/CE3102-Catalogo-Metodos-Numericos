function archivo_gausseana
  #Ejemplo para el metodo
  A= [2 -6 12 16; 1 -2 6 6;-1 3 -3 -7;0 4 3 -6]; #Matriz 4x4
  b= [70; 26; -30;-26];
  [x]=gausseana(A,b)#A matriz nxn, #b matriz de valores independientes
  
  A2=[4 -2 1; 20 -7 12; -8 13 17];
  b2= [11;70;17];
  #[x]=gausseana(A2,b2)
end


#Parametros de entrada
#A matriz nxn, #b matriz de valores independientes
#Salida
#x matriz  de las soluciones
function [x]=gausseana(A,b)#Metodo de eliminacion gausseana
  
  [n, m] = size(A);#Obtiene la cantidad de filas y columnas
  if n!=m #Comprueba que la matriz es cuadrada
    error("La matriz no es cuadrada");
  endif
  
  [A,b]=trian_Sup(A,b,n);#Se obtiene la matriz triangular superior y la nueva matriz de valores independientes
  
  if det(A)!=0#Se determina que la nueva matriz sea invertible
    
    [x]=hacia_Atras(A,b,n)#Se calcula la respuesta usando la susticion hacia atras
  else
    display("La matriz no tiene inversa")
  end
end

#Parametros de entrada
#A matriz nxn, #b matriz de valores independientes
#Salida
#A matriz  triangular superior, #b nueva matriz de valores independientes
#Transforma la matriz A a una triangular superior
function [A,b]=trian_Sup(A,b,n)
  M=[A b];
  for k=1: n-1#Columnas
    for i=k+1:n#Filas
      Mik=M(i,k)/M(k,k);#Valor      
      for j=k:n+1
        M(i,j)=M(i,j)-Mik*M(k,j);#Actualiza valores    
      endfor
    endfor
  endfor
  A=M(:,1:n);
  b=M(:,n+1);
end

#Parametros de entrada
#A matriz  triangular superior, #b nueva matriz de valores independientes
#Salida
#soluciones matriz de soluciones "x#"
#Susticion hacia atras
function [soluciones]=hacia_Atras(A,b,n)
  soluciones=zeros(n,1);#Lista de los valores de X la solucion
  i=n;
  while i>0   #filas
    sumatoria=0;
    for j=i+1:n#Columnas 
      sumatoria=sumatoria+(A(i,j)*soluciones(j));#Calcula la sumatoria de la funcion
    endfor
    x=(1/A(i,i))*(b(i)-sumatoria);#Termina el calculo de la X
    soluciones(i)=x;#Se agrega a la matriz de solucion
    i=i-1;
  end

  
end