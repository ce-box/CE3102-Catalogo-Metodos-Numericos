function archivo_gauss_seidel
  #Ejemplo para el metodo
  A3=[1 1 5;1 5 1; 5 1 1]
  b3=[7;7;7]
  [x,error]=gauss_seidel(A3,b3,100,10^-3)#A matriz nxn, #b matriz de valores independientes,#k numero de iteraciones,#tol Toleracion

  
end

#Parametros de entrada
#A matriz nxn, #b matriz de valores independientes
#k cantidad de iteraciones, #tol tolerancia 
#Salida
#x matriz  de las soluciones, #error array con los valores de error
function [x,error]=gauss_seidel(A,b,k,tol)
  [n, m] = size(A);#Obtiene la cantidad de filas y columnas
  if n!=m #Comprueba que la matriz es cuadrada
    error("La matriz no es cuadrada");
  endif
  [A,b]=matriz_dominante(A,b,n)#Comprueba que la matriz sea diagonalmente dominante, y si es posible la corrige 
  
  [L,U,D]=matrices_Jacobi(A,n);#Se obtienes las matrices L U D
  x=zeros(n,1);#Crea la matriz que tendra las respuestas de las iteraciones
  x1=x;#Toma el valor de la primera x(0)
  error=[];#Array donde se almacenaran los errores
  H=(L+D);#Matriz provisional para los calculos
  c=hacia_Adelante(H,b,n);#Matriz c => (L+D)c=b
  for i=1:k#Cantidad de iteraciones
    
    d=-U*x(1:n,i);#Matriz provisional para los calculos d(k)=-Ux(k)
    z=hacia_Adelante(H,d,n);# Hz(k)=d(k)
    x1=z+c;#Calcula el valor de las x
    x=[x x1];#Se alamarecen en la matriz de respuesta
    e_norma2=norm(A*x1-b);#Se calcula el error absoluto
    error=[error e_norma2];#Se almacena el error absoluto 
    if e_norma2<=tol#Comprueba el error absoluto, si es menor a la tolerancia se detiene
      break;
    endif
  endfor
  plot(1:i,error);#Se grafica error vs iteraciones
  title ("Gauss Seidel error vs iteraciones")
  xlabel("Iteraciones");
  ylabel("Error");
end


#Parametros de entrada
#A matriz nxn, #n largo del a matriz
#Salida
#L matriz tri.inferior , #U matriz tri.superior, #D matriz diagonal
#Funcion complementaria para obtener las 3 matrices extras
function [L,U,D]=matrices_Jacobi(A,n)#Metodo para obtener las matrices a usar L U D
  L=zeros(n,n);#Matrices nxn llenas de 0
  U=zeros(n,n);
  D=zeros(n,n); 
  for i=1:n#Filas
    for j=1:n#Columnas
      if j==i#Comprueba que sea la diagonal
        D(i,j)=A(i,j);#Almacena la el valor en la matriz diagonal
      elseif i<j
        U(i,j)=A(i,j);#Almacena la el valor en la matriz triangular superior
      elseif  i>j
        L(i,j)=A(i,j);#Almacena la el valor en la matriz triangular inferior
      end
    endfor
  endfor   
end


#Parametros de entrada
#A matriz  triangular inferior, #b nueva matriz de valores independientes, #n largo de la matriz
#Salida
#soluciones matriz de soluciones "x#"
#Sustitucion hacia adelante
function [soluciones]=hacia_Adelante(A,b,n)
  soluciones=zeros(n,1);#Matriz de solucines 
  soluciones(1)=b(1)/A(1,1);#Primera solucion
  i=2;
  while i<n+1#Mueve las filas   
    sumatoria=0;
    for j=1:i-1#Mueve columnas
      sumatoria=sumatoria+(A(i,j)*soluciones(j));#Realiza la sumatoria necesaria
    endfor
    x=((b(i)-sumatoria)/A(i,i));#Calcula el valor de x
    soluciones(i)=x;#Agrega x a la matriz de soluciones
    i=i+1;
  end

end

#Parametros de entrada
#A matriz nxn , #b nueva matriz de valores independientes, #n largo de la matriz
#Salida
#A matriz corregida, #b matriz valores independientes corregida
#Comprueba y corrige la matriz A, a diagonal dominante
function [A,b]=matriz_dominante(A,b,n)
  fila=[];#Almacena el numero de fila donde se detecto un error
  posicion=[]; #Almacena el indice donde se detecto un error
  for f=1:n#Mueve filas
    for c=1:n  #Mueve columnas   
          if abs(A(f,c)) == max(abs(A(f,:))) && f!=c#Comprueba si el numero es igual al mayor y no se encuenta en la diagonal
            fila=[fila f];#Se almacena el numero de fila
            posicion=[posicion c];#Se almacena el indice
            break;      
          endif
    endfor
  endfor
  if isempty(fila)#Se comprueba si la matriz fila esta vacia, si es asi esta es dominante
    display("dominante");
    
  else    #Si no es dominante
    for k=1:size(posicion,2)
      if size(find(posicion==posicion(k)),2)>1   #Comprueba si hay 2 o mas filas que tenga el mayor en la misma posicion
        error("No es posible acomodar la matriz")#si es asi se lanza error porque la matriz no podra corregirse
      endif
    endfor
    display("No es dominante, corrigiendo")
    N1=A(fila(1),:);     #Algortimo que cambia 2 filas tanto de A y b, para ir corregiendo la matriz
    N2=A(posicion(1),:);
    A(posicion(1),:)=N1;
    A(fila(1),:)=N2;
    N3=b(fila(1),:);
    N4=b(posicion(1),:);
    b(posicion(1),:)=N3;
    b(fila(1),:)=N4;
    [A,b]=matriz_dominante(A,b,n);
  end
end