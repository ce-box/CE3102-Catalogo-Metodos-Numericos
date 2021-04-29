function archivo_gauss_seidel
  #Ejemplo para el metodo

  
  A=[2 -10 3; 10 3 1; 0 -1 2];
  A2=[10 3 1;2 -10 3; 0 -1 2];
  b= [-5;14;14];
  b2= [14;-5;14];
  
  A3=[1 1 5;1 5 1; 5 1 1]
  b3=[7;7;7]
  [x]=gauss_seidel(A3,b3,5)#A matriz nxn, #b matriz de valores independientes

  
end

function [x]=gauss_seidel(A,b,k)
  [n, m] = size(A);#Obtiene la cantidad de filas y columnas
  if n!=m #Comprueba que la matriz es cuadrada
    error("La matriz no es cuadrada");
  endif
  [A,b]=matriz_dominante(A,b,n)
  
  [L,U,D]=matrices_Jacobi(A,n);
  x=zeros(n,1);
  x1=x;
  H=(L+D);
  c=hacia_Adelante(H,b,n);
  for i=1:k
    d=-U*x(1:n,i);
    z=hacia_Adelante(H,d,n);
    x1=z+c;
    x=[x x1];
    
  endfor

end

function [L,U,D]=matrices_Jacobi(A,n)
  L=zeros(n,n);
  U=zeros(n,n);
  D=zeros(n,n); 
  for i=1:n
    for j=1:n
      if j==i
        D(i,j)=A(i,j);
      elseif i<j
        U(i,j)=A(i,j);
      elseif  i>j
        L(i,j)=A(i,j);
      end
    endfor
  endfor   
end

#Sustitucion hacia adelante
function [soluciones]=hacia_Adelante(A,b,n)
  soluciones=zeros(n,1);
  variable=[];
  soluciones(1)=b(1)/A(1,1);
  i=2;
  while i<n+1   
    sumatoria=0;
    for j=1:i-1
      sumatoria=sumatoria+(A(i,j)*soluciones(j));
    endfor
    x=((b(i)-sumatoria)/A(i,i));
    soluciones(i)=x;
    %variable=[variable strcat("->x",num2str(i),": ",num2str(x),",")];
    i=i+1;
  end

end

function [A,b]=matriz_dominante(A,b,n)
  fila=[];
  posicion=[]; 
  for f=1:n
    for c=1:n     
          if abs(A(f,c)) == max(abs(A(f,:))) && f!=c
            fila=[fila f];
            posicion=[posicion c];
            break;      
          endif
    endfor
  endfor
  if isempty(fila)
    display("dominante");
    
  else    
    for k=1:size(posicion,2)
      if size(find(posicion==posicion(k)),2)>1   
        error("No es posible acomodar la matriz")
      endif
    endfor
    display("No es dominante corrigiendo")
    N1=A(fila(1),:);
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