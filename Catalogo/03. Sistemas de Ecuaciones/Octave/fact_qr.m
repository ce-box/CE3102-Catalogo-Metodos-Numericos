function archivo_qr
  
  %Valor de la matriz A
  A = [12 -51 4; 6 167 -68; -4 24 -41]
       
  %Valor de la tolerancia
  tol = 10 ^ -10
  
  %Se llama a la funcion fact_qr con A y tol
  [x] = fact_qr(A, tol)
  
end

function [x] = fact_qr(A, tol)
  
  len_A = size(A) %Se guarda el tamano de A
  last_A = zeros(len_A) %Se crea el valor anterior de A 

  while 1
    
      [Q, R] = qr(A) 
      last_A = A % Reasigna el ultimo A calculado
      A = R * Q %Reasigna A actual
      norm1 = norm(last_A - A) %Calcula la norma
      
      %Condicion de parada
      if norm1 < tol
          break;
      end
  end
  
  x = A
  
end