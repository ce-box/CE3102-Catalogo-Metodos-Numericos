function archivo_qr
  
  %Valor de la matriz A
  A = [12 -51 4; 6 167 -68; -4 24 -41]
       
  %vector de terminos independientes
  term_ind = []
  
  %Se llama a la funcion fact_qr con A y tol
  [x] = fact_qr(A, term_ind)
  
end

function [x] = fact_qr(A, ind)
  
  
  
end