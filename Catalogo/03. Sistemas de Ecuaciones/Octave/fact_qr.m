function archivo_qr
  
  %Valor de la matriz A
  A = [-1 -1 1; 1 3 3; -1 -1 5; 1 3 7]
       
  %vector de terminos independientes
  term_ind = ['x' 'y' 'z' ]
  
  %Se llama a la funcion fact_qr con A y tol
  [Q,R] = fact_qr(A, term_ind)
  
end

function [Q,R] = fact_qr(A, ind)
  
  final = length(A) %Cantidad de filas
  Ar = A %BackUp de A
  
  for cont = 1: final %Loop
    if cont == 1 %Obtiene H1
      
      vec = A(:,cont) %Obtiene el vector con la primera columna de A
      norma2 = norm(vec) %Calcula la norma2 del vector 
      [normCol] = normTrans(norma2, final-cont) %Se crea una columna con la normal
      u = vec - normCol %Calcula u
      
      [iden] = genIden(final-cont, final-cont) %Calcula identidad para matriz final-cont
      ut = u.' %Calcula transpuesta de u
      num = 2 * u * ut %Calcula 2*u*ut
      den = ut*u %Calcula ut*u
      H1 = iden - num/den %Calcula H1
      
      HA = (H1 * A) %Calcula HA = H1*A
      
      A = HA %Reasigna A
      
    endif
    
    
    if cont == 2 %Obtiene H2
      
      vec = A(:,cont) %Obtiene vector 
      norma2 = norm(vec) %Calcula la norma2 del vector
      [normCol] = normTrans(norma2, 3) %Hace columna con la norma
      u = vec + normCol %Calcula u
      
      [iden] = genIden(3, 3) %Genera matriz identidad 3x3
      ut = u.' %Calcula transpuesta de u 
      num = 2 * u * ut %Calcula 2 * u * ut
      den = ut*u %Calcula ut * u
      Htecho = iden - (num/den) %Calcula H2techo con lo antes obtenido
      
      [H2] = obtH2(final, Htecho) %Genera H2 apartir de H2techo
      
      H2H1A = H2 * H1 * Ar %Genera H2 * H1 * A
      
      A = H2H1A %Reasigna A
      
    endif
    
    if cont == 3 %Genera Q y R apartir de H1 H2 H3
      
      vec = [A(end-1,end) ; A(end, end)] %Genera nuevo vector
      norma2 = norm(vec) %Calcula la norma2 del vector
      [normCol] = normTrans(norma2, final-cont) %Crea columna con la norma2
      u = vec + normCol %Calcula u
      
      [iden] = genIden(1, 1) %Genera la matriz identidad 2x2
      ut = u.' %Calcula transpuesta de 2
      num = 2 * u * ut %Calcula 2*u*ut
      den = ut*u %Calcula ut * u
      Htecho = iden - (num/den) %Calcula H3techo a partir de lo obtenido
      
      [H3] = obtH3(final, Htecho) %Genera H3
      R = int32(H3*H2*H1*Ar) %Calcula R
      Q = H1*H2*H3 %Calcula Q
      
    endif
    
  endfor
  
end

function [normCol] = normTrans(norm, tot) %Genera una columna con la norma2
  
  col = zeros(tot+1,1) %Crea columna de ceros
  col(1,1) = norm %Anade la norma a la posicion 11
  
  normCol = col %Asigna columna
end

function [H2] = obtH2(tot, techo) %Retorna H2 apartir de H2 techo
  
  inicial = zeros(tot, tot) %Matriz de ceros
  inicial(1,1) = 1 %Anade un 1 a la posicion 11
  
  for j = 1: tot %Recorre filas
    for i = 1: tot %Recorre columnas
      if not(i == 1 || j == 1) %Verifica que no este en f1 ni c1
        inicial(j,i) = techo(j-1, i-1) %Asigna valor de H2techo
      endif
    endfor
  endfor
  
  H2 = inicial %Asigna H2 a la salida
  
end

function [H3] = obtH3(tot, techo) %Retorna H3 apartir de H3techo
  
  inicial = zeros(tot, tot) %Matriz de ceros
  inicial(1,1) = 1 %Anade un 1 a la posicion 11
  inicial(2,2) = 1 %Anade un 1 a la posicion 22
  
  for j = 1: tot %Recorre filas
    for i = 1: tot %Recorre columnas
      if not(i == 1 || j == 1) && not(i == 2 || j ==2) %Verifica que no este en f1 ni c1
        inicial(j,i) = techo(j-2, i-2) %Asigna valor de H2techo
      endif
    endfor
  endfor
  
  H3 = inicial %Asigna H3 a la salida
  
end
 
function [iden] = genIden(tot) %Genera matriz identidad
  I = zeros(tot+1,tot+1) %Genera matriz de ceros 
  for j=1: tot+1 %Recorre filas
    for i=1: tot+1 %Recorre columnas
      if j == i %Cuando j = i
        
        I(i,j) = 1 %Inserta un 1
      
      endif
    endfor
  endfor
  iden = I %Asigna iden a la salida
end