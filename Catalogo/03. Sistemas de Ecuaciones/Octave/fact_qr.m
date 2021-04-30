function archivo_qr
  
  %Valor de la matriz A
  A = [-1 -1 1; 1 3 3; -1 -1 5; 1 3 7]
       
  %vector de terminos independientes
  term_ind = []
  
  %Se llama a la funcion fact_qr con A y tol
  [Q,R] = fact_qr(A, term_ind)
  
end

function [Q,R] = fact_qr(A, ind)
  
  final = length(A)
  Ar = A
  
  for cont = 1: final
    if cont == 1
      
      vec = A(:,cont)
      norma2 = norm(vec)
      [normCol] = normTrans(norma2, final-cont)
      u = vec - normCol
      
      [iden] = genIden(final-cont, final-cont)
      ut = u.'
      num = 2 * u * ut
      den = ut*u
      iden
      H1 = iden - (num/den)
      
      HA = (H1 * Ar)
      
      A = HA
      
    endif
    
    if cont == 2
      vec = A(:,cont)
      norma2 = norm(vec)
      [normCol] = normTrans(norma2, 3)
      u = vec + normCol
      
      [iden] = genIden(3, 3)
      ut = u.'
      num = 2 * u * ut
      den = ut*u
      Htecho = iden - (num/den)
      
      [H2] = obtH2(final, Htecho)
      
      H2H1A = H2 * H1 * Ar
      
      A = H2H1A
      
    endif
    
    if cont == 3
      
      vec = [A(end-1,end) ; A(end, end)]
      norma2 = norm(vec)
      [normCol] = normTrans(norma2, final-cont)
      u = vec + normCol
      
      [iden] = genIden(1, 1)
      ut = u.'
      num = 2 * u * ut
      den = ut*u
      iden
      Htecho = iden - (num/den)
      
      [H3] = obtH3(final, Htecho)
      R = int32(H3*H2*H1*Ar)
      Q = int32(H1*H2*H3)
      
    endif
    
  endfor
  
end

function [normCol] = normTrans(norm, tot)
  col = zeros(tot+1,1)
  col(1,1) = norm
  
  normCol = col
end

function [H2] = obtH2(tot, techo)
  
  inicial = zeros(tot, tot)
  inicial(1,1) = 1
  
  for j = 1: tot
    for i = 1: tot
      if not(i == 1 || j == 1)
        inicial(j,i) = techo(j-1, i-1)
      endif
    endfor
  endfor
  
  H2 = inicial
  
end

function [H3] = obtH3(tot, techo)
  
  inicial = zeros(tot, tot)
  inicial(1,1) = 1
  inicial(2,2) = 1
  
  for j = 1: tot
    for i = 1: tot
      if not(i == 1 || j == 1) && not(i == 2 || j ==2)
        inicial(j,i) = techo(j-2, i-2)
      endif
    endfor
  endfor
  
  H3 = inicial
  
end
 
function [iden] = genIden(tot)
  I = zeros(tot+1,tot+1)
  for j=1: tot+1
    for i=1: tot+1
      if j == i
        
        I(i,j) = 1
      
      endif
    endfor
  endfor
  iden = I
end