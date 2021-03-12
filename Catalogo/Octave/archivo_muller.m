function archivo_muller
  %Ejemplo Numérico
  %P1: Definir parámetros de entrada
  f='sin(x)-(x/2)';
  x0=2;
  x1=2.2;
  x2=1.8;;
  iterMax=5;
  
  muller(f,x0,x1,x2,iterMax,-0.1);

end


function muller(f,x0,x1,x2,iterMax,tol)
  pkg load symbolic
  syms x
  f1=sym(f);%convierte el Texto en  simbolico
  f=matlabFunction(f1);
  e=[]
  con=0;
  while(con<iterMax)
    lista1=[]
    val=[];
    a=[];
    b=[];
    c=[0,0,1];
    a=[a ((x0-x2)^2)];
    a=[a ((x0-x2))];
    a=[a 1];
  
    b=[b ((x1-x2)^2)];
    b=[b ((x1-x2))];
    b=[b 1];
  
    lista1=[lista1;a];
    lista1=[lista1;b];
    lista1=[lista1;c];
  
    val=[val;f(x0)]
    val=[val;f(x1)]
    val=[val;f(x2)]
    
    r=resolver(lista1,val,x0,x1,x2)
    [t0,t1,t2]=masCercanos(x0,x1,x2,r)
    x0=t0;
    x1=t1;
    x2=t2;

    
    if(abs(f(r))<tol)
      break;
    else
      e=[e abs(f(r))]
      con=con+1;
    end
    clc;
  end
    clc;
    display(['-------------'])
    display(['la aproximacion R=' num2str(r)])
    display(['la aproximacion R=' num2str(con)])
    display(['la aproximacion R=' num2str(e)])
    plot(1:con,e);%Se muestra la grafica de ERROR VS Iteraciones de la funcion
    title ("Muller error vs iteraciones");
    xlabel("Iteraciones");
    ylabel("Error");
end
function [r]=resolver(matr,val,x0,x1,x2)
  resultado=[];
  resultado=inv(matr)*val;
  r=obtenerR(resultado(1),resultado(2),resultado(3),x2)
  
end

function [signo]=obtenerSig(b)
  signo=0;
  if b<0
    signo=-1;
  else
    signo=1;
  end
end

function [r]=obtenerR(a,b,c,x2)
  signo=obtenerSig(b)
  r=x2-((2*c)/(b+signo*sqrt((b^2)-4*a*c)))
end

function [t0,t1,t2]=masCercanos(x0,x1,x2,r)
  m0=abs(r-x0);
  m1=abs(r-x1);
  m2=abs(r-x2);
  if(m0<m2 && m1<m2)
    t0=x0;
    t1=x1;
    t2=r;
  elseif(m0<m1 && m2<m1)
    t0=x0;
    t1=x2;
    t2=r;
  elseif(m1<m0 && m2<m0)
    t0=x1;
    t1=x2;
    t2=r;
  end
end
  