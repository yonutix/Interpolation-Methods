function N = eval_interpolator_d(tip, eps)
    N=0;
    x = linspace( 0,300, 301 );
    x(1)=1;% x devine indicele fiecarui an
    
    nk=3;
    for i=4:150%calculez divizorii lui 300 
       if(mod(300,i)==0) 
          nk=[nk,i]; 
       end
    end
    h = 2*pi/300;%constanta din enunt
    sum = 0; ok = 0;%presupunem ca e convergenta;0 inseamna convergenta
    y = linspace( 0, 300, nk(2)+1 );%impart intervalul in nk(2)+1 puncte in care voi calcula polinomul de interpolare
    y(1)=1;%altfel y(1) ar fi 0
    %Lagrange
    if( tip == 1 )%daca se alege polinomul de interpolare Newton
        for i = 1:300
            sum = sum + abs(f(i) - lagrange( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 4 puncte de interpolare
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - lagrange( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - lagrange( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
    
%Newton
    if( tip == 2 )%daca se alege polinomul de interpolare Newton
        for i = 1:300
            sum = sum + abs(f(i) - newton( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 4 puncte de interpolare
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - newton( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - newton( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
    
%Spline liniar
    if( tip == 3 )%daca se alege spline liniar
        for i = 1:300
            sum = sum + abs(f(i) - liniar( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - liniar( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare

        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - liniar( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end

%Spline C2 Natural
    if( tip == 4 )%daca se alege spline-ul cubic natural
        for i = 1:300
            sum = sum + abs(f(i) - splinec2( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - splinec2( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - splinec2( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end

    %Spline cubic clasa C2 tensionat
    if( tip == 5 )%daca se alege spline-ul cubic Tensionat
        for i = 1:300
            sum = sum + abs(f(i) - splineTensionat( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - splineTensionat( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - splineTensionat( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
 end


%aproximare trigonometrica
    if( tip == 6 )%daca se alege polinomul de interpolare trigonometric 
        for i = 1:300
            sum = sum + abs(f(i) - trigonometrie( y(1:nk(2)+1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        
        y = linspace( 0, 300, nk(3)+1 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        y(1)=1;
        sum = 0;
        for i = 1:300
            sum = sum + abs( f(i) - trigonometrie( y(1:nk(3)+1), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        n=3;
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
               n=n+1;
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( 0, 300, nk(n)+1 );%se genereaza puncte de interpolare pentru noul Nk2
                y(1)=1;
                sum = 0;
                for i = 1:300
                    sum = sum + abs( f(i) - trigonometrie( y(1:nk(n)+1), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = n;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
end

function P = lagrange( x ,x0)
    n = length(x);
    P = 0;
    for k=1:n%suma de la 1 la k
        p1=1;%initializare produs numarator
        for i=1:n
            if (i~=k)
                p1=p1*(x0-x(i));
            end
        end%efectuare produs numarator
        p2=1;%initializare produs numitor
        for i=1:n
            if(i~=k)
                p2=p2*(x(k)-x(i));
            end%efectuare produs numitor
        end
        P=P+f(x(k))*p1/p2;%efectuare suma finala
    end
end

function b = newton(x,x0)
	 n=length(x);%lungimea vectorului x
	 c=difdiv(x,f(x)');
	 b=c(1);
	 p=1;
	for i=2:n
		p=p*(x0-x(i-1));
		b=b+p*c(i);
    end

end

function a=difdiv(x,y)
	n = length( x );
	for k = 1 : n-1
		y( k+1 : n ) = ( y( k+1 : n ) - y( k ) ) ./ ( x( k+1 : n ) - x( k ) );
    end
	a = y( : );
end

function rez = liniar( x,x0 )
    n=length(x);
    a=zeros(1,n-1);
    b=a;
    for i=1:n-1
        if (x0 <= x(i+1) && x0 >= x(i))
            a=(f(x(i+1))-f(x(i)))/(x(i+1)-x(i));
            b=(x(i+1)*f(x(i))-x(i)*f(x(i+1)))/(x(i+1)-x(i));
            break;
        end
    end
    rez=a*x0+b;
end

function rez = splinec2(x,x0)
    n = length(x);
    a = zeros(1,n);
    b = zeros(1,n-1);
    h = zeros(1,n-1);
    M = zeros(n,n);
    N = zeros(n,1);
    for i = 1:n
        a(i) = f(x(i));
    end
    for i = 1:n-1
        h(i) = x(i+1)-x(i);
    end
    for i = 2:n-1
        M(i,i-1) = h(i-1);
        M(i,i+1) = h(i);
        M(i,i) = 2*(h(i-1)+h(i));
    end
    M(1,1) = h(1);
    M(n,n) = h(n-1);
    for i = 2:n-1
        N(i,1) = (3*(a(i+1)-a(i))/h(i))-(3*(a(i)-a(i-1))/h(i-1));
    end
    c = inv(M)*N;
    for i = 2:n-1
        b(i) = (a(i+1)-a(i))/h(i) - (h(i)/3)*(2*c(i)+c(i+1));
    end
    b(1) = (a(2)-a(1))/h(1)-(h(1)/3)*(2*c(1)+c(2));
    for i = 1:n-1
        d(i) = (c(i+1)-c(i))/(3*h(i));
    end
    for i = 1:n-1
        if(x0 <= x(i+1) && x0 >= x(i))
           rez = a(i) + b(i) * (-x(i)+x0)+c(i)*(-x(i)+x0)^2+d(i)*(-x(i)+x0)^3; 
        end
    end
end

function rez = splineTensionat(x,x0)%aplic algoritmul si formlee din curs/seminar...
    n = length(x);
    a = zeros(1,n);
    b = zeros(1,n-1);
    c = zeros(1,n);
    M = zeros(n,n);
    N = zeros(n,1);
    fd1 = (f(x(2))-f(x(1)))/(x(2)-x(1));
    fdn = (f(x(n))-f(x(n-1)))/(x(n)-x(n-1));
    for i = 1:n
        a(i) = f(x(i));
    end
    h = x(2)-x(1);
    M(1,1) = 2*h;
    for i=2:n
        M(i,i-1) = h;
        M(i-1,i) = h; 
        if(i ~= n)
            M(i,i) = 4*h;
        end
    end
    M(n,n) = 2*h;
    for i = 2:n-1
        N(i,1) = (3*(a(i+1)-a(i))/h)-(3*(a(i)-a(i-1))/h);
    end
    N(1) = 3*(a(2)-a(1))/h-3*fd1;
    N(n) = 3*fdn-3*(a(n)-a(n-1))/h;
    c = inv(M)*N;
    for i = 1:n-1
        b(i) = (a(i+1)-a(i))/h - (h/3)*(2*c(i)+c(i+1));
    end
    for i = 1:n-1
        d(i) = (c(i+1)-c(i))/(3*h);
    end
    for i = 1:n-1
        if(x0 <= x(i+1) && x0 >= x(i))
           rez = a(i)+b(i)*(-x(i)+x0)+c(i)*(-x(i)+x0)^2+d(i)*(-x(i)+x0)^3; 
        end
    end
end

function rez = trigonometrie( x,x0 )
n=length(x);
 xi = linspace(-pi,pi,2*n+1);
            yi = f(1:300);
            sum1 = 0;
            for i=1:2*n+1
                sum1 = sum1 + yi(i);
            end;
           a0 = sqrt(2)/(2*n+1)*sum1;
            for i = 1:n
                sum1 = 0;
               sum2 = 0;
                for j=1:2*n+1
                    sum1 = sum1 + yi(j)*sin(xi(j)*i);
                    sum2 = sum2 + yi(j)*cos(xi(j)*i);
                end;
                b(i) = 2/(2*n+1)*sum1;
                a(i) = 2/(2*n+1)*sum2;
            end;
                sumsin = 0;
                sumcos = 0;
                for j=1:n
                    sumsin = sumsin + b(j)*sin(j*x0);
                    sumcos = sumcos + a(j)*cos(j*x0);
                end;
                sumt = a0/sqrt(2) + sumsin + sumcos;

rez=sumt;         
end

function rez = f( i )%valorile functie pot fi calculate doar in 300 de puncte
load sunspot.dat;
y=sunspot(1:300,2);
    rez=y(i);
end

