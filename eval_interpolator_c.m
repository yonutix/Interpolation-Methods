function N = eval_interpolator_c(tip, eps)
    %sectiune generala pentru toate cazurile
    x = linspace( -pi, pi, 1001 );%punctele in care se testeaza convergenta 
    Nk1 = 4;%numarul de puncte in care se va calcula prima data eroarea
    Nk2 = 8;%numarul de puncte in care se va calcula a 2-a oara eroarea
    h = 2 * pi / 1001;%constanta
    sum = 0; ok = 0;%presupunem ca e convergenta;0 inseamna convergenta
    y = linspace( -pi, pi, Nk1 );%punctele in care se va intepola prima data 
    ok=0;%presupunem metoda convergenta
    
    %Lagrange
    if( tip == 1 )%daca se alege polinomul de interpolare Lagrange 
        for i = 1:1001
            sum = sum + abs(  f(x(i)) - lagrange( y(1:Nk1), x(i))  )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end%insumarea erorilor in ficare punct da eroarea generala
        E1 = sqrt( h * sum );%eroarea pentru 4 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - lagrange( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end%a 2-a suma a erorilor
        E2 = sqrt( h * sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2 * Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i)) - lagrange( y(1:Nk2), x(i) ) )^2; 
                end%calculul erorii pentru Nk2 puncte de interpolare
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end

%newton
    if( tip == 2 )%daca se alege metoda Newton
        for i = 1:1001
            sum = sum + abs( f(x(i)) - newton( y(1:Nk1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h * sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - newton( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h * sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2 * Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i)) - newton( y(1:Nk2), x(i) ) )^2; 
                end
                E2 = sqrt( h * sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
    
%Spline liniear
    if( tip == 3 )%daca se alege spline linear
        for i = 1:1001
            sum = sum + abs( f(x(i)) - liniar(y(1:Nk1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end 
        E1 = sqrt( h * sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - liniar( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h * sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2*Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i) ) - liniar( y(1:Nk2), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
    
%Spline C2 Natural
    if(  tip == 4 )%daca se alege spline-ul natural
        for i = 1:1001
            sum = sum + abs(f(x(i)) - splinec2( y(1:Nk1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h * sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - splinec2( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2*Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i) ) - splinec2(y(1:Nk2), x(i)) )^2; 
                end
                E2 = sqrt( h * sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end
    %Spline cubic clasa C2 tensionat
    if( tip == 5 )%daca se alege polinomul de interpolare Lagrange 
        for i = 1:1001
            sum = sum + abs(f(x(i)) - splineTensionat( y(1:Nk1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - splineTensionat( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2*Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i) ) - splineTensionat( y(1:Nk2), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
        else N = inf;
        end
    end


%interpolare trigonometrica(in enunt se cerea aproximare dar nu prea avea sens asa ca am facut interpolare)
   if( tip == 6 )%daca se alege polinomul de interpolare Lagrange 
        for i = 1:1001
            sum = sum + abs(f(x(i)) - trigonometrie( y(1:Nk1), x(i)) )^2; %calculeaza suma din formula dupa care se calculeaza eroarea pentru 2 puncte de interpolare
        end
        E1 = sqrt( h*sum );%eroarea pentru 2 puncte de interpolare
        y = linspace( -pi, pi, Nk2 );%puncte echidistante generate pentru interpolarea in Nk2 puncte
        sum = 0;
        for i = 1:1001
            sum = sum + abs( f(x(i) ) - trigonometrie( y(1:Nk2), x(i)) )^2; %calculeaza suma din formula pentru eroarea generata de polinomul penttru Kn2 puncte de interpolare
        end
        E2 = sqrt( h*sum );%eroarea pentru Nk2 puncte  de interpolare
        if(E2 > E1) ok = 1; end %daca erorile nu sunt descrescatoare atunci metoda nu converge  
            while( abs( E2 - E1 ) > eps && E2 < E1)%cat timp valoarea absoluta a diferentei si cat timp sunt descrescatoare erorile...
                Nk1 = Nk2;
                Nk2 = 2*Nk2;%puterea lui 2 creste cu o unitate...
                E1 = E2;%^metoda e convergenta pentru Nk1 puncte de interpolare asa ca se incearca  Nk1*2
                y = linspace( -pi, pi, Nk2 );%se genereaza puncte de interpolare pentru noul Nk2
                sum = 0;
                for i = 1:1001
                    sum = sum + abs( f(x(i) ) - trigonometrie( y(1:Nk2), x(i) ) )^2; 
                end
                E2 = sqrt( h*sum );%^aplic formula...
                if(E2 > E1) ok = 1; break; end%daca eroarea creste metoda nu e convergenta
            end
        if(ok == 0) N = Nk1;%daca e convergeta functia returneaza numarul de puncte pentru care e convergenta,eltfel returneaza inf
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
	 c=difdiv(x,f(x));
	 b=c(1);
	 p=1;
	for i=2:n
		p=p*(x0-x(i-1));
		b=b+p*c(i);
    end

end

function a=difdiv(x,y)%cacularea diferentelor divizate pentru Newton
	n = length( x );
	for k = 1 : n-1
		y( k+1 : n ) = ( y( k+1 : n ) - y( k ) ) ./ ( x( k+1 : n ) - x( k ) );
    end
	a = y( : );
end

function rez = liniar( x,x0 )%se aplica metoda banala de calcul a ecuatiei unei drepte care trece prin 2 puncte 
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
function rez = splinec2(x,x0)%algoritmul din curs avea formulele diferite fata de cel din seminar asa ca am aplcat formula cum mi s-a parut mai logic
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

function rez = splineTensionat(x,x0)%am aplicat formulele din seminar/curs
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

function rez = trigonometrie( x,x0 )%Calculeaza polinomul trigonometric de interpolare...
    n=length(x);
    xi = linspace(-pi,pi,2*n+1);
    yi = f(xi);
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
function rez = f( x )%functia din enunt
    rez=exp(3*cos(x))/(2*pi*besseli(0,3));
end

