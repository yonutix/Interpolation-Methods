function  test_grafic()
x=linspace(-pi,pi,1001);%pentru functia continua
p=[1:300];%pentru functia discreta
q=f2(p);
hd=linspace(0,300,31);
hd(1)=1;
disp(hd);
h=linspace(-pi,pi,8);
g=linspace(-pi,pi,4);
j=linspace(-pi,pi,16);
dl=zeros(1,300);
dn=zeros(1,300);
dlin=zeros(1,300);
dspn=zeros(1,300);
dspt=zeros(1,300);
dt=zeros(1,300);
for i=1:300
   dl(i)=lagranged(hd,i);
   dn(i)=newtond(hd,i);
   dlin(i)=liniard(hd,i);
   dspn(i)=splinec2d(hd,i);
   dspt(i)=splineTensionatd(hd,i);
   dt(i)=trigonometrie(hd,i);
end%calculez fiecare polinom in cele 300 puncte
y=f(x);
l=zeros(1,1001);
lin=zeros(1,1001);
u=zeros(1,1001);
sn=zeros(1,1001);
st=zeros(1,1001);
nw=zeros(1,1001);
for i=1:1001
    l(i)=lagrange(h(1:8),x(i));%lagrange
    lin(i)=liniar(j(1:16),x(i));%liniar
    sn(i)=splinec2(h(1:8),x(i));%spline natural
    st(i)=splineTensionat(h(1:8),x(i));%spline tensionat
    u(i)=trigonometrie(g(1:4),x(i));%trigonometric
    nw(i)=newtond(h(1:8),x(i));%newton
end%calculez fiecre polinom in cele 1001 puncte
hold on; 
subplot(2,1,1);
plot(x,st,'c');
plot(x,lin,'y');
plot(x,y,'r');
plot(x,l,'m');
plot(x,u,'y');
plot(x,sn,'k');
plot(x,nw,'k');
hold off;
hold on;
subplot(2,1,2);
plot(p,q);
plot(p,dl,'c');
plot(p,dn,'r');
plot(p,dlin,'m');
plot(p,dspn,'k');
plot(p,dspt,'k');
plot(p,dt,'y');
hold off;
end

function rez = f( x )
    rez=exp(3*cos(x))/(2*pi*besseli(0,3));
end

function rez = f2(i)
load sunspot.dat;
y=sunspot(1:300,2);
    rez=y(i);
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

function P = lagranged( x ,x0)
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
        P=P+f2(x(k))*p1/p2;%efectuare suma finala
    end
end
function b = newtond(x,x0)
	 n=length(x);%lungimea vectorului x
	 c=difdiv(x,f2(x)');
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

function rez = liniard( x,x0 )
    n=length(x);
    a=zeros(1,n-1);
    b=a;
    for i=1:n-1
        if (x0 <= x(i+1) && x0 >= x(i))
            a=(f2(x(i+1))-f2(x(i)))/(x(i+1)-x(i));
            b=(x(i+1)*f2(x(i))-x(i)*f2(x(i+1)))/(x(i+1)-x(i));
            break;
        end
    end
    rez=a*x0+b;
end
function rez = splinec2d(x,x0)
    n = length(x);
    a = zeros(1,n);
    b = zeros(1,n-1);
    h = zeros(1,n-1);
    M = zeros(n,n);
    N = zeros(n,1);
    for i = 1:n
        a(i) = f2(x(i));
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

function rez = splineTensionatd(x,x0)%aplic algoritmul si formlee din curs/seminar...
    n = length(x);
    a = zeros(1,n);
    b = zeros(1,n-1);
    c = zeros(1,n);
    M = zeros(n,n);
    N = zeros(n,1);
    fd1 = (f2(x(2))-f2(x(1)))/(x(2)-x(1));
    fdn = (f2(x(n))-f2(x(n-1)))/(x(n)-x(n-1));
    for i = 1:n
        a(i) = f2(x(i));
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

function rez = trigonometried( x,x0 )
n=length(x);
 xi = linspace(-pi,pi,2*n+1);
            yi = f2(1:300);
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
function b = newtond(x,x0)
	 n=length(x);%lungimea vectorului x
	 c=difdivd(x,f2(x)');
	 b=c(1);
	 p=1;
	for i=2:n
		p=p*(x0-x(i-1));
		b=b+p*c(i);
    end

end

function a=difdivd(x,y)
	n = length( x );
	for k = 1 : n-1
		y( k+1 : n ) = ( y( k+1 : n ) - y( k ) ) ./ ( x( k+1 : n ) - x( k ) );
    end
	a = y( : );
end