%echelle
T = linspace(0,7,1400);
%période d'échantillonnage
dt = 0.005;

% ZMP en x
Zmpx = zeros(1,1400);
t = 0;
for i = 1:1400
    if(t < 2.5)
        Zmpx(i) = 0;
    else
        if(t < 3.25)
            Zmpx(i) = 0.3;
        else 
            if(t < 4.2)
                Zmpx(i) = 0.60;
            else
                if (t < 7)
                    Zmpx(i) = 0.85;
                end 
             end 
        end 
     end 
   t = t+dt; 
end

%ZMP en y
Zmpy = zeros(1,1400);
t = 0;
for i = 1:1400
    if(t < 1.80)
        Zmpy(i) = 0;
    else
        if(t < 2.7)
            Zmpy(i) = 0.1;
        else 
            if(t < 3.3)
                Zmpy(i) = -0.1;
            else 
                if(t < 4.2)
                    Zmpy(i) = 0.1;
                else
                    if (t < 5)
                        Zmpy(i) = -0.1;
                    else
                        if (t < 7)
                            Zmpy(i) = 0;
                        end
                    end
                    
                end
                
            end
        end
    end
    
   
    t = t+dt;
end

figure(1)
plot(T,Zmpx)
ylabel("x[m]")
xlabel("t(s)")
figure(2)
plot(T,Zmpy)
ylabel("y[m]")
xlabel("t(s)")

% Papier Kajita
A = [1 dt ((dt*dt)/2); 0 1 dt; 0 0 1];
B = [(dt*dt*dt)/6; (dt*dt)/2; dt];
C = [ 1 0 -0.814/9.81]; 
Q = C'*C;
R = 1e-6;
[K,P,CLC] = dlqr(A,B,Q,R);

% Partie 3 
%X
xk = zeros(3,1400);
uk = 0;
N = 160;
f = zeros(1,N);
for j = 1:1:N
    f(j) = inv((R+B'*P*B))*B'*(transpose(A-B*K)^(j-1))*C';
end

Pref = zeros(N,1);
Pk = zeros(1,1400); 
for k = 1:1:1400
    for i = 1:N
        if (k+i <= 1400)
            Pref(i,1) = Zmpx(k+i);
        else
            Pref(i,1) = Zmpx(1400);
        end
    end
    if(k ~= 1400)
    xk(:,k+1) = A*xk(:,k) + B*uk;
    end
    uk = -K*xk(:,k) + f*Pref;
    Pk(k) = C* xk(:,k); 
    
end 

figure(3)
hold on
plot(T,Pk)
plot(T,xk(1,:))

% Partie 3
% Y
xk = zeros(3,1400);
uk = 0;
N = 320;
f = zeros(1,N);
for j = 1:1:N
    f(j) = inv((R+B'*P*B))*B'*(transpose(A-B*K)^(j-1))*C';
end

Pref = zeros(N,1);
Pk = zeros(1,1400); 
for k = 1:1:1400
    for i = 1:N
        if (k+i <= 1400)
            Pref(i,1) = Zmpy(k+i);
        else
            Pref(i,1) = Zmpy(1400);
        end
    end
    if(k ~= 1400)
    xk(:,k+1) = A*xk(:,k) + B*uk;
    end
    uk = -K*xk(:,k) + f*Pref;
    Pk(k) = C* xk(:,k); 
    
end

figure(4)
hold on
plot(T,Pk)
plot(T,xk(1,:))


