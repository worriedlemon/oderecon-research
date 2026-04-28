load real_data\coffee.mat

% Find a cooling cup equation
t = t(1:79)' / 60;
T = coffee_temp(1:79)';
tair = t_amb;

%Find ordinary integration optimization problem solution
K0 = -0.02;
[Kopt, err1] = fminunc(@(K) CoffeeObj(@SimCup, K, t, T, tair),K0);
[tspan1, y1] = SimCup(Kopt,T(1), t(end), tair);

%Find ordinary 2 order integration optimization problem solution
K1 = -50;
K2 = -40;
K3 = 0;
K4 = 1;
[xopt2, err2] = fminunc(@(x) CoffeeObj(@SimCup2, x, t, T, tair), [K1 K2 K3 K4]);
[tspan2, y2] = SimCup2(xopt2, T(1), t(end), tair);

%Find fractional-order integration optimization problem solution
K0 = -0.2;
a0 = 1.0;
[xopt3, err3] = fminunc(@(x) CoffeeObj(@SimCupFraq, x, t, T, tair), [K0 a0]);
[tspan3, y3] = SimCupFraq(xopt3, T(1), t(end), tair);

figure(1);
plot(t, T, 'o');
hold on;
plot(tspan1, y1, '-');
plot(tspan2, y2 ,'-');
plot(tspan3, y3, '-', LineWidth=2);
xlabel('$t$, min', Interpreter='latex');
ylabel('$T$, $^\circ$C', Interpreter='latex');
legend('Experimental data', ...
    sprintf('First order integrator, $K$ = %.2f', Kopt), ...
    sprintf('Second order equation, $K_1$ = %0.2f, $K_2$ = %.2f, $K_3$ = %.2f, $K_4$ = %.2f', xopt2(1), xopt2(2), xopt2(3), xopt2(4)), ...
    sprintf('Fractional-order integrator, $K$ = %.2f, $\\alpha$ = %.2f', xopt3(1), xopt3(2)), ...
    Interpreter='latex');
 

disp('Error, 1st order integrator:');
disp(err1);
disp('Error, second order integrator:');
disp(err2);
disp('Error, fractional order integrator:');
disp(err3);


%% Objectvie function
function err = CoffeeObj(SimCup, K, t, T, tair)
%Objective function for optimization: quadratic error
    [tspan, y] = SimCup(K, T(1), t(end), tair);
    N = length(tspan);
    M = length(t);
    
    err = 0; nf = 0;
    %compare
    for i = 1:N
        for j = 1:M
            if(t(j) == tspan(i))
                err = err + (y(i) - T(j))^2;
                nf = nf + 1;
            end
        end
    end
    err = sqrt(err)/ nf;
end

function [tspan, y] = SimCup(K, T0, tmax, tair)
    %Simulate the cup behaviour
    h = 0.05;
    
    tau = T0;
    tspan = 0:h:tmax;
    N = length(tspan);
    y = zeros(1,N);
    y(1) = tau;
    
    for i = 2:N
        taufun = @(x)K*(x - tair);
        tau = RK4Step(taufun, h, tau);
        y(i) = tau;
    end
end

function [tspan, y] = SimCupFraq(x, T0, tmax, tair)
    %Simulate the cup behaviour with fractional order integration 
    h = 0.05;
    [tspan, y] = fractode(@(t,y)(x(1)*(y - tair)), 0:h:tmax, T0, x(2));
end

function [tspan, y] = SimCup2(x, T0, tmax, tair)
    %Simulate the cup behaviour
    h = 0.05;
    
    tau = T0;
    tspan = 0:h:tmax;
    N = length(tspan);
    y = zeros(1,N);
    
    
    K1 = x(1);
    K2 = x(2);
    K3 = x(3);
    K4 = x(4);
    r = 0;
    
    y(1) = tau;
    
    for i = 2:N
        taufun = @(x)(K3*(x) + K4*r);
        rfun = @(x)(K1*(tau - tair) + K2*x);
        tau = RK4Step(taufun, h, tau);
        r = RK4Step(rfun, h, r);
        
        y(i) = tau;
    end
end

function x = RK4Step(fun, h, xprev)
    b1 = 1/6;
    b2 = 1/3;
    b3 = 1/3;
    b4 = 1/6;

    k1 = fun(xprev);
    x2 = xprev + 0.5*h*k1;
    k2 = fun(x2);
    x3 = xprev + 0.5*h*k2;
    k3 = fun(x3);
    x4 = xprev + h*k3;
    k4 = fun(x4);
    x = xprev + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4);
end

