clear 
close all
clc
close all
disp('programmer:Seid Saeed Mirbagheri (400126116)')
bb=1; % Select accuracy
if bb==1
N=[33 41 49 57 65 ];
else
N=[129 137 145 153 161 ];
end

n_N=length(N);
Fac=zeros(4,n_N);
U_half=zeros(4,n_N);
H=zeros(1,n_N);

for ii=1:n_N    
%% Defining the problem domain

n_points = N(ii); % no_of_points
dom_length = 1;
h = dom_length/(n_points-1);
x = 0:h:dom_length; %X domain span
y = 0:h:dom_length; %Y domain span
Re = 1; %Reynolds number
nu = 1/Re;
Error1=zeros(n_points-1,n_points-1,5);
% Under-relaxation factors
alpha = 0.8;
alpha_p = 0.8;
U_inf= 1;
%% Initializing the variables
%Final collocated variables
u_final=zeros(n_points,n_points);
v_final=zeros(n_points,n_points);
p_final=zeros(n_points,n_points);
u_final(1,:) = U_inf;

%Staggered variables
u=zeros(n_points+1,n_points);
u_star=zeros(n_points+1,n_points);
d_e=zeros(n_points+1,n_points);
v=zeros(n_points,n_points+1);
v_star=zeros(n_points,n_points+1);
d_n=zeros(n_points,n_points+1);
p=zeros(n_points+1,n_points+1);
p_star=zeros(n_points+1,n_points+1);

p(n_points+1,n_points+1)=1;
p_star(n_points+1,n_points+1)=1;
pc=zeros(n_points+1,n_points+1);
b=zeros(n_points+1,n_points+1);
u(1,:)=2*U_inf;

u_new=zeros(n_points+1,n_points);
v_new=zeros(n_points,n_points+1);
p_new=zeros(n_points+1,n_points+1);
p_new(n_points+1,n_points+1)=1;
u_new(1,:)=2*U_inf;
tab_data=zeros(17,2);
%% Solving the governing equations
error = 1;
iterations = 0;
error_req = 1e-6; %final required error residual
figure(1); %for error monitoring

while error > error_req
    % x-momentum eq. - Interior
    for i = 2:n_points
        for j = 2:n_points - 1
            u_E = 0.5*(u(i,j) + u(i,j+1));
            u_W = 0.5*(u(i,j) + u(i,j-1));
            v_N = 0.5*(v(i-1,j) + v(i-1,j+1));
            v_S = 0.5*(v(i,j) + v(i,j+1));
            
            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;
            
            a_e = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;
            
            A_e = -h;
            d_e(i,j) = A_e/a_e;
            
            u_star(i,j) = (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j))/a_e + d_e(i,j)*(p(i,j+1) - p(i,j));
        end
    end
    
    % x-momentum eq. - Boundary
    u_star(1,:) = 2.*U_inf - u_star(2,:);
    u_star(n_points + 1,:) = -u_star(n_points,:);
    u_star(2:n_points,1) = 0;
    u_star(2:n_points,n_points) = 0;
    
    % y-momentum eq. - Interior
    for i = 2:n_points - 1
        for j = 2:n_points
            u_E = 0.5*(u(i,j) + u(i+1,j));
            u_W = 0.5*(u(i,j-1) + u(i+1,j-1));
            v_N = 0.5*(v(i-1,j) + v(i,j));
            v_S = 0.5*(v(i,j) + v(i+1,j));
            
            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;
            
            a_n = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;
            
            A_n = -h;
            d_n(i,j) = A_n/a_n;
            
            v_star(i,j) = (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j))/a_n + d_n(i,j)*(p(i,j) - p(i+1,j));
        end
    end
    
    % y-momentum eq. - Boundary
    v_star(:,1) = -v_star(:,2);
    v_star(:,n_points + 1) = -v_star(:,n_points);
    v_star(1,2:n_points) = 0;
    v_star(n_points,2:n_points) = 0;
    
    % Zeroing the corrections to begin with
    pc(1:n_points+1,1:n_points+1)=0;
    
    % Continuity equation a.k.a. pressure correction - Interior
    for i = 2:n_points
        for j = 2:n_points
            a_E = -d_e(i,j)*h;
            a_W = -d_e(i,j-1)*h;
            a_N = -d_n(i-1,j)*h;
            a_S = -d_n(i,j)*h;
            a_P = a_E + a_W + a_N + a_S;
            b(i,j) = -(u_star(i,j) - u_star(i,j-1))*h + (v_star(i,j) - v_star(i-1,j))*h;
            
            pc(i,j) = (a_E*pc(i,j+1) + a_W*pc(i,j-1) + a_N*pc(i-1,j) + a_S*pc(i+1,j) + b(i,j))/a_P;
        end
    end
    
    % Correcting the pressure field
    for i = 2:n_points
        for j = 2:n_points
            p_new(i,j) = p(i,j) + alpha_p*pc(i,j);
        end
    end
    
    % Continuity eq. - Boundary
    p_new(1,:) = p_new(2,:);
    p_new(n_points + 1,:) = p_new(n_points,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,n_points + 1) = p_new(:,n_points);
    
    % Correcting the velocities
    for i = 2:n_points
        for j = 2:n_points - 1
            u_new(i,j) = u_star(i,j) + alpha*d_e(i,j)*(pc(i,j+1) - pc(i,j));
        end
    end
    
    % x-momentum eq. - Boundary
    u_new(1,:) = 2*U_inf - u_new(2,:);
    u_new(n_points + 1,:) = -u_new(n_points,:);
    u_new(2:n_points,1) = 0;
    u_new(2:n_points,n_points) = 0;
    
    for i = 2:n_points - 1
        for j = 2:n_points
            v_new(i,j) = v_star(i,j) + alpha*d_n(i,j)*(pc(i,j) - pc(i+1,j));
        end
    end
    
    % y-momentum eq. - Boundary
    v_new(:,1) = -v_new(:,2);
    v_new(:,n_points + 1) = -v_new(:,n_points);
    v_new(1,2:n_points) = 0;
    v_new(n_points,2:n_points) = 0;
            
    
    % Continuity residual as error measure
    error = 0;
    for i = 2:n_points
        for j = 2:n_points
            error = error + abs(b(i,j));
        end
    end
    
    % Finishing the iteration
    u = u_new;
    v = v_new;
    p = p_new;
    iterations = iterations + 1;
    Error1(iterations,ii)=error;
end

% After the converged solution, we map the staggered variables to
% collocated variables
for i = 1:n_points
    for j = 1:n_points
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end


%% Contour and vector visuals.
x_dom = ((1:n_points)-1).*h;
y_dom = 1-((1:n_points)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);

figure(1);
%plot U velocity_X
subplot(2,3,ii)
contourf(X,Y,u_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf(' U velocity_X(m/s) for N= %d ' ,round(N(ii))))

figure(2);
%plot U velocity_X
subplot(2,3,ii)
contourf(X,Y,v_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf(' V velocity_Y(m/s) for N= %d ' ,round(N(ii))))

figure(3);
%Pressure(pa)
subplot(2,3,ii)
contourf(X,Y,p_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf('Pressure(pa) for N= %d ' ,round(N(ii))))

figure(4);
%plot Total velocity vector
subplot(2,3,ii)
hold on
grid on
quiver(X, Y, u_final, v_final, 5, 'r')
title(sprintf('Total velocity vector for N= %d ',round(N(ii))))

figure(5);
% Error 
ite=1:iterations;
figure(10);
subplot(3,3,ii)
plot(ite,Error1(:,ii))
hold on
xlabel('Iterations')
ylabel('Residual Error')
title(sprintf('Residual Error(Iterations) for N= %d ' ,round(N(ii))))


%% Centreline u validation - Comparison with benchmark for Re=100 and U_inf==2
if Re==100 && U_inf==1
figure(11);
subplot(3,3,ii)
plot(u_final(:,(n_points+1)/2),1-y, 'LineWidth', 1)
tab_data = xlsread('C:\Users\King\Desktop\data for run matlab\plot_u_y_Ghia100.csv','A2:B18');
y_ghia = tab_data(:,1);
u_ghia = tab_data(:,2);

figure(11); hold on
plot(u_ghia, y_ghia, 'o', 'LineWidth', 1)
title(sprintf('Residual Error(Iterations) for N= %d ' ,round(N(ii))))
xlabel('u')
ylabel('y')
legend('Numerical', 'Benchmark', 'location', 'southeast')
end

    H(ii)=h;
    Fac(:,ii)=[find(y==0.25);find(y==0.5);find(y==0.625);find(y==0.75)];
for i=1:4
    U_half(i,ii)=u_final((n_points+1)/2,Fac(i,ii));
end

end 


Er=Error(U_half,n_N);
R_e_delta_y=Error_Slope( Er,H,n_N );

%%%%%Draw the relative error value in terms of H
COL=['s','*','d','o'];
for i=1:4
figure(12)
subplot(1,2,1)
loglog(H(1:end-1),abs(Er(i,:)),sprintf(COL(i)),'linewidth',2);
title('Solution convergence )');
ylabel('successive error');
xlabel('h');
hold on
grid on

%%%%Draw the error slope in terms of h
subplot(1,2,2)
loglog(H(1:end-2),abs(R_e_delta_y(i,:)),'linewidth',2);
title('Error slope ');
ylabel('Slope');
xlabel('h');
hold on
grid on
end
legend('y=0.25','y=0.5','y=0.75','y=0.875')
