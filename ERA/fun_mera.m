function [A,B,C,D, eig_s, residue1, sysdisc,HSVs]=fun_mera(h,n,Ts,def);
% Eigensystem Realization Algorithm (ERA)
%
% Author: Samuel da Silva - UNICAMP
% e-mail: samsilva@fem.unicamp.br
% Date: 2006/10/20
%+++++++++++++++++++++++++++++++++++
% updated from SdS's code to handle multiple signals. 
% Author: Lingling Fan - USF
% email: linglingfan@usf.edu
% Date: 2/18/2019
% 
% [A,B,C,D]=era(h,n,N,Ts,def);
% 
% Inputs:
%    h: discrete-time impulse response
%    n: order of the system
%    N: number of samples to assembly the Hankel matrix
%    Ts: sample time
%    def: if = 1: the output will be the discrete-time state-space model
%         if = 2: the output will be the continuous-time state-space model
%          
% Otputs:
%    [A,B,C,D]: state-space model
%  
% Note: For now, it works to SISO systems and it is necessary the control toolbox
% 2/18/2019: the code will handle SIMO systems. 
%


% Hankel matrix
%H0 = hankel(h(2:N+1));            % k = 0
%H1 = hankel(h(3:N+2));            % k = 1;

% input data y will be n*N dimension, where n is the # of channels, N is #
% of samples. 
[n_ch, N1] = size(h);
N = N1-2; 
L = floor(N/3); 
H0=[]; H1=[];
for k = 1: N+1-L
    H0 = [H0; h(:,k+1:k+L)];
    H1 = [H1; h(:,k+2: k+1+L)]; 
end


% Factorization of the Hankel matrix by use of SVD
[U,Sigma,V] = svd(H0);   
[U1,S1,V1] = svd(H1);   
H1_v1 = U1(:,1:n)*S1(1:n, 1:n)*V1(:,1:n)'; 
% R and S are orthonormal and Sigma is a rectangular matrix

Sr = Sigma(1:n,1:n);            
Ur = U(:,1:n);
Vr = V(:,1:n);
Co = Sr^0.5*Vr';
Ob = Ur*Sr^0.5;
% The identified system matrix are:
A = Sr^-.5*U(:,1:n)'*H1_v1*V(:,1:n)*Sr^-.5;            % dynamic matrix
B = Co(:,1);                    % input matrix
C = Ob(1:n_ch,:);                    % output matrix
%D = h(1);                       % direct-transmission matrix
%B = Sr^(-.5)*Ur'*H0(:,1);

%C = H0(1:n_ch,:)*Vr*Sr^(-.5);
D = h(:,1); % first snapshot. 
HSVs = diag(Sigma);

sysdisc = ss(A,B,C,D,Ts);       % discrete-time system

if def == 2                            
    syscont = d2c(sysdisc,'zoh');       % Conversion of discrete LTI models to continuous time
    [A,B,C,D]=ssdata(syscont);          % continuous system
end

%--------------------------------------------------------------------------
% show signal match
% show eigenvalue
z = eig(A); 
eig_s = log(eig(A))/Ts; 
%figure(888); 
figure('name', 'eig Validation'); 
scatter(real(eig_s),imag(eig_s)/2/pi,'LineWidth',2); 
ylabel('Hz')
grid on; 

%show match by reconstructing signals. 
if(n_ch<= 5)
    row_plot = n_ch;
    col_plot = 1; 
else
    if(mod(sqrt(n_ch),1)>0)
        row_plot= floor(sqrt(n_ch))+1; 
    else
        row_plot=sqrt(n_ch);
    end
    if (mod(n_ch/row_plot,1)>0)
        col_plot = floor(n_ch/row_plot) + 1; 
    else 
        col_plot = n_ch/row_plot; 
    end
    % e.g., 6 signals: 3*2
    % e.g.; 7 signals: 3*3
end


ya =h.'; 
%% signal reconstruction 
for i1=1:N+2;
    for j1=1:n;
        Z(i1,j1)=z(j1)^(i1-1);
    end 
end
for i=1:size(ya,2) % for three signal, reconstruct
residue1(:,i) = pinv(Z)*ya(:,i);
y_hat(:,i)=Z*residue1(:,i);
end

%figure(998); 
figure('name', 'ERA Validation'); 
t1 = 0:Ts:(N+1)*Ts;
for i=1:n_ch
    subplot(row_plot, col_plot, i);      
    plot(t1, ya(:,i),'Linewidth',2); hold on;
    plot(t1, real(y_hat(:,i)),'r','Linewidth',1);  
    grid on
    legend('Original', 'ERA');
end


