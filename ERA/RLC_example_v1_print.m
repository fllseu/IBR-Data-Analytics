%% step 1: Generate time-domain data -- step response
clear; clc;
s= tf('s'); h = 0.001; t = 0:h:0.1;
vs = ones(length(t),1);                    % step excitation
R = 0.1; L = 0.5/377; C = 1/0.2/377;
sys = 1/(R*C*s+L*C*s^2+ 1);                % from vs to capacitor voltage 
vc = lsim(sys, vs, t); [n_vc, m_vc] = size(vc);
vc1 = vc + (rand(n_vc, m_vc) -0.5)*0.1/max(vc); % add noise

%% Step 2: run ERA to get necessary information: eigenvalues and residues
m = 3;                                    % order of the system  
y = vc1; [N, n_ch]=size(y); 
[A1,B1,C1,D1, eig_s,residue1, sysdisc, HSVs]=fun_mera(y',m ,h, 1);

%% Step 3: Assemble a transfer function 
res1 = residue1(:,1); 
s = tf('s'); Freq_s = 0; 
for i=1:length(eig_s)
    Freq_s = Freq_s + res1(i)/(s-eig_s(i));
end
G = Freq_s*s;     % consider the step excitation as the input





