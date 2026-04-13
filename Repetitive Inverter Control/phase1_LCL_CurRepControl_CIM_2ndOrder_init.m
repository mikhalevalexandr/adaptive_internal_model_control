

% clear; clc;

% Parameters
L1 = 150e-6;     % H inverter side inductance
R1 = 0.045;      % ohm
Rd = 1;          % ohm damping resistance
Cf = 22e-6;      % F
L2 = 450e-6;     % H grid side inductance
R2 = 0.135;      % ohm
Lg=0;            % H line inductance of the grid
Rg=0;
L2g=L2+Lg;
R2g=R2+Rg;

A = [ -(R1+Rd)/L1,     Rd/L1,       -1/L1;
        Rd/L2g,     -(R2g+Rd)/L2g,      1/L2g;
        1/Cf,        -1/Cf,         0    ];


B1 = [  1/L1  0;
       -1/L2g 0;
        0     0];

B2 = [ 1/L1;
        0;
        0 ];

C = [0, 1, 0]; 
D1 = [0 -1];     
D2 = 0;            

B = [B1, B2];      
D = [D1, D2];      
P = ss(A, B, C, D);

s = tf('s');

Gfffilter=33*(0.05*s+1)/((s+300)*(0.002*s+1));%feed-forward filter

%% grid plant (weak)
Lg=0.0;
Rg=0;
L2g=L2+Lg;
R2g=R2+Rg;
Agrid = [ -(R1+Rd)/L1,     Rd/L1,       -1/L1;
        Rd/L2g,     -(R2g+Rd)/L2g,      1/L2g;
        1/Cf,        -1/Cf,         0    ];

Pgrid = ss(Agrid, B2, C, D2);
Pgridtf=ss(Pgrid);
%% IMC Design
freq=50;
w=2*pi*freq;
w1=w;
w2=3*w;
w3=5*w;
freq0=50;
w0=2*pi*freq0;
w01=w0;
w02=3*w0;
w03=5*w0;


Pu = ss(A, B2, C, D2);
Putf=tf(Pu);
A=A;
B=B2;
C=C;
D=D2;
% Regulated output selection: E*y = theta
E  = 1;    
Ecomp = 0;      

p_e=size(E,1);
p=size(C,1);
n=size(A,1);
m=size(B,2);

wc=2*pi*2550;
freq=50;
tau=1/freq;
tau_d=tau-1/wc;


delta=exp(-tau_d*s);


Wtf=wc^2/(s^2+sqrt(2)*wc*s+wc^2);
W_ss=ss(Wtf);
Aw=W_ss.A;
Bw=W_ss.B;
Cw=W_ss.C;
Dw=W_ss.D;
Wtfd=Wtf*delta;

Mtf=inv(tf(1,1)-Wtf*delta);


M0_tf=1;
M0_ss=ss(M0_tf);
assignin('base','M0_ss',M0_ss);

Am0=M0_ss.A;
Bm0=M0_ss.B;
Cm0=M0_ss.C;
Dm0=M0_ss.D;
n_m0=size(Am0,1);

M1_tf=-wc^2/(s^2+sqrt(2)*wc*s+wc^2);
M1_ss=ss(M1_tf);
assignin('base','M1_ss',M1_ss);

Am1=M1_ss.A;
Bm1=M1_ss.B;
Cm1=M1_ss.C;
Dm1=M1_ss.D;
nu = size(Bm1,2);    % number of parameters
Bsharp = (Pu.B' * Pu.B)\Pu.B';  
N = null(B2');             % columns span orthogonal complement
Bperp = N.';              % (n-m) x n
% Check nonsingularity of [Bperp; B#]:
if abs(det([Bperp; Bsharp])) < 1e-10
    error('Singularity in [Bperp; B#]. Choose a different basis for B?.');
end

EPuinvtf=inv(E*Putf); 

Pubar=Pu;
Pubartf=Putf;
EPubarinvtf=EPuinvtf;

% Design compensators
Pi1_hat_tf = EPubarinvtf-(E*Putf)\M0_tf;
Pi2_hat_tf = (Ecomp*Pubartf)*EPubarinvtf-(Ecomp*Putf)/(E*Putf)*M0_tf;
Pi1_tf = EPubarinvtf-(E*Putf)\M0_tf-(E*Putf)\M1_tf*delta;
Pi2_tf = (Ecomp*Pubartf)*EPubarinvtf-(Ecomp*Putf)/(E*Putf)*M0_tf-(Ecomp*Putf)/(E*Putf)*M1_tf*delta;


Rubar=-1.1;

Rubartf=tf(Rubar);


Tdubartf=(eye(size(Pubartf*Rubartf))-Pubartf*Rubartf)\Pubartf;
ETdubarinvtf=inv(E*Tdubartf); 
Psi_tf = [Rubartf+(1-Rubartf*Putf)/(E*Putf)*E-(E*Tdubartf)\inv(M0_tf)*E (E*Tdubartf)\inv(M0_tf)*M1_tf;
          M0_tf\E                                              -inv(M0_tf)*M1_tf];
Psi_short_tf = [Rubartf (inv(Putf)-Rubartf)*M1_tf;
          eye(p_e) M1_tf];
Psi_tf=minreal(Psi_tf,1e-5);
Psi_short_tf=minreal(Psi_short_tf,1e-5);


%check column rank
SylvCheckMx=[Bperp*A;
             -E*C];
filter1order=wc/(s+wc);
filter2order=wc^2/(s^2+sqrt(2)*wc*s+wc^2);
bodeplot(filter1order,filter2order);
xi=44;
mu=0.26;

% Generalized plant for Hinf
Agen = [A    zeros(n,size(Aw,2));
        Bw*C Aw];
Bgen = [zeros(n,1) B1 B2;
        Bw*xi   Bw*D1 Bw*D2];
Cgen = [zeros(1,n)      Cw;
        zeros(1,n)      zeros(1,size(Cw,2));
        C              zeros(p,size(Cw,2))];
Dgen = [zeros(1,1) zeros(1,2) zeros(1,1);
        zeros(1,1) zeros(1,2) mu;
        xi D1 D2];
Pgen=ss(Agen, Bgen, Cgen, Dgen);
Pgen=minreal(Pgen);


[Contr2,~,~] = hinfsyn(Pgen, 1, 1);
Contr2tf=tf(Contr2);

%% time varying freq
freq0var=50;
w0var=2*pi*freq0var;
freq1var=49.5;
w1var=2*pi*freq1var;
freq2var=50.5;
w2var=2*pi*freq2var;

%% margins

%Hinfty
L = -Mtf*Contr2tf*Putf;
w = logspace(0,5,5000);                 % tune range
Ljw = squeeze(freqresp(L,w));           % complex values
dist = abs(1 + Ljw);                    % distance to -1
[minDist, idx] = min(dist)
w_crit = w(idx)

S = feedback(1,L);
Tcl = feedback(L,1);

w = logspace(0,5,5000);
Sw = squeeze(freqresp(S,w));
Tw = squeeze(freqresp(Tcl,w));

Ms = max(abs(Sw))   % sensitivity peak
Mt = max(abs(Tw))   % complementary sensitivity peak


%CIM
Control_CIM=lft(Psi_tf,delta);



L = -Control_CIM*Putf;
w = logspace(0,5,5000);                 % tune range
Ljw = squeeze(freqresp(L,w));           % complex values
dist = abs(1 + Ljw);                    % distance to -1
[minDist, idx] = min(dist)
w_crit = w(idx)




S = feedback(1,L);
Tcl = feedback(L,1);

w = logspace(0,5,5000);
Sw = squeeze(freqresp(S,w));
Tw = squeeze(freqresp(Tcl,w));

Ms = max(abs(Sw))   % sensitivity peak
Mt = max(abs(Tw))   % complementary sensitivity peak

nyquist(-Mtf*Contr2tf*Putf); grid on
figure
nyquist(-Control_CIM*Putf); grid on