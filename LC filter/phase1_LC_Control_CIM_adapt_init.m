
clear; clc; close all;

% Parameters (Table 3.1)
L1 = 150e-6;     % H inverter side inductance
R1 = 0.045;%0.1;      % ohm
Rd = 1;        % ohm damping resistance
Cf = 22e-6;      % F

% % State matrix A (3x3)
A = [ -(R1+Rd)/L1,       -1/L1;
        1/Cf,            0    ];

% B1 for w=[v_ref; v_g] (3x2)
B1 = [  1/L1;
        0];

% B2 for u (3x1)
B2 = [ 1/L1;
        0 ];

% Output y = e (1x3), with C1, D1, D2 as in the screenshot
C = [1, 0;
     0, 1];     % C1
D1 = [0;
      0];        % maps w to y
D2 = [0;
      0];             % maps u to y



s = tf('s');

%% IMC
Pu = ss(A, B2, C, D2);
A=A;
B=B2;
C=C;
D=D2;

Putf = tf(Pu);
p_des = [-50000, -10000];

K = -place(A,B2,p_des);
Acl = A + B2*K;
Bcl = B2;          % for exogenous channels we’ll add in the 2-DOF section
Ccl = C;
Dcl = D2;
eigAcl = eig(Acl)


E  = 1;         % (p_e = 1), measured y has p = 2
Ecomp = 0;      % E? (complement) – any matrix s.t. [E;Ecomp] is orthonormal-ish row-wise

p_e=size(E,1);
p=size(C,1);
n=size(A,1);
m=size(B,2);
% Internal model M(s) 1 freq
freq=1;
w=2*pi*freq;

% Internal model M(s) 2 freq
am=2*pi*8;

w1=w;
w2=3*w;

freq0=1;
w0=2*pi*freq0;
w01=w0;
w02=3*w0;
freq0var=freq0*3;
freq1var=1.2*freq0*3;

phi=(s+am)^2;

Theta=w2^2;
Theta0=w02^2;

M0_tf=s^2/phi;
M1_tf=1/phi;
Mtf = phi / ( (s^2 + w2^2));



% controllable companion
sys0_ccf = compreal(M0_tf,'c');
sys1_ccf = compreal(M1_tf,'c');

Am0_can = sys0_ccf.A;
Bm0_can = sys0_ccf.B;
Cm0_can = sys0_ccf.C;
Dm0_can = sys0_ccf.D;
M0_can=ss(Am0_can,Bm0_can,Cm0_can,Dm0_can);
M0_can_tf=tf(M0_can);
Am1_can = sys1_ccf.A;
Bm1_can = sys1_ccf.B;
Cm1_can = sys1_ccf.C;
Dm1_can = sys1_ccf.D;
M1_can=ss(Am1_can,Bm1_can,Cm1_can,Dm1_can);
M1_can_tf=tf(M1_can);


n_m0=size(Am0_can,1);
nu=size(Cm1_can,1);


Bsharp = (Pu.B' * Pu.B)\Pu.B';     % m x n Bsharp=[0.5 0]  (Bsharp * B = 1) 

A_psi=Am0_can-Bm0_can*Cm0_can;
C_psi=[-Cm0_can;
     Cm1_can];
B_psi=[Bm0_can Bm0_can*(K+Bsharp*A)-(Am0_can-Bm0_can*Cm0_can)*Bm0_can*Bsharp -Bm0_can];
D_psi=[eye(1) K+Cm0_can*Bm0_can*Bsharp -eye(1);
       zeros(nu,1) -Cm1_can*Bm0_can*Bsharp zeros(nu,1)];


SFPsi=ss(A_psi, B_psi, C_psi, D_psi);
SFPsi_min=minreal(SFPsi,1e-3);
Pi=(1-inv(Mtf))*(Bsharp*s-Bsharp*A);
v2u=Mtf/(1-Mtf*(K-Pi)*Putf);
dcg=dcgain(v2u);
v2umin=minreal(v2u,1e-6);


SFPsitf=tf(SFPsi_min);
SFPsitf_min=minreal(SFPsitf,1e-8);

tol_minreal = 1e-6;     % for cancellations

SFPsi_tf = minreal(tf(SFPsi), tol_minreal);

[num,den] = tfdata(SFPsi_tf);   % for SISO; for MIMO use 'cell'

[rows, cols] = size(num);
for i = 1:rows
    for j = 1:cols
        % Round to 2 decimal places (or use just round(x) for integers)
        num{i,j} = round(num{i,j}, 3);
        den{i,j} = round(den{i,j}, 3);
    end
end
SFPsi_tf = tf(num, den);

Q = 30000000*eye(n);  P = lyap(Acl',Q);