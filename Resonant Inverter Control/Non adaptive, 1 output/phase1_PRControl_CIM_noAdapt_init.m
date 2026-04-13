
clear; clc; close all;

% Parameters
L1 = 150e-6;     % H inverter side inductance
R1 = 0.045;      % ohm
Rd = 1;          % ohm damping resistance
Cf = 22e-6;      % F
L2 = 450e-6;     % H grid side inductance
R2 = 0.135;      % ohm
Lg=0.00;         % H line inductance of the grid
Rg=0;
L2g=L2+Lg;
R2g=R2+Rg;


A = [ -(R1+Rd)/L1,     Rd/L1,       -1/L1;
        Rd/L2g,     -(R2g+Rd)/L2g,      1/L2g;
        1/Cf,        -1/Cf,         0    ];

B1 = [  1/L1;
       -1/L2g;
        0];

B2 = [ 1/L1;
        0;
        0 ];

C = [0, 1, 0]; 
D1 = 0;      
D2 = 0;          

B = [B1, B2];      
D = [D1, D2];      
P = ss(A, B, C, D);

s = tf('s');



Gfffilter=33*(0.05*s+1)/((s+300)*(0.002*s+1));%feed-forward filter


%% CIM
Pu = ss(A, B2, C, D2);
A=A;
B=B2;
C=C;
D=D2;

Putf = tf(Pu);

E  = 1;    
Ecomp = 0;    

p_e=size(E,1);
p=size(C,1);
n=size(A,1);
m=size(B,2);
freq=50;
w=2*pi*freq;

am=2*pi;
kef=80;
phi=(s+am*kef)^6;
w1=w;
w2=3*w;
w3=5*w;

Mtf = phi / ( (s^2 + w1^2) * (s^2 + w2^2) * (s^2 + w3^2));

[Am,Bm,Cm,Dm] = ssdata(ss(Mtf));
M  = ss(Am,Bm,Cm,Dm);

nu = size(Bm,2);    % number of parameters
n_m=size(Am,1);

Bsharp = (Pu.B' * Pu.B)\Pu.B';     % m x n
N = null(B');             % columns span orthogonal complement
Bperp = N.';              % (n-m) x n
% Check nonsingularity of [Bperp; B#]:
if abs(det([Bperp; Bsharp])) < 1e-10
    error('Singularity in [Bperp; B#]. Choose a different basis for B?.');
end


% Solve generalized Sylvester for X 

Asylv=[Bperp;zeros(p_e,n)];
Bsylv=Am-Bm*Cm;

Csylv=-[Bperp*A;-E*C];
Dsylv=eye(size(Bsylv));
Esylv=[zeros(n-m,n_m);Cm];
X = generalized_sylvester(Asylv, Bsylv, Csylv, Dsylv, Esylv);
C0 = Bsharp * X * (Am - Bm*Cm) - Bsharp * A * X;

Pai1=ss(Am-Bm*Cm, Bm, C0, zeros(m,p_e));
Pai2=ss(Am-Bm*Cm, Bm, Ecomp*C*X+Ecomp*D*C0, zeros(m,p_e));

A_np = A + X*Bm*E*C; 
B_np = B;
C_np = C;
D_np = D;
Pubar = ss(A_np, B_np, C_np, D_np);
Pubartf=tf(Pubar);

Rubar=-5;

Rubartf=tf(Rubar);
eig(A_np+B*Rubar*C_np)

Aim = Am;
Bim = Bm*E;
C1  = C*X+Ecomp'*Ecomp*D*C0;
C2  = -C0;  
Ctot=[C1;C2];
D1  = eye(size(C1,1),size(Bim,2));                     
D2=zeros(size(C2,1),size(Bim,2));     
Dtot=[D1;D2];                  

Gamma = ss(Aim, Bim, Ctot, Dtot);        

R= [Rubar eye(m)]*Gamma;
opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';
opts.FreqUnits = 'Hz';
bodeplot(R*Putf,opts);
grid on
