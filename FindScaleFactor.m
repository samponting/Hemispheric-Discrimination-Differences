CF = csvread('linss10e_1.csv');

Scone = CF(:,4);

Lw = 0.692839;
Mw = 0.349676;
Lum = (Lw.*CF(:,2)) + (Mw.*CF(:,3));

SonLplusM = Scone ./ Lum;
max(SonLplusM);
Sscaled = SonLplusM / max(SonLplusM);


Sw = 1/max(SonLplusM);


%%
% Equal energy white
EEW = ones([441 1]);

% Find cone responses 
L = EEW'*CF(:,2);
M = EEW'*CF(:,3);
S = EEW'*CF(:,4);

% Find MB axes
r = (Lw*L)./((Lw*L)+(Mw*M));
b = (Sw*S)./((Lw*L)+(Mw*M));

X = (S)./((Lw*L)+(Mw*M));

Sw = 1/X;
