clear
clc
close all
%%
kB=1.38064852*10^(-23);
T=273.15;
Order=5;
E = -10^Order:(10^(Order-3)):10^Order-1;

%% ------------> Dagger
MU_dagger = 3;
al_Dagger = 50;

%% ------------> A
MU_A = 5;
al_A = 50;

%% ------------> B
MU_B = 3;
al_B = 50;

%%
alphaD = 1.648773*(10^(-41))* al_Dagger;
alphaB = 1.648773*(10^(-41))* al_B;
alphaA = 1.648773*(10^(-41))* al_A;

muD = 3.33564*(10^(-30))* MU_dagger;
muA = 3.33564*(10^(-30))* MU_A;
muB = 3.33564*(10^(-30))* MU_B;
%%
Delta= alphaD - (alphaA + alphaB);
EX=exp((Delta* E)/(2* kB * T));

sD = sinh((muD*E) /(kB * T));
sA = sinh((muA*E) /(kB * T));
sB = sinh((muB*E) /(kB * T));

%%
K_rel=((E.*muA.*muB)./(kB.*T .*muD)).* (EX) .*  ((sD) ./( sA .* sB ))
plot(E,K_rel);
