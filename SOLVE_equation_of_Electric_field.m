clear
clc
close all
kB=1.38064852*10^(-23);
T=273.15;
%% ------------> Dagger
MU_dagger = 2.758627;
al_Dagger = 86.104;

%% ------------> A
MU_A = 0.0557;
al_A = 47.596667;

%% ------------> B
MU_B = 2.434532;
al_B = 32.971667;

%%
alphaD = 1.648773*(10^(-41))* al_Dagger;
alphaB = 1.648773*(10^(-41))* al_B;
alphaA = 1.648773*(10^(-41))* al_A;

muD = 3.33564*(10^(-30))* MU_dagger;
muA = 3.33564*(10^(-30))* MU_A;
muB = 3.33564*(10^(-30))* MU_B;
%%
syms E;
Delta= alphaD - (alphaA + alphaB);
EX=exp((Delta* E)/(2* kB * T));

sD = sinh((muD*E) /(kB * T));
sA = sinh((muA*E) /(kB * T));
sB = sinh((muB*E) /(kB * T));

K_rel= 1.001 ==((E.*muA.*muB)./(kB.*T .*muD)).* (EX) .*  ((sD) ./( sA .* sB ));
%%
s=vpasolve(K_rel,E)
AU=s/(5.14220674763*10^11)