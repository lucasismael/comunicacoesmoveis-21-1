clear all;clc;
% Carregando o canal real
load('Prx_Real_2021_1.mat')
vtPrxdBm = dPrx;
vtDist = dPath;
% Tamanho da janela
dW = 5;
% Transforma pot�ncia em mWatts
vtPtrxmW = 10.^(vtPrxdBm/10);
nSamples = length(vtPtrxmW);
% Vetores para canal estimado
vtDesLarga = [];
vtDesPequeEst = [];
%
% C�lculo do desvanecimenro lento e r�pido
dMeiaJanela = round((dW-1)/2);  % Meia janela
ij = 1;
for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
    % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
    vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
    % Desvanecimento de pequena escala [dB]
    vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
    ij = ij + 1;
end
%
% C�lculo da envolt�ria normalizada (para efeitos de c�lculo do fading)
indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
%vtPrxW = ((10.^(vtPrxdBm(indexes)./10))/1000);
vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
envNormal = sqrt(vtPtrxmWNew)./sqrt(desLarga_Lin);
%
% Ajuste no tamanho dos vetores devido a filtragem
vtDistEst = vtDist( dMeiaJanela+1 : nSamples-dMeiaJanela );
vtPrxdBm = vtPrxdBm( dMeiaJanela+1 : nSamples-dMeiaJanela );
%
% C�lculo reta de perda de percurso
vtDistLog = log10(vtDist);
vtDistLogEst = log10(vtDistEst);
% C�lculo do coeficientes da reta que melhor se caracteriza a perda de percurso
dCoefReta = polyfit(vtDistLogEst,vtPrxdBm,1); 
% Expoente de perda de percurso estimado
dNEst = -dCoefReta(1)/10;
%disp(['Estima��o dos par�metros de larga escala (W = ' num2str(dW) '):'])
%disp(['   Expoente de perda de percurso estimado n = ' num2str(dNEst)]);
% Perda de percurso estimada para os pontos de medi��o
vtPathLossEst = polyval(dCoefReta,vtDistLogEst);
% Sombreamento
vtShadCorrEst = vtDesLarga - vtPathLossEst';
% Calcula a vari�ncia do sombreamento estimado
stdShad = std(vtShadCorrEst);
meanShad = mean(vtShadCorrEst);
%disp(['   Desvio padr�o do sombreamento estimado = ' num2str(stdShad)]);
%disp(['   M�dia do sombreamento estimado = ' num2str(meanShad)]);
vtPathLossEst = - vtPathLossEst;
%Full Signal
plot(vtDistLogEst,vtPrxdBm, '-','linewidth',2);hold on;
%Pathloss + Shadowing
vtwoFading = vtPrxdBm-vtDesPequeEst';
plot(vtDistLogEst,vtwoFading, '-','linewidth',2);
%Pathloss only
vtwoFadingShad = vtPrxdBm-vtDesPequeEst'-vtShadCorrEst';
plot(vtDistLogEst,vtwoFadingShad,'-','linewidth',2);
%Plot configs
title("Canal Real: Influ�ncia dos efeitos do canal (W =" + dW +" )");
xlabel('log_{10}(d)');
ylabel('Pot�ncia [dBm]');
legend('Prx canal completo', 'Prx (perda de percurso + sombreamento)','Prx (somente perda de percurso)');
axis([0.5 1.3 -100 -30]);
grid;