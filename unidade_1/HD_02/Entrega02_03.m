clear all;clc;
% Carregando o canal real
load('Prx_Real_2021_1.mat')
vtPrxdBm = dPrx;
vtDist = dPath;
% Tamanho da janela
vtW = [2 5 10];
for dW = 1:length(vtW)
    % Transforma potência em mWatts
    vtPtrxmW = 10.^(vtPrxdBm/10);
    nSamples = length(vtPtrxmW);
    % Vetores para canal estimado
    vtDesLarga = [];
    vtDesPequeEst = [];
    %
    % Cálculo do desvanecimenro lento e rápido
    dMeiaJanela = round((vtW(dW)-1)/2);  % Meia janela
    ij = 1;
    for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
        % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
        vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
        % Desvanecimento de pequena escala [dB]
        vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
        ij = ij + 1;
    end
    %
    % Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
    indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
    %vtPrxW = ((10.^(vtPrxdBm(indexes)./10))/1000);
    vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
    desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
    envNormal = sqrt(vtPtrxmWNew')./sqrt(desLarga_Lin);
    %
    % Ajuste no tamanho dos vetores devido a filtragem
    vtDistEst = vtDist( dMeiaJanela+1 : nSamples-dMeiaJanela );
    vtPrxdBm = vtPrxdBm( dMeiaJanela+1 : nSamples-dMeiaJanela );
    %
    % Cálculo reta de perda de percurso
    vtDistLog = log10(vtDist);
    vtDistLogEst = log10(vtDistEst);
    % Cálculo do coeficientes da reta que melhor se caracteriza a perda de percurso
    dCoefReta = polyfit(vtDistLogEst,vtPrxdBm,1); 
    % Expoente de perda de percurso estimado
    dNEst = -dCoefReta(1)/10;
    % Perda de percurso estimada para os pontos de medição
    vtPathLossEst = polyval(dCoefReta,vtDistLogEst);
    % Sombreamento
    vtShadCorrEst = vtDesLarga - vtPathLossEst';
    % Calcula a variância do sombreamento estimado
    stdShad = std(vtShadCorrEst);
    meanShad = mean(vtShadCorrEst);

    fitWindow = fitmethis([envNormal],'figure','off','output','off');
    fprintf('-------------------------------------------------------------------------------------\n')
    fprintf('-------------------------------------Estimativas-------------------------------------\n')
    fprintf('\n \tJanela\tn\t\tDesvio-Padrão\tMédia\n')
    fprintf('\t%d \t%f \t%f \t%f\n',vtW(dW), dNEst, stdShad, meanShad)
    fprintf('-----------------------------------Melhores PDFs-------------------------------------\n')
    fprintf('\n \tPDFs\t\tName\t\tLogL\n')
    fprintf('\t%-10s \t%-10s \t%-10.3e\n',"Primeira",fitWindow(1).name,fitWindow(1).LL)
    fprintf('\t%-10s \t%-10s \t%-10.3e\n',"Segunda",fitWindow(2).name,fitWindow(2).LL)
    %fprintf('-------------------------------------------------------------------------------------\n')
end
