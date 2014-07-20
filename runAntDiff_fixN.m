L = 500;
n = 1000;
j = 10;
AllH = [];

% Map
m = graphMap;
m.setMap ( 13, 1, 9 );
    
for i = 1:1:n
    %% Usar o AntDiff_fixN
    % AntDiff_fixN
    adiffusion = AntDiff_fixN;
    adiffusion.setAntDiff_fixN( 100, L );
    
    % run
    adiffusion.run_AntDiff_fixN ( m );
    
    % Obter o caminho mais utilizado
    % Não há como definir uma vez que não há objetivo
    
    [HT] = janelaEntropia (j, adiffusion.trajeto);
    AllH(:,end+1) = HT;
    
    i    
end

figure ('name', 'Entropia média', 'color', 'w')
plotshaded ( AllH', 'r' );
hold off;
axis ([-1*floor(log10(L)) L 2.1 3.7]);
