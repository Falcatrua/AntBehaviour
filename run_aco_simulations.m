L = 500;
n = 1000;
count = 0;
count2 = 0;
j = 10;
FoundPaths = sparse([]);
Fcount = [];
NonFoundPaths = sparse([]);
NFcount = [];

% Map
m = graphMap;
m.setMap ( 13, 1, 9 );
    
% Comparar com o Dijkstra
map_clone = sparse ( m.adj );
[d, pathD ]= graphshortestpath( map_clone, m.source, m.feeder,...
    'Directed', 'false', 'Method', 'Dijkstra' );

for i = 1:1:n
    %% Usar o ACO
    % ACO
    a = ACO;
    a.setAco ( 1, L );
    
    % Map
    m = graphMap;
    m.setMap ( 13, 1, 9 );

    
    % run
    a.run_ACO ( m );
    
    % Obter o caminho mais utilizado
    [path, msg] = caminhoMaisUsado ( m );
   
    % Obter curva de crescimento numérico
    for ii = 1:1:size(a.trajeto,2)
        ac(ii,1) = length(find(a.trajeto(:,ii)));
    end
    
    non = 1;
    [HT] = janelaEntropia (j, a.trajeto);
    if length(path) == length(pathD)
        %if path == pathD
        count = count+1;
        %display ('Os caminhos são iguais');
        %FoundPaths(:,count) = HT;
        for ii = 1:1:length(HT)
            FoundPaths(ii,count) = HT(ii);            
            %i
        end
        Fcount(count,:) = ac;
        non = 0;
        %end
    end
    if non == 1
        count2 = count2 +1;
        NonFoundPaths(:,count2) = HT;
        NFcount(count2,:) = ac;
    end
    
    i    
end

figure ('name', 'Entropia média', 'color', 'w')
if ( size(FoundPaths,2) ~= 0 )
    plotshaded ( FoundPaths', 'b' );
end
hold on;
if ( size(NonFoundPaths,2) ~= 0 )
    plotshaded ( NonFoundPaths', 'r' );
end
l = legend ('caminho certo', ' ', 'caminho errado', ' ');
set (l, 'box' ,'off');
axis ([-1*floor(log10(L)) L 2.1 3.7]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');

if ( size(FoundPaths,2) ~= 0 )
    plotshaded ( Fcount, 'c' );
end
hold on;
if ( size(NonFoundPaths,2) ~= 0 )
    plotshaded ( NFcount, 'm' );       
end




