classdef graphMap < handle
    properties
        nodesN = []; % número de nós
        adj = []; % matriz de adjacências
        edgesPh = [];%mapa com valores de feromônios;
        edgesLen = [];%mapa com valores de distâncias;
        nodesC = []; %contagem do número de visitas em cada nó
        feeder = []; % nó com alimento
        source = []; % nó do ninho
        gr;
    end

    methods
        % Settings
        function setMap ( graphMap, n_nodes, source, feeder )
            graphMap.nodesN = n_nodes;
            graphMap.edgesPh = ones(graphMap.nodesN,graphMap.nodesN);%mapa com valores de feromônios;
            graphMap.edgesLen = ones(graphMap.nodesN,graphMap.nodesN)*2;%mapa com valores de distâncias;
            graphMap.adj = graphMap.mazeAdjacencyMatrix;
            graphMap.edgesPh(graphMap.adj==0)=0;
            graphMap.edgesLen(graphMap.adj==0)=0;
            graphMap.nodesC = zeros(graphMap.nodesN,1); %contar o número de visitas em cada nó
            graphMap.feeder = feeder;
            graphMap.source = source;
            
            graphMap.adj = graphMap.symm ( graphMap.adj );
        end
        % graph generation
        function [maze_vertex] = mazeAdjacencyMatrix ( graphMap )
            maze_vertex = [0 1 1 1 1 1 1 0 0 0 0 0 0;...
                1 0 1 0 0 0 1 1 0 0 0 0 0;...
                1 1 0 1 0 0 0 0 1 0 0 0 0;...
                1 0 1 0 1 0 0 0 0 1 0 0 0;...
                1 0 0 1 0 1 0 0 0 0 1 0 0;...
                1 0 0 0 1 0 1 0 0 0 0 1 0;...
                1 1 0 0 0 1 0 0 0 0 0 0 1;...
                0 1 0 0 0 0 0 0 1 0 0 0 1;...
                0 0 1 0 0 0 0 1 0 1 0 0 0;...
                0 0 0 1 0 0 0 0 1 0 1 0 0;...
                0 0 0 0 1 0 0 0 0 1 0 1 0;...
                0 0 0 0 0 1 0 0 0 0 1 0 1;...
                0 0 0 0 0 0 1 1 0 0 0 1 0];
        end

        % Utils
        function [m_out] = symm ( graphMap,  m )
            m_out = (tril(m,-1)' + tril(m));
        end

        % Plottings
        function plot ( graphMap )
            nodes_ids = cellfun (@num2str,num2cell([1:size(graphMap.adj,1)]),'un',0);
            graphMap.gr = biograph(graphMap.adj, nodes_ids, 'LayoutType', 'radial');

            graphMap.edgesPh = (graphMap.edgesPh-min(min(graphMap.edgesPh(graphMap.edgesPh>0))))./...
                (max(max(graphMap.edgesPh(graphMap.edgesPh>0)))-min(min(graphMap.edgesPh(graphMap.edgesPh>0))));
            graphMap.edgesPh(graphMap.edgesPh<0)=0;
            for i = 1:1:numel(graphMap.gr.edges)
                nodes = regexp(graphMap.gr.edges(i).ID,' -> ', 'split');
                nn = cellfun(@str2num,nodes);
               graphMap.gr.edges(i).Weight = graphMap.edgesPh(nn(1),nn(2));
                set (graphMap.gr.edges(i), 'LineColor', [0.2 0.2 0.2]);
                if (graphMap.edgesPh(nn(1),nn(2))>0)
                    set (graphMap.gr.edges(i), 'LineColor', 1-([graphMap.edgesPh(nn(1),nn(2))]*[0 1 1]));
                    set (graphMap.gr.edges(i), 'LineWidth', 1.5);
                else
                    set (graphMap.gr.edges(i), 'LineColor', [0.85 0.85 0.85]);
                    set (graphMap.gr.edges(i), 'LineWidth', 0.3);
                end
            end
            graphMap.nodesC = (graphMap.nodesC-min(graphMap.nodesC))./(max(graphMap.nodesC)-min(graphMap.nodesC));
            for i = 1:1:numel(graphMap.gr.nodes)
                node = str2double(graphMap.gr.nodes(i).ID);
                set (graphMap.gr.nodes(node), 'Color', [1 1 1],...
                    'LineWidth', 5*graphMap.nodesC(i)+0.1, 'Shape', 'circle',...
                    'Size', ((graphMap.nodesC(i))*10+8)*[1 1]);
                set (graphMap.gr.nodes(node), 'LineColor', 1-(graphMap.nodesC(i)*[1 .5 .5]));
            end

            set(graphMap.gr,'LayoutScale',0.75);
            view(graphMap.gr);
        end
        
        function structurePlot ( graphMap )
            nodes_ids = cellfun (@num2str,num2cell([1:size(graphMap.adj,1)]),'un',0);
            graphMap.gr = biograph(graphMap.adj, nodes_ids, 'LayoutType', 'radial');
            for i = 1:1:numel(graphMap.gr.edges)
                nodes = regexp(graphMap.gr.edges(i).ID,' -> ', 'split');
                nn = cellfun(@str2num,nodes);
                set (graphMap.gr.edges(i), 'LineColor', [0.85 0.85 0.85]);
                set (graphMap.gr.edges(i), 'LineWidth', 1);
            end
            graphMap.nodesC = (graphMap.nodesC-min(graphMap.nodesC))./(max(graphMap.nodesC)-min(graphMap.nodesC));
            for i = 1:1:numel(graphMap.gr.nodes)
                node = str2double(graphMap.gr.nodes(i).ID);
                set (graphMap.gr.nodes(node), 'Color', [1 1 1],...
                    'LineWidth', 3, 'Shape', 'circle',...
                    'Size', [10 10]);
                set (graphMap.gr.nodes(node), 'LineColor', [1 .5 .5]);
            end
            set(graphMap.gr,'LayoutScale',0.75);
            view(graphMap.gr);
        end
        
        function nodesCountOnlyPlot ( graphMap )
            nodes_ids = cellfun (@num2str,num2cell([1:size(graphMap.adj,1)]),'un',0);
            graphMap.gr = biograph(graphMap.adj, nodes_ids, 'LayoutType', 'radial');
            for i = 1:1:numel(graphMap.gr.edges)
                nodes = regexp(graphMap.gr.edges(i).ID,' -> ', 'split');
                nn = cellfun(@str2num,nodes);
                set (graphMap.gr.edges(i), 'LineColor', [0.85 0.85 0.85]);
                set (graphMap.gr.edges(i), 'LineWidth', 1);
            end
            graphMap.nodesC = (graphMap.nodesC-min(graphMap.nodesC))./(max(graphMap.nodesC)-min(graphMap.nodesC));
            for i = 1:1:numel(graphMap.gr.nodes)
                node = str2double(graphMap.gr.nodes(i).ID);
                set (graphMap.gr.nodes(node), 'Color', [1 1 1],...
                    'LineWidth', 3, 'Shape', 'circle',...
                    'Size', ((graphMap.nodesC(i))*10+8)*[1 1]);
                set (graphMap.gr.nodes(node), 'LineColor', 1-(graphMap.nodesC(i)*[1 .5 .5]));
            end
            set(graphMap.gr,'LayoutScale',0.75);
            view(graphMap.gr);
        end
        

    end
end