classdef AntDiff_fixN < handle
    properties
        ants = struct;
        iterations = []; % numero de iterações
        antsN = [];
        trajeto = sparse([]); % guarda o trajeto de cada formiga ao longo da simulação
        counter = 0;
    end
    methods
        function run_AntDiff_fixN (AntDiff_fixN, graphMap)
            AntDiff_fixN.setAnts ( AntDiff_fixN.antsN );
            AntDiff_fixN.initializePosition ( graphMap.source );
            AntDiff_fixN.trajeto = sparse([]);
            
            for it = 1:1:AntDiff_fixN.iterations

                for i = 1:1:size(AntDiff_fixN.ants.pos(:,2),1)
                    AntDiff_fixN.trajeto(i,it) = AntDiff_fixN.ants.pos(i,2);
                end
                
                % Sort new position for outgoing ants
                [ possible_moves, moves_p ] = AntDiff_fixN.lookForPossibleMoves ( graphMap );
                AntDiff_fixN.ants.pos(:,3) = AntDiff_fixN.selectNextPosition ( moves_p, possible_moves );

                % move all ants to new positions
                AntDiff_fixN.updateAntsPosition();
                
            AntDiff_fixN.counter = AntDiff_fixN.counter + 1;
            end            
        end
        
        % Settings
        function setAntDiff_fixN ( AntDiff_fixN, ants_number, iterations )
            AntDiff_fixN.antsN = ants_number;
            AntDiff_fixN.iterations = iterations;
        end
        
        function setAnts ( AntDiff_fixN, ants_number )
            AntDiff_fixN.ants.pos = zeros (ants_number,3)-1;
            AntDiff_fixN.ants.tourL = zeros (ants_number,1);            
        end
        function updateAntsPosition ( AntDiff_fixN )
            AntDiff_fixN.ants.pos(:,1) = AntDiff_fixN.ants.pos(:,2);
            AntDiff_fixN.ants.pos(:,2) = AntDiff_fixN.ants.pos(:,3);
            AntDiff_fixN.ants.pos(:,3) = zeros(size(AntDiff_fixN.ants.pos,1),1)-1;
            %l = pos(:,2).*size(edgesLen,1) - (size(edgesLen,1)-pos(:,1));
            %tourl = tourl + edgesLen(l);
        end
        function [nextP] = selectNextPosition ( AntDiff_fixN, moves_p, possible_moves )
            moves_p = cumsum(moves_p')' ./  ( sum(moves_p')' * ones(1,size(moves_p,2)) );
            moves_p(possible_moves==0) = 0;
            r = rand(size(moves_p,1),1);
            (((moves_p > r*ones(1,size(moves_p,2)))')'==1);
            nextP = repmat([1:1:size(moves_p,2)],size(moves_p,1),1);
            nextP = (nextP.* (((moves_p > r*ones(1,size(moves_p,2)))')'==1));
            nextP(nextP==0)=NaN;
            nextP = min(nextP')';
        end
        function [possible_moves, moves_p] = lookForPossibleMoves ( AntDiff_fixN, graphMap )
            possible_moves = zeros ( size(AntDiff_fixN.ants.pos,1), size(graphMap.adj,2) );
            moves_p = possible_moves;
            for a = 1:1:size(AntDiff_fixN.ants.pos,1)
                [cx] = AntDiff_fixN.getAntCPos ( a );
                [px] = AntDiff_fixN.getAntPPos ( a );

                % get possible moves
                [cx] = AntDiff_fixN.getAntCPos ( a );
                [px] = AntDiff_fixN.getAntPPos ( a );
                possible_moves(a,:) = graphMap.adj(cx,:);
                possible_moves(a,cx) = 0; %probabilidade de ficar no mesmo ponto
                %if ( sum (possible_moves(a,:)) > 2 )
                %    possible_moves(a,px) = 0; %allow return to previous only if its in a dead end node
                %end
                
                % calc probability of each
                if (sum(possible_moves(a,:))>0)
                    moves_p(a,:) = possible_moves(a,:);
                    moves_p(a,:) = moves_p(a,:)./sum(moves_p(a,:));
                end                
            end
            
        end
        % Utils
        function initializePosition ( AntDiff_fixN, source_node )
            AntDiff_fixN.ants.pos(:,1) = source_node;
            AntDiff_fixN.ants.pos(:,2) = source_node;
            AntDiff_fixN.ants.pos(:,3) = -1;
        end
        function [x] = getAntCPos ( AntDiff_fixN, i )
            x = AntDiff_fixN.ants.pos(i,2);
        end
        function [x] = getAntPPos ( AntDiff_fixN, i )
            x = AntDiff_fixN.ants.pos(i,1);
        end
        function printAnts ( AntDiff_fixN )
            AntDiff_fixN.ants.pos            
            AntDiff_fixN.ants.mode
        end

    end
end
        
        

        
    