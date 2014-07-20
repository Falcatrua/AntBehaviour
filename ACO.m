classdef ACO < handle
    properties
        ants = struct;
        rho = []; % taxa de evaporação de feromônio
        Q = 20; % quantidade de feromônio base a ser depositada
        alpha = 1;
        beta = 0.5;
        antsN = []; % numero de formigas
        iterations = []; % numero de iterações

        trajeto = sparse([]); % guarda o trajeto de cada formiga ao longo da simulação
        counter = 0;
    end
    methods
        % Run
        function run_ACO ( ACO, graphMap )
            ACO.setAnts ( 1 );
            ACO.initializePosition ( graphMap.source );
            ACO.trajeto = sparse([]);

            for it = 1:1:ACO.iterations

                for i = 1:1:size(ACO.ants.pos(:,2),1)
                    ACO.trajeto(i,it) = ACO.ants.pos(i,2);
                end

                % Sort new position for outgoing ants
                %[ possible_moves, moves_p ] = lookForPossibleMoves ( ants, map, aco );
                %ants.pos(:,3) = selectNextPosition ( moves_p, possible_moves );
                recruit_count = 0;

                % Set mode
                for a = 1:1:size(ACO.ants.mode,1)
                    ACO.ants.path{a} = ACO.ants.path{a}(:);
                    graphMap.nodesC(ACO.ants.pos(a,2)) = graphMap.nodesC(ACO.ants.pos(a,2))+1;
                    if (ACO.ants.mode(a,1)==0)
                        ACO.ants.path{a}(end+1) = ACO.ants.pos(a,2); % updating path so far
                    end

                    % set mode to backwards if ant just got to the feeder
                    if (ACO.ants.pos(a,2)==graphMap.feeder)
                        ACO.ants.mode(a,1) = 1;
                        % generate ants deterministic backwards path without loops
                        ACO.ants.path{a} = ACO.backPath ( ACO.ants.path{a} );
                        ACO.ants.tourL(a,1) = size ( ACO.ants.path{a},1 );%por enquanto, na real tem que ser a soma de map.edgesLen do path
                    end
                    if ( ACO.ants.pos(a,2) == graphMap.source )
                        ACO.ants.mode(a,1) = 0;
                        ACO.ants.path{a} = [graphMap.source];
                        recruit_count = recruit_count+1;
                    end

                end

                % Sort new position for outgoing ants
                [ possible_moves, moves_p ] = ACO.lookForPossibleMoves ( graphMap );
                ACO.ants.pos(:,3) = ACO.selectNextPosition ( moves_p, possible_moves );

                % edges Pheromone update
                ACO.updateEdgesPhMarking ( graphMap  );
                graphMap.edgesPh = graphMap.symm ( graphMap.edgesPh );

                % move all ants to new positions
                ACO.updateAntsPosition();



                %s = sprintf ('t = %d; ants = %d; rc = %d; lambda = %f', it, size(ants.pos(:,2),1), recruit_count, exp( -0.25 * recruit_count ));
                %display(s);
                pvs(it,:) = [it, size(ACO.ants.pos(:,2),1), recruit_count, exp( -0.25 * recruit_count )];
                if (recruit_count > 0)
                    for i = 1:1:recruit_count
                        lambda = exp( -0.25 * recruit_count );
                        ACO.recruit (lambda, graphMap.source);
                    end
                end
            ACO.counter = ACO.counter + 1;
            end
            
        end
        % Settings
        function setAco ( ACO, ants_number, iterations )
            ACO.Q = 20;
            ACO.rho = 0.05; %pheromone evaporation rate
            ACO.alpha = 1;
            ACO.beta = 0.5;
            ACO.antsN = ants_number;
            ACO.iterations = iterations;
        end
        function setAnts ( ACO, ants_number )
            ACO.ants.pos = zeros (ants_number,3)-1;
            ACO.ants.tourL = zeros (ants_number,1);
            ACO.ants.mode = zeros(ants_number,1); %0 = forward mode, 1 = backward mode;
            ACO.ants.path = cell (ants_number,1);%percurso até o alimento para poder gerar o retorno deterministico
        end
        % Algorithm
        function updateEdgesPhMarking ( ACO, graphMap )
            for a = 1:1:size(ACO.ants.pos(:,1),1)
                if (ACO.ants.mode(a,1)==1)
                    e = ACO.ants.pos(a,2); i = ACO.ants.pos(a,3);
                    graphMap.edgesPh(e,i) = graphMap.edgesPh(e,i) + (ACO.Q /ACO.ants.tourL(a,1));
                    % tem que simetrizar agora porque a função symm não vai pegar
                    % pelo maior valor
                    graphMap.edgesPh(i,e) = graphMap.edgesPh(i,e) + (ACO.Q /ACO.ants.tourL(a,1));
                    %display ('marcar:');
                    %[e i (ACO.Q /ACO.ants.tourL(a,1));...
                    %    graphMap.edgesPh(e,i) graphMap.edgesPh(i, e) 0] %evento de marcação
                end
            end
            graphMap.edgesPh = graphMap.edgesPh.*(1-ACO.rho);
        end
        function [nl_path] = backPath ( ACO, path )
            if (size(path,2)>1)
                sprintf ('\n\nError\nbackPath:: path has more than one column');
            end
            path = path(:);
            loop = 1;
            nl_path = path;
            while ( loop > 0 )
                loop = 0;
                for i = 1:1:length(path)
                    for j = i+1:1:length(path)
                        if (path(i)==path(j))
                            loop = 1;
                            nl_path = [path(1:i); path(j+1:end)];
                        end
                    end
                end
                path = nl_path;
            end
        end
        function updateAntsPosition ( ACO )
            ACO.ants.pos(:,1) = ACO.ants.pos(:,2);
            ACO.ants.pos(:,2) = ACO.ants.pos(:,3);
            ACO.ants.pos(:,3) = zeros(size(ACO.ants.pos,1),1)-1;
            %l = pos(:,2).*size(edgesLen,1) - (size(edgesLen,1)-pos(:,1));
            %tourl = tourl + edgesLen(l);

        end
        function [nextP] = selectNextPosition ( ACO, moves_p, possible_moves )
            moves_p = cumsum(moves_p')' ./  ( sum(moves_p')' * ones(1,size(moves_p,2)) );
            moves_p(possible_moves==0) = 0;
            r = rand(size(moves_p,1),1);
            (((moves_p > r*ones(1,size(moves_p,2)))')'==1);
            nextP = repmat([1:1:size(moves_p,2)],size(moves_p,1),1);
            nextP = (nextP.* (((moves_p > r*ones(1,size(moves_p,2)))')'==1));
            nextP(nextP==0)=NaN;
            nextP = min(nextP')';
        end
        function [possible_moves, moves_p] = lookForPossibleMoves ( ACO, graphMap )
            possible_moves = zeros ( size(ACO.ants.pos,1), size(graphMap.adj,2) );
            moves_p = possible_moves;
            for a = 1:1:size(ACO.ants.pos,1)
                [cx] = ACO.getAntCPos ( a );
                [px] = ACO.getAntPPos ( a );

                %for out going ants
                if (ACO.ants.mode(a)==0)
                    % get possible moves
                    [cx] = ACO.getAntCPos ( a );
                    [px] = ACO.getAntPPos ( a );
                    possible_moves(a,:) = graphMap.adj(cx,:);
                    possible_moves(a,cx) = 0;
                    if ( sum (possible_moves(a,:)) > 2 )
                        possible_moves(a,px) = 0; %allow return to previous only if its in a dead end node
                    end

                    % calc probability of each
                    if (sum(possible_moves(a,:))>0)
                        moves_p(a,:) = (graphMap.edgesPh(cx,:).^ACO.alpha);
                        moves_p(a,:) = moves_p(a,:).*possible_moves(a,:); %exclude possible_moves==0
                        moves_p(a,:) = moves_p(a,:)./sum(moves_p(a,:));
                    end

                else % returning ants (deterministic, as in path)
                    possible_moves(a,ACO.ants.path{a}(end-1)) = 1;
                    moves_p(a,ACO.ants.path{a}(end-1)) = 1;
                    if (size(ACO.ants.path{a},1)>1)
                        ACO.ants.path{a} = ACO.ants.path{a}(1:size(ACO.ants.path{a},1)-1);
                    else
                        ACO.ants.path{a} = [];
                        ACO.ants.mode(a) = 0;
                    end

                end


            end

        end
        % Utils
        function initializePosition ( ACO, source_node )
            ACO.ants.pos(:,1) = source_node;
            ACO.ants.pos(:,2) = source_node;
            ACO.ants.pos(:,3) = -1;
        end
        function [x] = getAntCPos ( ACO, i )
            x = ACO.ants.pos(i,2);
        end
        function [x] = getAntPPos ( ACO, i )
            x = ACO.ants.pos(i,1);
        end
        function printAnts ( ACO )
            ACO.ants.pos
            ACO.ants.path
            ACO.ants.mode
        end
        function addAnts ( ACO, n_to_add, source_node )
            ACO.ants.pos = [ ACO.ants.pos; [ones(n_to_add,2)*source_node, zeros(n_to_add,1)-1] ];
            ACO.ants.tourL = [ ACO.ants.tourL; zeros(n_to_add,1)];
            ACO.ants.mode = [ ACO.ants.mode; zeros(n_to_add,1)];
            ACO.ants.path = [ ACO.ants.path; cell(n_to_add,1)];
        end
        function recruit ( ACO, lambda, source_node)
            n_to_add = poissrnd ( lambda );
            ACO.addAnts ( n_to_add, source_node );
        end

    end
end