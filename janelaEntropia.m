function [H] = janelaEntropia (t, trajeto)
%%janela temporal deslizante de entropia para a trajetória
j = t;
H = [];
for i = 1:1:size(trajeto,2)-j
    c = []; g = [];
    g = grp2idx ( nonzeros(trajeto(:,i:i+j)) );
    for e = 1:1:max(g)
        c(e,1) = sum (g==e);
    end
    c = c./sum(c);
    c = -log2(c).*c;
    c(isnan(c)) = 0;
    H(i) = sum(c);
end

end

%{
H = [];
HT = [];

for i = 1:1:size(trajeto,2)-t
    c = []; g = []; h = [];
    for a = 1:1:size(trajeto,1);
        g = grp2idx ( trajeto(i:i+t,a) );
        for e = 1:1:max(g)
            c(e,1) = sum (g==e);
        end
        c = c./sum(c);
        c = -log2(c).*c;
        h = sum(c);
        H(i,a) = h;
    end
    c = []; g = []; h = []; T = [];
    T = trajeto(i:i+t,:);
    T = T(:);
    g = grp2idx ( T );
    for e = 1:1:max(g)
        c(e,1) = sum (g==e);
    end
    c = c./sum(c);
    c = -log2(c).*c;
    h = sum(c);
    HT(i,1) = h;
    
    %plotHT (HT)
end
%}
%{
function plotHT (HT)
    figure; hold on;
    for i = 1:1:numel(sim)
        n = length(HT);
    for ii = 1:1:n
        if (ii<n)
            plot (i,sim{i}.ht(ii),'.','color',1-(ii/n)*[1 1 1],'markerSize',10);
        else
            plot (i,sim{i}.ht(ii),'r.','markerSize',10);
        end
    end
end
hold off;
xlim([0 numel(sim)+1]); 
clear ('i','n','ii');
%}    