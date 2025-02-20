function acde_clustering()
    close all; clc;

    load fisheriris; 
    data = meas;

    option.npop   = 50;   
    option.maxit  = 500;  
    option.Kmax   = 5;    
    option.Beta   = 0.8;  
    option.Cr     = 0.7;  
    option.eps    = 1e-6; 
    option.gamma  = 10;
    [option.s, option.d] = size(data);

    num_runs = 10;
    % array
    all_final_centers   = cell(num_runs, 1);
    all_final_clusters  = cell(num_runs, 1);
    all_final_activeKs  = zeros(num_runs, 1);

    for run_id = 1:num_runs
        pop = initialize_population(option, data);

        for it = 1:option.maxit
            new_pop = pop;
            for i = 1:option.npop
                U = mutation(pop, i, option);
                Child = crossover(pop(i), U, option);
                Child = repair_chromosome(Child, option, data);
                Child = fitness(Child, option, data);
                if Child.cost < pop(i).cost
                    new_pop(i) = Child;
                end
            end
            pop = new_pop;
        end

        [~, best_idx] = min([pop.cost]);
        best_sol = pop(best_idx);

        clusters = assign_clusters(data, best_sol.centers);
        empty_mask = false(size(best_sol.centers,1),1);

        for k = 1:size(best_sol.centers,1)
            if sum(clusters == k) == 0
                empty_mask(k) = true;
            end
        end

        final_centers = best_sol.centers(~empty_mask,:);
        old2new = cumsum(~empty_mask);
        new_clusters = arrayfun(@(x) old2new(x), clusters);

        all_final_centers{run_id}  = final_centers;
        all_final_clusters{run_id} = new_clusters;
        all_final_activeKs(run_id) = size(final_centers, 1);

        fprintf('Run %d finished. Active Clusters = %d\n', ...
            run_id, all_final_activeKs(run_id));
    end

    figure('Name','ACDE Clustering Results (10 Runs)');
    for run_id = 1:num_runs
        subplot(2,5,run_id);
        gscatter(data(:,1), data(:,2), all_final_clusters{run_id});
        hold on;
        centers_now = all_final_centers{run_id};
        scatter(centers_now(:,1), centers_now(:,2), 150, 'kx','LineWidth',2);
        title(['Run ' num2str(run_id) ' | K=' num2str(all_final_activeKs(run_id))]);
        hold off; grid on;
    end
end

% done
function pop = initialize_population(option, data)
    vector_length = option.Kmax*(option.d+1);
    % array
    pop = repmat(struct('vector',[], 'centers',[], 'thresholds',[], ...
                        'activeK',0, 'cost',inf), option.npop,1);

    for i = 1:option.npop
        vec = zeros(vector_length,1);
        for k = 1:option.Kmax
            for dim = 1:option.d
                mn = min(data(:,dim));
                mx = max(data(:,dim));
                vec((k-1)*(option.d+1) + dim) = mn + rand()*(mx - mn);
            end
            vec(k*(option.d+1)) = rand();
        end
        pop(i).vector = vec;
        pop(i) = fitness(pop(i), option, data);
    end
end

% done
function U = mutation(pop, i, option)
    r = randsample(setdiff(1:option.npop, i), 3, false);     
    U.vector = pop(r(1)).vector + option.Beta*(pop(r(2)).vector - pop(r(3)).vector);
end

% done
function Child = crossover(parent, U, option)
    Child.vector = parent.vector;
    n = option.Kmax * option.d +1;
    j_rand = randi(n);
    for j = 1:n
        if rand() < option.Cr || j == j_rand
            Child.vector(j) = U.vector(j);
        end
    end
end

% done
function Child = repair_chromosome(Child, option, data)
    for k = 1:option.Kmax
        for dim = 1:option.d
            idx = (k-1)*(option.d+1)+dim;
            mn = min(data(:,dim));
            mx = max(data(:,dim));
            if Child.vector(idx) < mn, Child.vector(idx) = mn; end
            if Child.vector(idx) > mx, Child.vector(idx) = mx; end
        end
        tdx = k*(option.d+1);
        if Child.vector(tdx) < 0, Child.vector(tdx) = 0; end
        if Child.vector(tdx) > 1, Child.vector(tdx) = 1; end
    end
    
    thr = Child.vector((1:option.Kmax)*(option.d+1));
    on = (thr > 0.5);
    nOn = sum(on);
    if nOn == 0
        r = randperm(option.Kmax,2);
        Child.vector(r*(option.d+1)) = 0.6 + 0.4*rand(2,1);

    elseif nOn == 1
        offidx = find(~on);
        if ~isempty(offidx)
            choice = offidx(randi(numel(offidx)));
            Child.vector(choice*(option.d+1)) = 0.6 + 0.4*rand();
        end
    end
end

% done
function ind = fitness(ind, option, data)
% matrix
    c = decode_centers(ind.vector, option.Kmax, option.d);
    t = decode_thresholds(ind.vector, option.Kmax, option.d);
    mask = (t > 0.5);
    cA = c(mask,:);
    nA = sum(mask);

    if nA == 0
        ind.cost = 1e10; 
        ind.centers = []; 
        ind.thresholds=[]; 
        ind.activeK=0;
        return
    end

    clus = assign_clusters(data, cA);
    for kk = 1:nA
        if sum(clus==kk) < 2
            cA(kk,:) = mean(data,1);
        end
    end

    dbv = db_index(data, clus, cA);
    cm = 1/(dbv+option.eps);
    ind.cost = -cm + option.gamma*(option.Kmax - nA);
    
    ind.centers=cA; 
    ind.thresholds=t; 
    ind.activeK=nA;
end

% done
function val = db_index(data, clusters, centers)
    K = size(centers,1);
    if K < 2, val = 100; return; end
    % array
    s = zeros(K,1);
    for k = 1:K
        pts = data(clusters==k,:);
        if ~isempty(pts)
            s(k) = mean(sqrt(sum((pts - centers(k,:)).^2,2)));
        end
    end
    % matrix
    r = zeros(K); 
    for i=1:K
        for j=i+1:K
            d = norm(centers(i,:)-centers(j,:));
            if d~=0
                r(i,j) = (s(i)+s(j))/d;
                r(j,i) = r(i,j);
            end
        end
    end
    val = mean(max(r,[],2)); 
    if isnan(val), val=0; end
end

% done
function clusters = assign_clusters(data, centers)
    if isempty(centers) 
        clusters = []; return;
    end
    %array
    [~, clusters] = min(pdist2(data,centers),[],2);
end

% done
function c = decode_centers(v, Kmax, d)
    c = zeros(Kmax,d);
    for k=1:Kmax
        st=(k-1)*(d+1)+1; en=(k-1)*(d+1)+d;
        c(k,:)=v(st:en)';
    end
end

% done
function t = decode_thresholds(v,Kmax,d)
    t=zeros(Kmax,1);
    for k=1:Kmax
        t(k)=v((k-1)*(d+1)+d+1);
    end
end
