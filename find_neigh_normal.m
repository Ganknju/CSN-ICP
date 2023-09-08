function [vi] = find_neigh_normal(normal_base,vi0,tree_source,tree_normal,radius)
tree = KDTreeSearcher(tree_source);
vi = zeros(length(vi0),1);
parfor i=1:length(vi0)
    query = [tree_source(vi0(i),1) tree_source(vi0(i),2) tree_source(vi0(i),3)];
    [idx,~] = rangesearch(tree, query, radius);
    idxs = idx{1};
    center_normal = zeros(length(idxs),3);
    center_normal(:,1) = normal_base(i,1) * ones(length(idxs),1);
    center_normal(:,2) = normal_base(i,2) * ones(length(idxs),1);
    center_normal(:,3) = normal_base(i,3) * ones(length(idxs),1);
    bias_tree_normal = tree_normal(idxs,:) - center_normal;
    neigh_normal_bias = bias_tree_normal(:,1).^2 + bias_tree_normal(:,2).^2 + bias_tree_normal(:,3).^2;
   [~,mini] = min(neigh_normal_bias);
    if (~isempty(neigh_normal_bias))
    vi(i,1) = idxs(mini);
    else
    vi(i,1) = vi0(i);
    end
end
end