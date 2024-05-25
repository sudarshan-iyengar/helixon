% Function to split positive and negative pressures
function [p_pos_full, p_neg_full] = split_pos_neg(p)
    p_pos_full = zeros(size(p));
    p_neg_full = zeros(size(p));

    p_pos_indices = find(p > 0);
    p_neg_indices = find(p < 0);

    p_pos_full(p_pos_indices) = p(p_pos_indices);
    p_neg_full(p_neg_indices) = p(p_neg_indices);

    p_neg_full = abs(p_neg_full);
end
