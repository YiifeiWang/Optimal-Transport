function [x] = LP_l1_proj(x, mu)
% This function projects a given sequence x0 to positive l1 ball, 
% under the constraint that the sum of the sequence is mu(>0).

[m, n] = size(x);
xx = sort(x, 'descend');
mul = cumsum(ones(m, 1));
xx_cum = cumsum(xx);
judge = sum((mul .* xx - xx_cum + mu) > 0); 
judge = max(judge, 1);
index = sub2ind([m, n], judge, (1:n));
mid = (xx_cum(index) - mu)./judge;
% diff = sum(max(x+mid, 0)) - mu;
% while norm(start, mid) > 0
%     start = (diff>0) .* start + (diff<0) .* mid;
%     last = (diff<0) .* last + (diff>0) .* mid;
%     mid = (start+last)/2;
%     diff = sum(max(x+mid, 0)) - mu;
% end
x = max(x-mid, 0);

end

