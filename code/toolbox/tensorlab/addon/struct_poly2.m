function [x,state] = struct_poly2(z, task, t)
%STRUCT_BORI short discription
%   long description
%   
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2015/11/10   NV      Initial version
    
    state = [];
    
    if isempty(task)
        x = nan(size(t,1), size(z,1));
        for r = 1:size(z,1)
            x(:,r) = polyval(z(r,:), t(:,r));
        end
    elseif ~isempty(task.r)
        x = nan(size(t,1), size(z,1));
        for r = 1:size(z,1)
            x(:,r) = polyval(task.r(r,:), t(:,r));
        end
    elseif ~isempty(task.l)
        x = nan(size(z));
        for r = 1:size(z,1)
            x(r,:) = task.l(:,r).'*conj(bsxfun(@power, t(:,r), size(z,2)-1:-1:0));
        end
    end
end
