function [Sy, Sz , ExtChi, ExtTheta] = order_parameters(Q)
% Order parameters from order matrix Q ([qx2, qy2, qy2, qxqy, qyqz, qzqx])
    ibuf = repmat({':'}, 1, numel(size(Q))-1); % python >>> matlab
    Sy = sqrt((Q(ibuf{:},1)-Q(ibuf{:},3)).^2+(2*Q(ibuf{:},6)).^2);
    Sz = sqrt((Q(ibuf{:},1)-Q(ibuf{:},2)).^2+(2*Q(ibuf{:},4)).^2);
    ExtChi = 0.5*atan2(2*Q(ibuf{:},4), (Q(ibuf{:},1)-Q(ibuf{:},2)));
    ExtTheta = pi/2-0.5*atan2(2*Q(ibuf{:},6), (Q(ibuf{:},1)-Q(ibuf{:},3)));
end