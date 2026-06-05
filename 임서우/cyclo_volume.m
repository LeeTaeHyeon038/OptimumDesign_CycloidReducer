function V = cyclo_volume(x, p)
% CYCLO_VOLUME  Volume of the cycloid gear, proxy for reducer size (Wang Eq.7).
%   V = disk (within dedendum circle) - centre hole - pin holes + teeth
%   x = [Dp drp B D K1 Dw dsw]
    Dp = x(1); drp = x(2); B = x(3); D = x(4); K1 = x(5); dsw = x(7);

    et = K1*Dp/p.zp;        % recurring term  K1*Dp/zp

    V = (pi/4)*B*( (Dp - et - drp)^2 ...           % disk
                 - p.zw*(dsw + 2*p.Delta2 + et)^2 ...   % zw pin holes
                 - D^2 ) ...                            % centre hole
        + et*p.zc*B;                                    % teeth
end
