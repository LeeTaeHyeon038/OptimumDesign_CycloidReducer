function V = GearObj_fmincon(x)
    zs=15; zk=16; zi=8; N=1000;
    e=x(1); q=x(2); m=x(3); h=x(4); Rs=x(5); Rh=x(6);
    r = e * zk / m;
    V = GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
    if ~isfinite(V) || V < 1000
        V = 1e9;
    end
end