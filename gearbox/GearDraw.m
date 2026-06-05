% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function y=GearDraw(xfin, xp)
    zs=15;
    zk=16;
    t=linspace(0, 2*pi, 600);
    xc=(xfin(1)*(zs+1)/xfin(3))*cos(t)+xfin(1)*cos(zk*t)-xfin(2)*cos(t+atan(sin(zs*t)./(1/xfin(3)+cos(zs*t))));
    yc=(xfin(1)*(zs+1)/xfin(3))*sin(t)+xfin(1)*sin(zk*t)-xfin(2)*sin(t+atan(sin(zs*t)./(1/xfin(3)+cos(zs*t))));
    plot(xc, yc, '-b');
    for i=1:8
        alpha=(i-1)*2*pi/8;
        xpp=xfin(5)*cos(t);
        ypp=xfin(5)*sin(t)+xfin(7);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '--b');
    end
    for i=1:8
        alpha=(i-1)*2*pi/8;
        xpp=xfin(6)*cos(t);
        ypp=xfin(6)*sin(t)+xfin(7);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '-b');
    end
    for i=1:16
        alpha=(i-1)*2*pi/16;
        xpp=xfin(2)*cos(t);
        ypp=xfin(2)*sin(t)+xfin(1)*(zs+1)/xfin(3);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '-b');
    end  
    
    %Draw gear before optimization
    xc=(xp(1)*(1+zs)/xp(3))*cos(t)+xp(1)*cos(16*t)-xp(2)*cos(t+atan(sin(zs*t)./(1/xp(3)+cos(zs*t))));
    yc=(xp(1)*(1+zs)/xp(3))*sin(t)+xp(1)*sin(16*t)-xp(2)*sin(t+atan(sin(zs*t)./(1/xp(3)+cos(zs*t))));
    hold on;
    plot(xc, yc, '-r');
    for i=1:8
        alpha=(i-1)*2*pi/8;
        xpp=xp(5)*cos(t);
        ypp=xp(5)*sin(t)+xp(7);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '--r');
    end
    for i=1:8
        alpha=(i-1)*2*pi/8;
        xpp=xp(6)*cos(t);
        ypp=xp(6)*sin(t)+xp(7);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '-r');
    end
    for i=1:16
        alpha=(i-1)*2*pi/16;
        xpp=xp(2)*cos(t);
        ypp=xp(2)*sin(t)+xp(1)*(zs+1)/xp(3);
        A=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[xpp;ypp];
        xpp=A(1,:);
        ypp=A(2,:);
        hold on;
        plot(xpp, ypp, '-r');
    end        
end