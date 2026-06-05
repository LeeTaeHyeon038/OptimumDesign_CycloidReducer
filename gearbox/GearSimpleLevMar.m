% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function xnew=GearSimpleLevMar(nmax, obj, jac, x0)
    k=10;
    n=0;
    err=1e-10;
    lambda=1000;
    jacTjac=jac(x0)'*jac(x0);
    grad=jac(x0)'*obj(x0);
    xprev=x0;
    xnext=x0-(jacTjac+lambda*diag(diag(jacTjac)))\grad;
    while n<nmax && norm(xprev-xnext)>1e-6
        fnext=0.5*norm(obj(xnext))*norm(obj(xnext));
        fprev=0.5*norm(obj(xprev))*norm(obj(xprev));
        if (fprev-fnext)> err
            n=n+1;
            xprev=xnext;
            lambda=lambda/k;
            jacTjac=jac(xprev)'*jac(xprev);
            grad=jac(xprev)'*obj(xprev);
            xnext=xprev-(jacTjac+lambda*diag(diag(jacTjac)))\grad;
        else
            n=n+1;
            lambda=lambda*k;
            jacTjac=jac(xprev)'*jac(xprev);
            grad=jac(xprev)'*obj(xprev);
            xnext=xprev-(jacTjac+lambda*diag(diag(jacTjac)))\grad;
        end
    end 
    xnew=xnext;
end