% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function xnew=GearSteepestDescent(nmax, obj, jac, x0)
    n=0;
    err=1e-10;
    lambda=1e-13;
    grad=jac(x0)'*obj(x0);
    xprev=x0;
    xnext=x0-lambda*grad;
    while n<nmax && norm(xprev-xnext)>1e-6
        fnext=0.5*norm(obj(xnext))*norm(obj(xnext));
        fprev=0.5*norm(obj(xprev))*norm(obj(xprev));
        if (fprev-fnext)> err
            n=n+1;
            xprev=xnext;
            grad=jac(xprev)'*obj(xprev);
            xnext=xprev-lambda*grad;
        else 
            grad=jac(xprev)'*obj(xprev);
            xnext=xprev-lambda*grad;
            break;
        end
    end 
    xnew=xnext;
end