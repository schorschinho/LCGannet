function [xCurrent,Resnorm,fval,exitflag,extra_arguments,output,lambda,Jx] =LM(obj,xGuess,lb,ub,opt)
%%  LEVENBERG-MARQUARD ALGORITHM WITH BROYDEN UPDATES FOR THE JACOBIAN
%
%   LevenbergMarquardt is similar to lsqnonlin with the levenberg-marquardt
%   algorithm with the three main advantages:
%
%   1)  the jacobian can be updated using the Broyden method which 
%       minimizes function evaluations
%   2)  the variables are transformed to implement box-constraints, and
%   3)  function arguments can be passed on.
%
%   Some functionality is very similar, others differ slightly. 
%
%   syntax:
%
%   xCurrent=LevenbergMarquardt(obj,xGuess);
%   xCurrent=LevenbergMarquardt(obj,xGuess,lb);
%   xCurrent=LevenbergMarquardt(obj,xGuess,lb,ub);
%   xCurrent=LevenbergMarquardt(obj,xGuess,lb,ub,opt);
%
%   where obj and xGuess are the objective function and the inital guess.
%   lb and ub are the box constraints and can be left empty, or set to
%   infinity or -infinity, respectively, if no constraints are binding. opt
%   allows to pass on other arguments as explained below.
%
%   example: 
%
%   as an example we can fit a 2-degree taylor on a sin function in the
%   unit interval:
%     xGuess=ones(1,3);
%     x=linspace(0,1)';
%     curve=@(para)para(1)+para(2).*x+para(3).*x.^2;
%     obj=@(para)[sin(x)-curve(para)];
%   lets exploit some functionalities:
%     opt.Display='iter';
%     opt.title='curve fitting for sin function';
%   lets assume we are interested in how good the fit looks outside the
%   fitted range (1,2)
%     opt.pltfn=@(para)[sin(x+1)-para(1)+para(2).*(x+1)+para(3).*(x+1).^2]; 
%   lets use the jacobian method that exploits romberg extrapolation (more
%   precise generally)
%     opt.Jacobian='romberg';
%     xCurrent =LevenbergMarquardt(obj,xGuess,[],[],opt);
%     figure
%     plot(x,[sin(x) curve(xCurrent)])
%   lets restrict all parameters to be positive and shut off the broyden
%   type updates.
%     opt.Broyden_updates='off';
%     xCurrent2=LevenbergMarquardt(obj,xCurrent,zeros(1,3),[],opt);
%     plot(x,[sin(x) curve(xCurrent) curve(xCurrent2)])
%     legend('sin','curve fit','curve fit with restrictions')
%
%   There is another example attached to this file that deals explicitly
%   with argument passing.
%
%   author: Alexander.Dentler at gmail.com, October 9th 2015, all errors
%   are yours to deal with. 
%
%   version 3

%%  INITIAL PARAMETER READ IN
xForm=xGuess;
n = numel(xGuess);
xGuess=xGuess(:);
t1=cputime;
if exist('lb','var')
    if isempty(lb)
        lb=-inf(size(xForm));
    end
else
    lb=-inf(size(xForm));
end
if exist('ub','var')
    if isempty(ub)
        ub=inf(size(xForm));
    end
else
    ub=inf(size(xForm));
end
%%  DEFAULT PARAMETERS
%   display
Display='limit';             %   different degrees of notifications, default gives final output
title=[];                   %   user can add a title that will show up before the iteration (and renewals of the header) and on the figure name
IterDispRenewal=30;         %   renew header of iterative display after so many steps
OutputFcn=[];               %   collects user specified output function which is either a function handle or a cell array with functions in each cell
PlotIterations=0;           %   create plots at each iteration
pltfn=[];                   %   user can supply one function which takes x as an argument and plots it values
%   passing arguments
extra_arguments=cell(1,0);  %   arguments created with each iteration that is passed as argument in the next iteration (for nested iterations)
%   iteration and function count
MaxIter=200*n;              %   maximal no of iterations
MaxFunEvals=200*n;          %   maximal no of function calls
%   acceptance tolerances
AccTol=0;                   %   breakup optimization if Resnorm is below this value (accept function value)
FooTol=1e-7;                %   breakup optimization if gradients are below this value (first order optimality condition)
RelFooTol=1e-6;             %   breakup optimization if relative gradients are below this value (first order optimality condition)
IncrTol=1e-6;               %   new evaluated point shows enough improvement to be fully accepted and lambda does not decrease
TolFun=NaN;                 %   breakup optimization if absolute Resnorm improvement falls below this level
RelTolFun=NaN;              %   breakup optimization if relative Resnorm improvement falls below this level
TolX=1e-7;                  %   breakup optimization if all absolute changes in parameters fall below this level
RelTolX=1e-5;               %   breakup optimization if all relative changes in parameters fall below this level
%   Levenberg-Marquardt parameters
Jacobian='off';             %   user supplies jacobian as second output argument of fun if set to 'on', otherwise 'off' (default)
MaxStepNo=2;                %   no of steps taken to find jacobian when "Jacobiab" is set to 'limit' or 'extrapolation'
FinDiffRelStep=eps^(1/2);   %   multiplies step size of non-supplied Jacobian finite difference step
TypicalX=ones(n);                 %   scales step size of non-supplied Jacobian finite difference step
DerivativeCheck='off';      %   if set to 'on' it checks Jacobian in first step
ScaleProblem='none';    %   if set to 'Jacobian' the problem is rescaled, by 'none' its not.
InitDamping=1e-2;           %   initial dampening
MinDamping=1e-7;            %   minimal dampening
MaxDamping=1e20;             %   maximal dampening
FactDamping=10;             %   increases or decreases dampening in loop
MaxEigTol=1e-6;             %   if largest eigenvalue becomes smaller than this value we attempt to use contraction mapping
Broyden_updates='off';       %   set to 'on' it gives Broyden updates for the Jacobian for every 2*n steps, set to 'off' it requires updates in each iteration, 
%                               set to an integer thats the amount of steps until we update
conservative_updates=1;     %   set to 1 it will only enforce the tolerances for foo, stepsize or eigenvalue when we just updated Jacobian
%%  DYNAMIC READ IN OF STRUCTURE THAT GIVES OPTIONS
fval=[];J=[];
if nargin>=5 % optional parameter structure opt has been provided
    if isstruct(opt)
        list=who;
        %%  check specifically if bounds are also given in option structure
        %   and check if they are identical, otherwise give error
        if isfield(opt, 'lb')
            if ~isequal(opt.lb(:),lb(:)) && ~isequal(-inf(numel(xGuess),1),lb(:))
                error('Lower bounds given are not identical to lower bounds in option structure.')
            end
        end
        if isfield(opt, 'ub')
            if ~isequal(opt.ub(:),ub(:)) && ~isequal(inf(numel(xGuess),1),ub(:))
                error('Upper bounds given are not identical to upper bounds in option structure.')
            end
        end
        %%  dynamic read in
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
end
%%  ERROR CHECKING AND PREP
%   boundary error check
if (numel(lb)>0 || numel(lb)>0) && (numel(xGuess)~=numel(lb) || numel(xGuess)~=numel(ub)) || any(lb(:)>=ub(:))
    error('Number of bounds does not correspond to number of variable to optimize, or bounds are reversed.')
end
lb=lb(:);
ub=ub(:);
if any(xGuess<lb) || any(xGuess>ub)
    warning('The guess is outside the specified domain. We correct for this. But don''t let that happen again.')
    xGuess(xGuess<lb)=lb(xGuess<lb)+eps;
    xGuess(xGuess>ub)=lb(xGuess>ub)-eps;
end
%   transformations
yGuess=bnd2unbnd(xGuess,lb,ub);
transformback=@(y)unbnd2bnd(y,lb,ub);
%   objective
objx=@(x,vara)obj(reshape(x,size(xForm)),vara{:});
objy=@(y,vara)obj(reshape(transformback(y),size(xForm)),vara{:});
%   Jacobian matter
switch Jacobian
    case 'limit'
        Jacobian_method=2;
    case {'extrapolation','romberg'}
        Jacobian_method=3;
    case 'on'
        Jacobian_method=4;
    otherwise
        Jacobian_method=1;
end
option_Jacobian.Jacobian_method=Jacobian_method;
option_Jacobian.FinDiffRelStep=FinDiffRelStep;
option_Jacobian.TypicalX=TypicalX;
option_Jacobian.MaxStepNo=MaxStepNo;
option_Jacobian.lb=lb;
option_Jacobian.ub=ub;
if ~isscalar(Broyden_updates)
    if ischar(Broyden_updates)
        switch Broyden_updates
            case 'off'
                Broyden_updates=0;
            otherwise
                Broyden_updates=2*n;
        end
    end
end
%   dampening
lambda=InitDamping;lambda_old=InitDamping;
%   display
switch Display
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
    case 'simplex'
        prnt = 4;
    otherwise
        prnt = 1;
end
option_Jacobian.dsp=2*(prnt==3)+(prnt==4);
prnt_fun=@(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,how)fprintf('%4.0f  %4.0f  % 8.2g   % 8.2g   % 8.2g   % 8.2g   % 8.2g      % 8.2g          % 8.2g     % 8.2g % 8.2g % 8.2g % 8.2g   %s  \n',x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,how);
txt=['\n                                                 Norm   Rel norm   First-order   Rel first-order      Largest                              Comment\n',...
        '  I#    F#      f(x)      Df(x)   relDf(x)    of step    of step    optimality        optimality   eigenvalue   lambda      rho    ratio          '];
header_wo_title=sprintf(txt);
if ~isempty(title)
    header = sprintf(horzcat('\n',title,txt));
else
    header = header_wo_title;
end
%%  STARTING VALUES
flag_nochange=1;
h=NaN(n,1);E=NaN(n,1);
iteration=0;funccount=1;Jacobian_counter=1;how='Initial evaluation';
xLMstep=NaN(n,1);yLMstep=NaN(n,1);newtonianstep=NaN(n,1);xgradient=NaN(n,1);
rho=NaN;ratio=NaN;fun_improv=NaN;rel_fun_improv=NaN;maxFoo=NaN;maxrelFoo=NaN;
maxabsstep=NaN;maxabsrelstep=NaN;MaxEigJJ=NaN;
%   first evaluation
if ~isempty(fval) && (~isempty(J) || (isempty(J) && Jacobian_method<4))
    Resnorm=sum(fval.^2);
else
    [fval,Resnorm,J,extra_arguments]=eval_fun(objy,yGuess,Jacobian_method,extra_arguments);
end
option_Jacobian.m=numel(fval);
%   initial display
if prnt>1
    disp(header)
    prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
        MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
    prnt_first={iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxabsstep,maxabsrelstep,maxFoo,maxrelFoo,...
        MaxEigJJ,lambda,rho,ratio,how};
    prnt_fun(prnt_first{:});
end
%   good staring values?
if any(~isfinite(fval))
    disp('staring values: ')
    disp(num2str(unbnd2bnd(yGuess,lb,ub)))
    disp('fval: ')
    disp(num2str(fval))
    error('Some starting value is not finite.')
end
%   derivative check
if strcmp(DerivativeCheck,'on') && Jacobian_method==4
    option_Jacobian_check=option_Jacobian;
    option_Jacobian_check.Jacobian_method=1;
    J_check=eval_Jacobian(@(x)objx(x,extra_arguments),yGuess,fval,option_Jacobian_check);
    dh=D_unbnd2bnd(yGuess,lb,ub);
    err=bsxfun(@times,(J-J_check),1./dh);
    relerr=err./bsxfun(@plus,abs(transformback(yGuess)),abs(transformback(yGuess)))'/2;
    disp(horzcat('Derivative check gives largest error with ',...
        num2str(max(err(:))),' and largest relative error with ',num2str(max(relerr(:)))))
end
%%  DYNAMIC LOOP TO OPTIMIZE
%   initial call of output functions and plotting function
OtptFcnVl.xInit=xGuess;OtptFcnVl.fvalInit=fval;
OtptFcnVl.pltfn=pltfn;OtptFcnVl.transformback=transformback;
OtptFcnVl.xForm=xForm;OtptFcnVl.title=title;
OtptFcnVl.ResnormHist=Resnorm;

OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
OtptFcnVl.Jacobian_counter=Jacobian_counter;OtptFcnVl.fval=fval;
OtptFcnVl.Resnorm=Resnorm;OtptFcnVl.J=J;OtptFcnVl.h=h;OtptFcnVl.E=E;
OtptFcnVl.xLMstep=xLMstep;OtptFcnVl.yLMstep=yLMstep;
OtptFcnVl.newtonianstep=newtonianstep;OtptFcnVl.xgradient=xgradient;
OtptFcnVl.extra_arguments=extra_arguments;

OtptFcnVl.rho=rho;OtptFcnVl.ratio=ratio;OtptFcnVl.lambda=lambda;
OtptFcnVl.fun_improv=fun_improv;OtptFcnVl.rel_fun_improv=rel_fun_improv;
OtptFcnVl.maxFoo=maxFoo;OtptFcnVl.maxrelFoo=maxrelFoo;
OtptFcnVl.maxabsstep=maxabsstep;OtptFcnVl.maxabsrelstep=maxabsrelstep;
OtptFcnVl.MaxEigJJ=MaxEigJJ;
%   callup functions
stop=eval_outputfun(OutputFcn,transformback(yGuess),OtptFcnVl,'init');
if PlotIterations
    OtptFcnVl=eval_PlotIterations(transformback(yGuess),OtptFcnVl,'init');
end
%   loop
howJ='user-supplied Jacobian';
while Resnorm>AccTol && funccount< MaxFunEvals && iteration < MaxIter && stop==false
    %%  INITIALIZE ITERATION
    iteration=iteration+1;
    maxabsstep=NaN;maxabsrelstep=NaN;MaxEigJJ=NaN;
    rho=NaN;ratio=NaN;fun_improv=NaN;rel_fun_improv=NaN;
    %%  evaluate Jacobian at current point
    dh=D_unbnd2bnd(yGuess,lb,ub);
    if isempty(J) && Jacobian_method<4
        [J,h,E,func_evals_Jacobian]=eval_Jacobian(@(x)objx(x,extra_arguments),yGuess,fval,option_Jacobian);
        funccount=funccount+func_evals_Jacobian;
        Jacobian_counter=1;
        howJ='full Jacobian update';
    end
    how=howJ;
    if prnt==4
        disp('Jacobian:')
        disp(num2str(bsxfun(@times,J,1./dh))) %     transform back to x space
        if Jacobian_method==2 || Jacobian_method==3
            disp('step size for Jacobian:')
            disp(num2str(h))
            disp('Error size:')
            disp(num2str(E))
        end
    end
    %%  evaluate first order optimatlity
    xgradient=bsxfun(@times,-J,1./dh)'*fval;
    relFoo=xgradient./(1e-12+abs(transformback(yGuess)));
    maxFoo=max(abs(xgradient));
    maxrelFoo=max(abs(relFoo));
    if (maxFoo<FooTol || maxrelFoo<RelFooTol) && (Jacobian_counter==1 || ~conservative_updates)
        if maxFoo<FooTol && maxrelFoo<RelFooTol
            how='absolute and relative first order optimatility';
        elseif maxFoo<FooTol
            how='absolute first order optimatility';
        else
            how='relative first order optimatility';
        end
        break
    end
    %%  set up surrogate model
    JJ=J'*J;
    ygradient=-J'*fval;
    switch ScaleProblem
        case 'Jacobian'
            T=diag(diag(JJ));
        otherwise
            T=eye(size(JJ,1));
    end
    LMStepFun=@(lambda)pinv((JJ+lambda*T),eps)*ygradient;
    MaxEigJJ=max(eig(JJ));
    if MaxEigJJ<MaxEigTol && (Jacobian_counter==1 || ~conservative_updates)
        how='largest eigenvalue';
        break
    end
    %%  step size and stopping criteria for step size
    %   if step size is too small, we get out
    yLMstep=LMStepFun(lambda);
    newtonianstep=LMStepFun(0);
    xLMstep=transformback(yGuess+yLMstep)-transformback(yGuess);
    maxabsstep=norm(xLMstep,2);
    maxabsrelstep=norm(xLMstep./(1e-12+abs(transformback(yGuess))),2);
    if (maxabsstep<TolX || maxabsrelstep<RelTolX) && (Jacobian_counter==1 || ~conservative_updates)
        if maxabsstep<TolX && maxabsrelstep<RelTolX
            how='absolute and relative step size';
        elseif maxabsstep<TolX
            how='absolute step size';
        else
            how='relative step size';
        end
        break
    elseif (maxabsstep<TolX || maxabsrelstep<RelTolX)
        Jacobian_counter=Broyden_updates;
    end
    %%  evaluation of new point
    [LM_fval,LM_Resnorm,LM_J,LM_extra_arguments]=eval_fun(objy,yGuess+yLMstep,Jacobian_method,extra_arguments);
    funccount=funccount+1;
    if prnt==4
        disp('No, step size, New vs old parameters:')
        disp(num2str([ (1:n)' xLMstep transformback(yGuess)+xLMstep transformback(yGuess)]))
        disp(horzcat('Difference in norm of residuals :',num2str(LM_Resnorm-Resnorm)))
        disp(horzcat('New norm of residuals           :',num2str(LM_Resnorm)))
        disp(horzcat('Old norm of residuals           :',num2str(Resnorm)))
    end
    %%  evaluate dampening
    fun_improv=LM_Resnorm-Resnorm;
    rel_fun_improv=LM_Resnorm/Resnorm-1;
    rho=(Resnorm-LM_Resnorm)/(2*yLMstep'*(lambda*yLMstep+ygradient));
    ratio=(Resnorm-sum((fval+J*yLMstep).^2))/(Resnorm-LM_Resnorm);
    lambda_old=lambda;
    if rho>IncrTol || TolFun>fun_improv || RelTolFun>rel_fun_improv
        %   note
        how=horzcat('*',how);
        %   good evaluation, decrease dampening
        lambda=max(lambda/FactDamping,MinDamping);
        flag_nochange=0;
        %%  Jacobian update
        if Jacobian_method==4
            %%  user supplied function also updated Jacobian
            J=LM_J;
            howJ='user-supplied Jacobian';
        else
            %%  update jacobian or destroy it
            Jacobian_counter=Jacobian_counter+1;
            if Jacobian_counter>=Broyden_updates || isempty(J)%(2*n) || ~Broyden_updates
                J=[];
                howJ='full Jacobian update';
            else
                J=J+((LM_fval-fval-J*yLMstep)*yLMstep')/(yLMstep'*yLMstep);
                howJ='Broyden-type update';
            end
        end
        %%  other values
        fval=LM_fval;
        Resnorm=LM_Resnorm;
        yGuess=yGuess+yLMstep;
        extra_arguments=LM_extra_arguments;
    elseif lambda==MaxDamping && (Jacobian_counter<=2 || ~conservative_updates)
        how='dampening';
        break
    else
        %   bad evaluation, increase dampening
        if Jacobian_counter>1 && Jacobian_method<4 && lambda==MaxDamping
            J=[];
            lambda=InitDamping;
            howJ='full Jacobian update';
        elseif Jacobian_counter>1 && Jacobian_method<4
            %   quick dampening as we use a Broyden updated jacobian which
            %   is quick, but suboptimal so we dont want to waste time with
            %   large steps that are imprecise.
            lambda=min(lambda*(FactDamping^2),MaxDamping);
            howJ='quick dampening';
        else
            lambda=min(lambda*FactDamping,MaxDamping);
            howJ='soft dampening';
        end
    end    
    if Resnorm<=AccTol
        break
    end
    %%  iterative display for LM
    if prnt>2
        if mod(iteration,IterDispRenewal)==0 || IterDispRenewal==1
            disp(header)
            prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
                MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
        end
        prnt_fun(iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxabsstep,maxabsrelstep,maxFoo,maxrelFoo,...
            MaxEigJJ,lambda_old,rho,ratio,how);
    end
    %%  OUTPUT FUNCTIONS
    OtptFcnVl.ResnormHist=[OtptFcnVl.ResnormHist Resnorm];    
    OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
    OtptFcnVl.Jacobian_counter=Jacobian_counter;OtptFcnVl.fval=fval;
    OtptFcnVl.Resnorm=Resnorm;OtptFcnVl.J=J;OtptFcnVl.h=h;OtptFcnVl.E=E;
    OtptFcnVl.xLMstep=xLMstep;OtptFcnVl.yLMstep=yLMstep;
    OtptFcnVl.newtonianstep=newtonianstep;OtptFcnVl.xgradient=xgradient;    
    OtptFcnVl.rho=rho;OtptFcnVl.ratio=ratio;OtptFcnVl.lambda=lambda_old;
    OtptFcnVl.fun_improv=fun_improv;OtptFcnVl.rel_fun_improv=rel_fun_improv;
    OtptFcnVl.maxFoo=maxFoo;OtptFcnVl.maxrelFoo=maxrelFoo;
    OtptFcnVl.maxabsstep=maxabsstep;OtptFcnVl.maxabsrelstep=maxabsrelstep;
    OtptFcnVl.MaxEigJJ=MaxEigJJ;
    OtptFcnVl.extra_arguments=extra_arguments;
    stop=eval_outputfun(OutputFcn,transformback(yGuess),OtptFcnVl,'iter');
    if PlotIterations
        OtptFcnVl=eval_PlotIterations(transformback(yGuess),OtptFcnVl,'iter');
    end
end
%%  FINAL ASSIGNMENT
xCurrent=reshape(transformback(yGuess),size(xForm));
OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
if Resnorm<=AccTol
    how='FULL CONVERGENCE';
    exitflag=1;
elseif fun_improv>-TolFun || rel_fun_improv>-RelTolFun ...
        || maxFoo<FooTol || maxrelFoo<RelFooTol || maxabsstep<TolX ...
        || maxabsrelstep<RelTolX ||  MaxEigJJ<MaxEigTol ...
        || lambda==MaxDamping || rho<IncrTol
    exitflag=0;
elseif funccount==MaxFunEvals || iteration==MaxIter
    how='did not converge';
    exitflag=-1;
elseif stop
    how='Output function terminated algorithm';
    exitflag=1;
else
    exitflag=1;
end
if nargout>5
    output=OtptFcnVl;
    output.how=how;
    output.nochange=flag_nochange;
    output.maxchange=max(abs(OtptFcnVl.xInit-xCurrent(:)));
    output.xCurrent=xCurrent;
    output.time_used_in_sec=cputime-t1;
    if Jacobian_method<4
        output.algorithm='levenberg-marquardt with Broyden rank-1 updates for Jacobian';
    else
        output.algorithm='levenberg-marquardt with user-supplied Jacobian';
    end
    output.ErrJacobian=E;
    output.StepsizeJacobian=h;
    output.extra_arguments=extra_arguments;
    if nargout>7
        if isempty(J) && Jacobian_method<4
            [J,h,E,func_evals_Jacobian]=eval_Jacobian(@(x)objx(x,extra_arguments),yGuess,fval,option_Jacobian);
            funccount=funccount+func_evals_Jacobian;
        end
        dh=D_unbnd2bnd(yGuess,lb,ub);
%         Jx=bsxfun(@times,J,1./dh); 
        Jx = J;
        output.firstorder=(Jx'*fval);
        output.firstorderopt=max(abs((Jx'*fval)));
        output.Jacobian=Jx;
        output.dh=dh;
    end
end
%   close things
eval_outputfun(OutputFcn,transformback(yGuess),OtptFcnVl,'done');
if PlotIterations
    eval_PlotIterations(transformback(yGuess),OtptFcnVl,'done');
end
if prnt>0    
    if prnt==3
        if iteration>1
            disp(header)
        end
        prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
            MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
        prnt_fun(prnt_first{:})    
    elseif prnt==1
        disp(header_wo_title)
        %prnt_fun(prnt_first{:})   
    end
    prnt_fun(iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxabsstep,maxabsrelstep,maxFoo,maxrelFoo,...
        MaxEigJJ,lambda_old,rho,ratio,how);
end
end

function    [fval,Resnorm,J,extra_arguments]=eval_fun(objy,yGuess,Jacobian_method,extra_arguments)
%%  evaluate
J=[];
argu_out=cell(1,1+(Jacobian_method==4)+numel(extra_arguments));
[argu_out{:}]=objy(yGuess,extra_arguments);
fval=argu_out{1};
%   Jacobian
if Jacobian_method==4;
    J=argu_out{2};
end
%   if we got extra arguments
if numel(extra_arguments)>0
    [extra_arguments{:}]=deal(argu_out{2+(Jacobian_method==4):end});
end
Resnorm=sum(fval.^2);
end

function    [J,h,E,func_evals_Jacobian]=eval_Jacobian(objx,yGuess,fval,option_Jacobian)
%%  prep structure
lb=option_Jacobian.lb;
ub=option_Jacobian.ub;
xGuess=unbnd2bnd(yGuess,lb,ub);
Jacobian_method=option_Jacobian.Jacobian_method;
if any(lb==xGuess) || any(ub==xGuess)
    Jacobian_method=1;
end
option_Jacobian.f_0=fval;
%%  evaluate Jacobian
switch Jacobian_method
    case 1
        E=[];
        [J,h,func_evals_Jacobian]=jacobiansimple(objx,xGuess,option_Jacobian);
    case 2
        [J,h,func_evals_Jacobian,E]=jacobianlim(objx,xGuess,option_Jacobian);
    case 3
        [J,h,func_evals_Jacobian,E]=jacobianext(objx,xGuess,option_Jacobian);
end
%%  transform back to y space
df=D_unbnd2bnd(yGuess,lb,ub);
J=bsxfun(@times,J,df);
end

function    stop=eval_outputfun(OutputFcn,xGuess,optimValues,state)
stop=false;
if ~isempty(OutputFcn)
    if isa(OutputFcn,'function_handle')
        stop = OutputFcn(xGuess,optimValues,state);
    elseif isa(OutputFcn,'cell')
        stop=false(numel(OutputFcn),1);
        for i1=1:numel(OutputFcn)
            stop(i1)= OutputFcn{i1}(xGuess,optimValues,state);
        end
        stop=any(stop);
    end
end
end

function    OtptFcnVl=eval_PlotIterations(xGuess,OtptFcnVl,state)
xForm=OtptFcnVl.xForm;
switch state
    case 'init'
        fig=figure;
        %   optional figure name
        if ~isempty(OtptFcnVl.title)
            set(fig,'Name',OtptFcnVl.title)
        end
        subp(1)=subplot(3,2,1,'Parent',fig);title('current x');hold('on');
        bar(1:numel(xGuess),xGuess,'Parent',subp(1));
        
        subp(2)=subplot(3,2,2,'Parent',fig);title('current fval');hold('on');
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval,'Parent',subp(2));
        xlim(subp(2),[1 numel(OtptFcnVl.fval)])
        
        subp(3)=subplot(3,2,3,'Parent',fig);title('change from initial x');hold('on');
        subp(4)=subplot(3,2,4,'Parent',fig);title('change from initial fval');hold('on');
        subp(5)=subplot(3,2,5,'Parent',fig);title('attempted step');hold('on');
        
        if ~isempty(OtptFcnVl.pltfn)            
            xForm(:)=xGuess;
            subp(6)=subplot(3,2,6,'Parent',fig);title('user-supplied');hold('on');
            temp=OtptFcnVl.pltfn(xForm);
            plot(temp,'Parent',subp(6));
            xlim(subp(6),[1 size(temp,1)])
        else
            subp(6)=subplot(3,2,4,'Parent',fig);title('history residual norm');hold('on');
        end
        OtptFcnVl.fig=fig;
        OtptFcnVl.subp=subp;
        pause(.001)
    case 'iter'
        for i1=1:numel(OtptFcnVl.subp)
            delete(get(OtptFcnVl.subp(i1), 'Children'));
        end
        %   parameter values
        bar(1:numel(xGuess),xGuess,'Parent',OtptFcnVl.subp(1));
        xlim(OtptFcnVl.subp(1),[.5 numel(xGuess)+.5])
        %   residuals
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval,'Parent',OtptFcnVl.subp(2));
        xlim(OtptFcnVl.subp(2),[1 numel(OtptFcnVl.fval)])
        %   change in parameter values from start
        bar(1:numel(xGuess),xGuess-OtptFcnVl.xInit,'Parent',OtptFcnVl.subp(3));
        xlim(OtptFcnVl.subp(3),[.5 numel(xGuess)+.5])
        %   change in residuals from start
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval-OtptFcnVl.fvalInit,'Parent',OtptFcnVl.subp(4));
        xlim(OtptFcnVl.subp(4),[1 numel(OtptFcnVl.fval)])
        %   attempted step
        bar(1:numel(xGuess),OtptFcnVl.xLMstep(:),'Parent',OtptFcnVl.subp(5));
        xlim(OtptFcnVl.subp(5),[.5 numel(xGuess)+.5])
        %   user supplied function
        if ~isempty(OtptFcnVl.pltfn)            
            xForm(:)=xGuess;
            temp=OtptFcnVl.pltfn(xForm);
            plot(temp,'Parent',OtptFcnVl.subp(6));
            xlim(OtptFcnVl.subp(6),[1 size(temp,1)])
        else
            plot(OtptFcnVl.ResnormHist,'Parent',OtptFcnVl.subp(6));
        end
        OtptFcnVl.fig=OtptFcnVl.fig;
        OtptFcnVl.subp=OtptFcnVl.subp;
        pause(.001)
    case 'done'
        % Cleanup of plots, guis, or final plot
        close(OtptFcnVl.fig);
    otherwise
end
end

function    X=unbnd2bnd(Y,lb,ub)
%%  transforms variable from -infinity to infinity to domain [lb,ub]
%   complements bnd2unbnd
X=NaN(size(Y));
for i1=1:numel(Y)
    if isfinite(lb(i1)) && isfinite(ub(i1)) % lower and upper bound
        X(i1) = (lb(i1)+ub(i1))/2+ (ub(i1)-lb(i1))/2*sin(2*Y(i1)/(ub(i1)-lb(i1)));
        %         X(i1) = max(lb(i1),min(ub(i1),(sin(Y(i1))+1)/2.*(ub(i1) - lb(i1)) + lb(i1)));
        %         X(i1)=((atan(Y(i1))+pi/2)/pi).*(ub(i1)-lb(i1))+lb(i1);
    elseif isfinite(lb(i1)) && ~isfinite(ub(i1)) %  just lower bound
        X(i1)= lb(i1)-1 + sqrt(Y(i1).^2+1);
        %         X(i1)= lb(i1) + Y(i1).^2;
        %         X(i1)=exp(Y(i1))+lb(i1)-eps;
    elseif ~isfinite(lb(i1)) && isfinite(ub(i1)) %  just upper bound
        X(i1)= ub(i1)+1 - sqrt(Y(i1).^2+1);
        %         X(i1)=ub(i1)-Y(i1).^2 ;
        %         X(i1)=-exp(-Y(i1))+ub(i1)+eps;
    else %  no bounds
        X(i1)=Y(i1);
    end
end
end

function    Y=bnd2unbnd(X,lb,ub)
%%  transforms variable from [lb,ub] to domain -infinity to infinity
%   complements unbnd2bnd
Y=NaN(size(X));
for i1=1:numel(X)
    if isfinite(lb(i1)) && isfinite(ub(i1)) %   bounded on both ends
        Y(i1)=(ub(i1)-lb(i1))/2*asin((2*X(i1)-(ub(i1)+lb(i1)))/(ub(i1)-lb(i1)));
        %         Y(i1) = 2*pi+asin(max(-1,min(1,2*(X(i1) - lb(i1))/(ub(i1)-lb(i1)) - 1)));
        %         Y(i1)=tan((X(i1)-lb(i1))./(ub(i1)-lb(i1))*pi-pi/2);
    elseif isfinite(lb(i1)) && ~isfinite(ub(i1)) %  just lower bound
        Y(i1)=sqrt((lb(i1) -1 - X(i1)).^2-1);
        %         Y(i1) = sqrt(X(i1)-lb(i1));
        %         Y(i1)=log(X(i1)-lb(i1)+eps);
    elseif ~isfinite(lb(i1)) && isfinite(ub(i1)) %  just upper bound
        Y(i1)=sqrt((ub(i1) +1 - X(i1)).^2-1);
        %         Y(i1) = sqrt(ub(i1) - X(i1));
        %         Y(i1)=-log(-X(i1)+ub(i1)+eps);
    else %  no boundaries
        Y(i1)=X(i1);
    end
end
end

function    df=D_unbnd2bnd(Y,lb,ub)
%%  transform back to y space
df=NaN(1,numel(Y));
for i1=1:numel(Y)
    if isfinite(lb(i1)) && isfinite(ub(i1)) % lower and upper bound
        df(i1) = cos(2*Y(i1)/(ub(i1)-lb(i1)));
    elseif isfinite(lb(i1)) && ~isfinite(ub(i1)) %  just lower bound
        df(i1)= Y(i1)/sqrt(Y(i1).^2+1);
    elseif ~isfinite(lb(i1)) && isfinite(ub(i1)) %  just upper bound
        df(i1)= -Y(i1)/sqrt(Y(i1).^2+1);
    else %  no bounds
        df(i1)=1;
    end
end
end

