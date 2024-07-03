function [results,fid] = rotortrim(data,ff,fid)
    options = optimset('Display','off');

    fig = figure(fid); fid = fid+1;
    fig.Name = ['Rotorstate'];
    tlo = tiledlayout('flow');
    

    for kk = 1: length(ff.theta0)
        % RESULTS FROM THE TRIMMER
        data.mu = ff.mu(kk); 

        data.theta0 = ff.theta0(kk);
        
        data.mux = data.mu*cos(ff.alfaH(kk));
        data.muz = data.mu*sin(ff.alfaH(kk));
        data.v = ff.V(kk);
        data.alfaH = ff.alfaH(kk);
        data.alfaD = ff.alfaD(kk);
        data.theta1c = 0;
        data.theta1s = -ff.B1(kk);
        % data.Beta0 = ff.Beta0(kk)*pi/180;
        lambdai = ff.lambdai(kk);
        cT_0 = ff.cT(kk);

        tol = 1e-4;
        err = 2*tol;
        niter = 0 ;
        maxiter = 200;
        
        % Equation between BET and Momentum theory to obtain the induced velocity
        figure(200);
        ax(1) = subplot(3,1,1);
        xlabel(ax(1),'n')
        ylabel(ax(1),'error \%')
        ax(2) = subplot(3,1,2);
        xlabel(ax(2),'n')
        ylabel(ax(2),'$c_T$')
        ax(3) = subplot(3,1,3);
        xlabel(ax(3),'n')
        ylabel(ax(3),'$\lambda_i$')
        h(1) = animatedline(ax(1),'Marker','o');
        h(2) = animatedline(ax(2),'Marker','o');
        h(3) = animatedline(ax(3),'Marker','o');
        while err>tol && niter<maxiter
            % Option 2:
            x0 = lambdai;
            lambda = fsolve(@(x) infloweq(x,data,cT_0),x0,options);
            lambdai = lambda-data.mu*sin(data.alfaD);
            % Increase inflow by a value 1/B to acount for tiplosses
            % lambdai = lambdai./data.B;
            
            forces=BET(data,lambdai);
            cT_1 = forces.CT;
            err = abs(cT_1-cT_0);
            
            addpoints(h(1),niter,err)
            addpoints(h(2),niter,cT_1)
            addpoints(h(3),niter,lambdai)
            drawnow limitrate
            % cT_0 = cT_0+(cT_1-cT_0)*0.1;
            cT_0 = cT_1;
            niter = niter+1;
        end
        
        if niter==maxiter
            count = count+1;
            err_vec(count,:) = [err,kk];
        end
        
        Tbet(kk) = forces.T; 
        Qbet(kk) = forces.Q;
        Beta0(kk) = forces.Beta0;
        Beta02(kk) = forces.Beta02;
        CTbet(kk) = forces.CT;
        CQbet(kk) = forces.CQ;
        
        info.mux = data.mux;
        info.muz = data.muz;
        info.R = data.R;
        info.Omega =data.Omega;
        info.y0 = data.e;
        info.Beta0 = forces.Beta0*pi/180;
        info.Beta1c = forces.a1;
        rotorstate(info,kk,'ff',tlo);
    end
    results.Tbet = Tbet;
    results.Qbet = Qbet;
    results.Beta0 = Beta0;
    results.cT = CTbet;
    results.cQ = CQbet;
end