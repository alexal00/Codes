function fid = plotting(results,data, fid,leg,name)
CTbet = results.cT;
CQbet = results.cQ;
Beta0 = results.Beta0;
% error_Tmom =zeros(length(theta0),1);
% error_Qmom =zeros(length(theta0),1);
error_Tbet =zeros(length(data.theta0),1);
error_Qbet =zeros(length(data.theta0),1);
lt = length(data.theta0);
for l=1:lt
    % error_Tmom(l)= abs(data.cT(l)-CTmom(l));
    % error_Qmom(l)= abs(data.cQ(l)-CQmom(l));
    error_Tbet(l)= abs(data.cT(l)-CTbet(l));
    error_Qbet(l)= abs(data.cQ(l)-CQbet(l));
    % error_Beta0(l)=abs(data.Beta0(l)-Beta0(l));
end

fig =figure(fid); fid = fid+1;
fig.Name = ['cT' name];
% errorbar(1:7,CTmom,error_Tmom); hold on
errorbar(1:lt,CTbet,error_Tbet','DisplayName',leg); hold on
plot(1:lt,data.cT)
grid on
legend('Interpreter','latex','Location','best')
xlabel('Experiment')
ylabel('$c_T$','Interpreter','latex')

fig = figure(fid); fid = fid+1;
fig.Name = ['cQ' name];
errorbar(1:lt,CQbet,error_Qbet','DisplayName',leg); hold on
plot(1:lt,data.cQ)
grid on
legend('Interpreter','latex','Location','best')
xlabel('Experiment')
ylabel('$c_Q$','Interpreter','latex')

% fig = figure(20+fid);
% fig.Name = ['Beta0' name];
% errorbar(1:7,Beta0,error_Beta0,'DisplayName',leg); hold on
% % plot(1:7,data.Beta0); hold on
% % plot(1:7,Beta02); hold on
% xlabel('Experiment')
% ylabel('$\beta_0$','Interpreter','latex')
% legend('Interpreter','latex','Location','best')
% grid on
end
