close all; clear all; clc
%% IB damping function
setPlot
fid = 1;
step = 1;
%% Boomerang configuration
% From Sela & Rosen, the boomeran configuration is defined as perfectly
% symmetric without chordwise offsets
% a = linspace(0,2,21);
% b = linspace(2,0,21);

h = [0.5 1 2];
options = optimset('Display','off');
R = 5.5;
e = 0.04*R;
blades = [3 4 5];
d = 0;
for Nb = blades
    dpsi = 2*pi/Nb;
    for ii = 1:Nb
        psi(ii,1) = (ii-1)*dpsi;
    end
    T = rotmat(Nb,psi);
    ln = 50;
    E_0 = zeros(ln,1);
    E_c = zeros(ln,1);
    C_std = T'*eye(Nb)*T;
    % C_0std = C_std(1,1);
    % C_cstd = C_std(2,2)*2;
    col = 'Collective, $h/e$=';
    cyc = 'Cyclic, $h/e$=';
    
    
    % tlo = tiledlayout('flow');
    o = linspace(0,1,ln);
    disp(['Number of blades:' int2str(Nb)])
    % for di = d
    fh = figure(fid); fid = fid+1;
    fh.Name = ['DampingEffibNb' int2str(Nb)];% 'd' num2str(di)];
    ax = gca;
    hold(ax,"on");
    % disp(['Chordiwse offset:' num2str(di)])
    for ii = 1:length(h)
        disp(['h/e:' num2str(h(ii))])
        cb = e.*h(ii).*(1+o);
        ca = e.*h(ii).*(1-o);
        % [~,~,C_Rstd] = ibfun(e,0,0,2*e*h(ii),2*e*h(ii),Nb);
        % C_Rstd = C_Rstd(1,1)*1/(4*(1-cos(dpsi)))*eye(Nb);
        % C_std = T'*C_Rstd*T;
        % C_0std = C_std(1,1);
        % C_cstd = C_std(2,2)*2;
        for jj = 1:length(o)
            % [Kxilad0,Kxilbd0,~] = ibfun(e,0,0,ca(jj),cb(jj),Nb);
            [Kxila,Kxilb,~] = ibfun(e,e*d,e*d,ca(jj),cb(jj),Nb);
            
            % if di~=0
            %     E_0(jj) = (Kxilb+Kxila)^2/(Kxilbd0^2);
            %     E_c(jj) = (Kxilb^2+Kxila^2+2*Kxilb*Kxila*cos(dpsi))/Kxilbd0^2;
            % else
            E_0(jj) = (Kxilb+Kxila)^2/Kxilb^2;
            E_c(jj) = (Kxilb^2+Kxila^2+2*Kxilb*Kxila*cos(dpsi))/Kxilb^2;
            % end
        end
        % [x,idx] = sort(x,1,"ascend");
        % C_0 = C_0(idx);
        % C_c = C_c(idx);

        plot(ax,o,E_0,'Color',Color{1},'LineStyle',LineStyle{ii},'DisplayName',[col num2str(h(ii))])
        plot(ax,o,E_c,'Color',Color{2},'LineStyle',LineStyle{ii},'DisplayName',[cyc num2str(h(ii))])
    end
    grid(ax,"on")
    % title(ax,['Chordwise offset d=' num2str(di)])
    xlabel(ax,'$o/(2h)$')
    % if di==0
    ylabel(ax,'$E_{fc}=C_{\xi}^{ib}/C_{\xi}^{std}$')
    % else
    %     ylabel(ax,'$E_{fc}=C_{\xi}^{ib}/C_{\xi}^{std}|_{d=0}$')
    % end
    legend(ax,"Location","best")

    % end

    % Only symmetrical configurations
    fh = figure(fid); fid = fid+1;
    fh.Name = ['DampingEffCycibNb' int2str(Nb)];% 'd' num2str(di)];
    ax = gca;
    hold(ax,"on");

    d_vec = linspace(0,2,ln);
    spanoff = [0.25 0.5 1 2];
    for kk = 1:length(spanoff)
        c_a = spanoff(kk);
        for jj = 1:length(d_vec)
        di = d_vec(jj);
        [Kxilad,Kxilbd,~] = ibfun(e,e*di,e*di,c_a,c_a,Nb);
        [Kxilar,Kxilbr,~] = ibfun(e,0,0,c_a,c_a,Nb);
        E_cd(jj) = (Kxilbd^2+Kxilad^2+2*Kxilbd*Kxilad*cos(dpsi))/(Kxilbr^2+Kxilar^2+2*Kxilbr*Kxilar*cos(dpsi));
        end
        cyc = 'Cyclic, $c_{a,b}/e$=';

        plot(ax,d_vec,E_cd,'Color',Color{kk},'LineStyle',LineStyle{ii},'DisplayName',[cyc num2str(c_a)])
    end
    grid(ax,"on")
    % title(ax,['Chordwise offset d=' num2str(di)])
    xlabel(ax,'$\bar{d}/e$')
    % if di==0
    ylabel(ax,'$E_{fc}=C_{\xi}^{ib}/C_{\xi}^{ib}|_{d=0}$')
    % else
    %     ylabel(ax,'$E_{fc}=C_{\xi}^{ib}/C_{\xi}^{std}|_{d=0}$')
    % end
    legend(ax,"Location","best")
end
%%
FolderName=strcat(pwd,'\images');
FoldeName=char(FolderName);
% FolderName = '\Images';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, fullfile(FolderName, [FigName, '.eps']),'epsc');    %<---- 'Brackets'
end
