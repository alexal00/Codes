%**********************************************************************
%                   Ground Resonance
%       Lag damping sizing using the Deutsch criterion
%
%
% Functions called: gr_eq.m, vg_plot
%**********************************************************************
function data=ae71(data,mode)
% Load data from data structure containing information about the rotor
Cx =data.c_x; Cy = data.c_y;
b = data.n_b; %m.r. n of blades  
omega = data.omega;

% limits of ground resonance analysis
OM1 = .4*omega; % 40% of nominal rotor omega
OM2 = 1.2*omega; % 120% of nominal rotor omega
OM_max = round(1.5*omega,-1); % max omega for plot results

%% Airframe natural frequencies
OMEGAx2 = data.k_x/data.m_x; 
OMEGAy2 = data.k_y/data.m_y; 
OMEGAx = sqrt(OMEGAx2); % natural frequency of longitudinal mode, dimensional [rad/s]
OMEGAy = sqrt(OMEGAy2); % natural frequency of lateral mode, dimensional [rad/s]


%% Define rotor lag dynamics parameters
S_b = data.S_b; % blade lag static moment [kg*m]
m_b = data.m_b; % blade mass [kg]
I_b = data.I_b; % blade lag moment of inertia [kg*m^2]

% non-dimensional rotating natural frequency of lag motion (with offset and spring) -> (ecc*R*Scsi)/Ibeta = 3*e/(2*(1-e))
% data.NUcsi2 = (data.ecc*data.R*data.Scsi)/data.Ibeta+data.Kcsi/(data.Ibeta*OMEGA^2);
% Southwell coefficients
data.NUcsi_K2 = (data.e*data.S_b)/data.I_b; % formula for K2 valid only for articulated rotor!!1
data.NUcsi_K1 = (data.k_xi/data.I_b);

data.NUcsi2 = data.NUcsi_K2+data.NUcsi_K1/(omega^2);
NUxi = sqrt(data.NUcsi2);

disp (['Non-dim. rotating natural frequency of lag motion (per-rev) = ' num2str(NUxi) ]);

if NUxi < 1
    disp ('Soft in-plane rotor (lag frequency < 1 per-rev)');
else
    disp ('Stiff in-plane rotor (lag frequency > 1 per-rev)');  
end

%% Deutsch criterion

if NUxi < 1
    if mode == 0
        deutsch = 1;
        damp = 'standard b2h';
    elseif mode == 1
        damp = 'inter-blade';
        deutsch = 2*(1-cos(2*pi/data.n_b));
    elseif mode == 2
        damp = 'inter-2-blade';
        deutsch = 2*(1-cos(2*2*pi/data.n_b));
    end
    % Lag damping necessary for support longitudinal mode
    Cxix = (b/4*(1-NUxi)/NUxi*data.S_b^2)/(data.c_x/OMEGAx2)*1/deutsch; %[N*m*s/rad]

    % Lag damping necessary for support lateral mode
    Cxiy = (b/4*(1-NUxi)/NUxi*data.S_b^2)/(data.c_y/OMEGAy2)*1/deutsch; %[N*m*s/rad]
    
    % select the higher damping between lateral and longitudinal results
    Cxi_min_req = max(Cxix,Cxiy);

    disp (['Deutsch criterion for ' damp]);
    disp (['Minimum lag damping required to stabilize landing gear longitudinal mode [N*m*s/rad] = ' num2str(Cxix) ]);
    disp (['Minimum lag damping required to stabilize landing gear lateral mode [N*m*s/rad] = ' num2str(Cxiy)]);
    disp (['Minimum lag damping required to avoid ground resonance [N*m*s/rad] = ' num2str(Cxi_min_req)]);
    disp (' ' );
    % Use the value given by Hammond (modified or not by the input)
    c_xi = data.c_xi;
    disp (['Lag damping used =' num2str(c_xi)])

else
    
    disp ('Deutsch criterion always satisfied (system always stable, no ground resonance), both lag and support damping not required');

end
    
%% Eigenvalues computing in case of soft in-plane rotor
% Vector of Omegas to investigate GR stability
Omega = [linspace(0.01,round(omega*1.5),400)]';
% Creater marker cell for plotting
% Marker = {'o','+','*','.','x','s','d','^'};
% For a soft-in-plane rotor plot the omega-f and omega-zeta diagrams of
% stbaility
if NUxi < 1
    % specified lag damping - Ccsi
    data.c_xi = c_xi;
    % Calculate eigenvalues and obtain maximum eigenvalue
    [p_max,p2]=gr_eq(data,Omega,mode);

    %% Plot the evoluction of the maximum eigenvalue with time
    % plot results
    % fighandle=figure;
    % fighandle.Name = ['Deutsch' int2str(mode)];
    % plot(Omega*60/(2*pi),p_max,'linewidth',1.5)
    % hold on
    % limit=axis;
    % plot([OM1 OM1]*60/(2*pi),[limit(3) limit(4)],'r--',[omega omega]*60/(2*pi),[limit(3) limit(4)],'k--',[OM2 OM2]*60/(2*pi),[limit(3) limit(4)],'r--')
    % xlabel('$\Omega$  [rpm]')
    % ylabel('Max(Re($\omega$))')
    % ylim([min(p_max) 0])
    % % title('Soft in-plane rotor, specified lag damping')
    % grid on
    % legend('solution','40\% 1/rev','1/rev','120\% 1/rev','Location','best')

    %% Omega-f and omega-zeta plots
    fighandle=figure;
    fighandle.Name = ['Vg_ccsimin' int2str(mode)];
    vg_plot(p2,Omega');
    % limit=axis;
    % plot([OM1 OM1]*60/(2*pi),[limit(3) limit(4)],'r--',[omega omega]*60/(2*pi),[limit(3) limit(4)],'k--',[OM2 OM2]*60/(2*pi),[limit(3) limit(4)],'r--')
    
    %% Root-locs of all eigenvalues with time
    % % Root locus
    % fighandle=figure;
    % fighandle.Name = ['RL_ccsi' int2str(mode)];
    % for kk=1:size(p2,1)
    % scatter(round(real(p2(kk,:)),3),imag(p2(kk,:)),'Marker',Marker{kk},'DisplayName',strcat('Eigenvalue ',num2str(kk))); hold on
    % end
    % legend
    % xlabel('Real'); ylabel('Imag');
    % ylim([-20 20])
    % grid on

end

%% Plot the coleman diagram
% %% Coleman diagram
% % To obtain the Coleman diagram the damping must be set to zero
% %%%%%%%%%% coupled case %%%%%%%%%% 
% data.c_xi = 0;
% data.c_x = 0; data.c_y = 0;
% [~,p4] = gr_eq(data,Omega,mode);
% 
% % Root locus
% % fighandle=figure;
% % fighandle.Name = ['RL_ccsi0' int2str(mode)];
% % for kk=1:size(p4,1)
% % scatter(round(real(p4(kk,:)),3),imag(p4(kk,:)),'Marker',Marker{kk},'DisplayName',strcat('Eigenvalue ',num2str(kk))); hold on
% % end
% % legend('show','Location','best')
% % xlabel('Real'); ylabel('Imag');
% % ylim([-20 20])
% % grid on
% 
% skip_flag = 0;
% aux = 0;
% aux2 = 0;
% fighandle = figure;
% fighandle.Name = ['Coleman' int2str(mode)];
% for i = 1:length(Omega)
% 
%     p = p4(:,i);
%     check = abs(real(p)) > 1e-3;
%     aux3 = 1; aux4 = 1;
%     if ~isempty(find(check, 1))
%         aux2 = aux2+1;
%         aux = aux+1;
%     else
%         aux = aux+1;
%     end
%     for j = 1:length(p)
%         if real(p(j)) > 1e-3 % roots in ground resonance instability range
%             p_skip(aux3,aux2) = p(j);
%             Omega_skip(i) = Omega(i);
%             skip_flag = 1;
%             aux3 = aux3+1;
%         else
%             p_ok(aux4,1) = p(j); % roots in stability range
%             Omega_ok(aux) = Omega(i);
%             aux4 = aux4+1;
% %             sol_col(j,i)=Omega(i)*abs(imag(p(j)));
%         end
%     end
%     plot(Omega(i)*60/(2*pi),Omega(i)*abs(imag(p_ok)),'b+'); hold on
%     clear("p_ok")
% end
% %%%%%%%%%% uncoupled case %%%%%%%%%%
% 
% data.S_b = 0;
% % The lag frequency os calculated without introducing this uncoupling
% [~,p_un] = gr_eq(data,Omega,mode);
% 
% % for j=1:size(p4,1)
% % plot(Omega_ok,Omega*abs(imag(p_ok(j,:))),'+'); hold on
% % end
% for j = 1:size(p_un,1) 
%     plot(Omega*60/(2*pi),Omega'.*abs(imag(p_un(j,:))),'kx','MarkerSize',2);hold on
% end
% limit=axis;
% 
% % plot on Coleman diagram the limits of instability range
% 
% if skip_flag == 1
% 
%     index_omega_skip = find(Omega_skip);
%     index_omega_skip_2 = find(Omega_skip==0);
% 
%     kk = 0;
%     for ii=1:length(index_omega_skip)-1
%          if index_omega_skip(ii+1) == index_omega_skip(ii)+1
%              kk=kk+1;
%              index_omega_skip_3(kk) = ii;
%              flag_middle(kk) = 0;
%          else
%              kk=kk+1;
%              index_omega_skip_4(kk) = ii;
%              flag_middle(kk) = 1;
%          end
%     end
% 
% 
%     hold on
%     plot(Omega_skip(index_omega_skip(1))*ones(1,2)*60/(2*pi),[limit(3) limit(4)],'r--') 
%     hold on
%     plot(Omega_skip(index_omega_skip(end))*ones(1,2)*60/(2*pi),[limit(3) limit(4)],'r--') 
% 
%     if find(flag_middle == 1) ~= 0
%         middle_index_omega_skip_1 = index_omega_skip(find(index_omega_skip_4));
%         middle_index_omega_skip_2 = index_omega_skip(find(index_omega_skip_4)+1);
%         hold on
%         plot(Omega_skip(middle_index_omega_skip_1)*ones(1,2)*60/(2*pi),[limit(3) limit(4)],'r--')
%         hold on
%         plot(Omega_skip(middle_index_omega_skip_2)*ones(1,2)*60/(2*pi),[limit(3) limit(4)],'r--')
%     end
% 
% end
% xlabel('$\Omega$  [rpm]')
% ylabel('$\omega$  [rad/s]')
% ylim([0 round(OM_max,0)])
% % title('Coleman diagram')
% 
% % Update the values for the output structure
data.S_b = S_b; data.c_xi = c_xi; data.c_x = Cx; data.Cy = Cy;
end




        
        

