%*****************************************************************************************************************************************************************************************************************************
% Ground Resonance analysis solver
%
% 4 DoFs
% 2 longitudinal and lateral in-plane motion of the rotor hub
% 2 cyclic lag DoFs

%  [ b/2*Ibeta       0          0       (b/2)*Scsi  ]   [csi1c_dotdot]   [ Ibeta*Ccsistar       2*Ibeta*Omega       0           0   ]   [csi1c_dot]   [ Ibeta*(NUcsi2-1)      Ibeta*Ccsistar         0            0      ]   [csi1c]   [0]
%  | 0           b/2*Ibeta  -(b/2)*Scsi     0       | * |csi1s_dotdot| + | -2*Ibeta*Omega       Ibeta*Ccsistar      0           0   | * |csi1s_dot| + | -Ibeta*Ccsistar       Ibeta*(NUcsi2-1)       0            0      | * |csi1c] = |0|
%  | 0         -(b/2)*Scsi      Mx          0       |   |xh_dotdot   |   | 0                        0           Mx*Cxstar       0   |   |xh_dot   |   | 0                           0            Mx*omegax2       0      |   |xh   |   |0|   
%  [(b/2)*Scsi     0            0           My      ]   [yh_dotdot   ]   [ 0                        0               0     My*Cystar ]   [xh_dot   ]   [ 0                           0                0        My*omegay2 ]   [xh   ]   [0]
%
%*****************************************************************************************************************************************************************************************************************************

function [aut_max,P]=gr_eq(data,omega,mode)
Toll = 1e-5;
nItMax = 10;
% Adimensionalization factor used for the equations
adim = (0.5*data.I_b*data.n_b);

S_bad=data.S_b/data.I_b;
m_xad=(data.m_x+data.n_b*data.m_b)/adim;
m_yad=(data.m_y+data.n_b*data.m_b)/adim; 



M=eye(4);
M(4,1)=S_bad;
M(1,4)=S_bad;
M(3,2)=-S_bad;
M(2,3)=-S_bad;
M(3,3)=m_xad;
M(4,4)=m_yad;

sz = size(M,1);
p = zeros(2*sz,length(omega));
aut_max = zeros(length(omega),1);
EV = zeros(2*sz,2*sz,length(omega)); 

I=eye(sz);
Z=zeros(sz);

AA = [M Z; Z I];

AA_inv = inv(AA);

    for k=1:length(omega)   
      
      om=omega(k);
      k_xad=data.k_x/(adim*om^2);
      k_yad=data.k_y/(adim*om^2);
      c_xad=data.c_x/(adim*om);
      c_yad=data.c_y/(adim*om);
      % NUcsi2=(data.ecc*data.R*data.Scsi)/data.Ibeta+...
      %     data.Kcsi/(i^2*data.Ibeta);
      NUcsi2 = data.NUcsi_K2+data.NUcsi_K1/om^2;
      
      if mode == 0
        deutsch = 1;
      elseif mode == 1
        deutsch = 2*(1-cos(2*pi/data.n_b));
      elseif mode == 2
        deutsch = 2*(1-cos(2*2*pi/data.n_b));
      end

      c_xiad=deutsch*data.c_xi/(data.I_b*om);

      C=zeros(4);
      C(1,1)=c_xiad;
      C(2,1)=-2;
      C(1,2)=2;
      C(2,2)=c_xiad;
      C(3,3)=c_xad;
      C(4,4)=c_yad;
    
      K=zeros(4);
      K(1,1)=NUcsi2-1;
      K(2,1)=-c_xiad;
      K(1,2)=c_xiad;
      K(2,2)=NUcsi2-1;
      K(3,3)=k_xad;
      K(4,4)=k_yad;
    
      B = [C K; -I Z];
     
      A = -AA_inv*B;
      % The eigenvalues are obtained as lambda/Omega, beware!
      % That is why they are called p as the output, since they can be
      % considered as the n.d.s laplace variable
      p2(:,k) = eig(A);
      if k==1
      [EV(:,:,1), eg] = eig(A); 
      p(:,1) = diag(eg);
      else
        for j = 1:2*sz
        eg0 = p(j,k-1);
        ev0 = EV(:,j,k-1);
        nIt = 0;
            while (norm([eg0;ev0] - [p(j,k); EV(:,j,k)]) > Toll) && (nIt < nItMax) 
                A2 = [eye(2*sz)*ev0, eye(2*sz)*eg0 - A;
                     0,  2*ev0'];
                B = [-(eye(2*sz)*eg0 - A)*ev0;
                     1 - ev0'*ev0];
                X = A2\B;
                p(j,k) = eg0;
                EV(:,j,k) = ev0;
                eg0 = eg0 + X(1);
                ev0 = ev0 + X(2:end);
                nIt = nIt+1;
            end
        if nIt==nItMax
            p(j,k)=p2(j,k);
        else
        p(j,k) = eg0;
        end

        EV(:,j,k) = ev0;
        end
      end
      aut_max(k,1)=max(real(p(:,k)));
    end

P = p;


out = 1;
end
