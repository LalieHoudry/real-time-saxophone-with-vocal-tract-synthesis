%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%          Synthèse en TR de saxophone avec CV à 1 mode
%                + diagramme de bifurcation display
%
%    Lalie Houdry
%        16 juin 2026
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear functions
clear
close all hidden

%% Initialization of the audio output object

Fs = 8000 ; % Hz, sample frequency
Nbuf = 1024 ; % samples

ADW = audioDeviceWriter('SampleRate',Fs,'BufferSize',Nbuf) ; % ,'Driver','ASIO',"Device",'Focusrite USB ASIO'

%% Paramètres 

tvec = (0:Nbuf)/Fs ;                    
params.NmodesSax = 4 ; params.NmodesCV = 1 ;                  % nmodes
params.Nmodes = params.NmodesSax + params.NmodesCV ;
tCelsius = 20 ;                                               % temperature
c = 331.45 * sqrt((273.15+tCelsius)/273.15);
rho = 1.292*273.15/(273.15+tCelsius);
params.R1 = 8e-3 ; Zcs = rho*c/(pi*params.R1^2) ;             % R1
NsBD = 8*Nbuf ; 

params.epsi = 1e-3 ;                                          % param de régularisation
X0 = zeros(params.Nmodes+2,1) ; odeOptions = odeset('AbsTol',1e-2,'RelTol',1e-2);
params.om_r = 2*pi*2000 ; params.q_r = 0.5 ; params.dyn = 1 ; % params d'anche
fmin = 1 ; fmax = 4000 ;
fmin = 1 ; fmax = 4000 ; fvec = fmin:fmax ;
svec = (2*1i*pi).*fvec ;

%% Init
% Control parameters
gamma = 0.2 ;
zeta = 0.4 ;
% modal sax
load("ComplexModalParameters_18_YAS.mat","Cn","sn","notenames") ; 
params.Cks = Cn.' ; params.sks = sn.' ; 
fZmod = @(Zc,Ck,sk) Zc.*sum(Ck./(svec-sk) + conj(Ck)./(svec-conj(sk)),1);
Ck_sax = params.Cks(1:params.NmodesSax,2) ; sk_sax = params.sks(1:params.NmodesSax,2) ; 
Zmod_sax = fZmod(Zcs,Ck_sax,sk_sax) ;
% conduit vocal
Rm = 17.5e-3 ; params.Zms = params.R1^2/Rm^2 ;
w1 = 2*pi*200 ; xi1 = 0.2 ; A1 = 500 ; % params modaux du conduit vocal
fZmod_CV = @(Rm,An,wn,xin) rho*c/Rm^2.*sum(svec.*An./(svec.^2+wn.^2+2.*svec.*xin.*wn),1) ; 
Zmod_CV = fZmod_CV(Rm,A1,w1,xi1) ;
sk_CV = -2*xi1*w1 + w1*sqrt(1-xi1^2) ;
Ck_CV = (A1/2)*(1+xi1/sqrt(1-xi1^2)) ; 

%% Construction of the user interface

uiBuildControlFigure
set(groot,'defaulttextinterpreter','latex');

%% Execution of the synthesis loop

while 1

    if refreshBDb.Value
        if exist('BDax','var')
            delete(BDax);
            delete(BDbar);
        end
        BDax = uiaxes(controlfig,'Position',[margin/2 dl 6*dl 5*dl],'xlim',[0 2],'ylim',[1 1000],'zlim',[0 1],'CLim',[0 imag(sk_sax(end))./(2*pi)]) ; 
        BDax.XLabel.String = '$\gamma$'; BDax.XGrid = 'on'; % BDax.XTick = 0:0.1:gammalim(2); % BDax.XTickLabel(~~mod(10*BDax.XTick,5)) = {''};
        BDax.YLabel.String = '$f_{CV}$'; BDax.YGrid = 'on'; % BDax.YTick = 0:0.1:zetalim(2); % BDax.YTickLabel(~~mod(10*BDax.YTick,5)) = {''};
        BDax.Colormap = [0 0 0; hsv]; BDax.ZGrid = 'on';
        % BDax.CLim = [0.4*min(imag(params.sks)) 1.1*max(imag(params.sks))];
        bd = scatter3(BDax,NaN,NaN,NaN,20,NaN,'filled');
        BDbar = colorbar(BDax,'Ticks',imag(sk_sax)./(2*pi),'TickLabels',{'$s_{1,s}$','$s_{2,s}$','$s_{3,s}$','$s_{4,s}$'},"TickLabelInterpreter","latex") ;
        refreshBDb.Value=0; % Housekeeping

        sBD = zeros(NsBD,1);
        pgs = zeros(NsBD/Nbuf,1);
        pws = zeros(NsBD/Nbuf,1);
    end

    if playb.Value 
       [~,Xs] = ode45(@(t,Y) odefun_CV(Y,t,gamma,zeta,Ck_sax,sk_sax,Ck_CV,sk_CV,params), tvec, X0, odeOptions) ;
       X0 = Xs(end,:) ;
       ps = 2*sum(real(Xs(2:end,1:params.NmodesSax)),2) ; 
       s = ps*0.5 ; ADW([s(:) s(:)]) ;
       % playlamp.Color = min(1,std(ps))*[0 1 0]; % Light up in green depending of the sound level

       % à modifier 
       if drawBDb.Value 
            sBD = [sBD(Nbuf:end); s(:)];
            pgs = [pgs(2:end); gamma];
            pws = [pws(2:end); w1/(2*pi)];
            if (length(unique(pgs))==1)&&(length(unique(pws))==1)&&(~any((bd.XData==gamma)&(bd.YData==w1))) % Control parameters are static, but not on a point that was already drawn
                A = rms(sBD-mean(sBD));
                NTs = diff(find(diff(sign(sBD))>0));
                f0 = Fs/mean(NTs) * (A>1e-3);
                if isnan(f0), f0=0; end
                bd.XData = [bd.XData gamma];
                bd.YData = [bd.YData w1];
                bd.ZData = [bd.ZData A];
                bd.CData = [bd.CData f0];
            end
        end
    end

    if refreshimpb.Value
        Zmod_sax = fZmod(Zcs,Ck_sax,sk_sax) ; Zmod_CV = fZmod_CV(Rm,A1,w1,xi1) ;
        plot(zmodax,fvec(fmin:fmax),abs(Zmod_sax(fmin:fmax)/Zcs),fvec(fmin:fmax),abs(Zmod_CV(fmin:fmax)/Zcs),LineWidth=2) ; legend(zmodax,{"$Z_s$", "$Z_m$"}, "Interpreter", "Latex");
        Zmslbl.Text = ['$\frac{Z_{c,m}}{Z_{c,s}} = $ ' num2str(params.Zms)] ; 
        refreshimpb.Value = 0 ;
        % figure(22) ; clf ; 
        % plot(fvec,abs(Zmod_CV)/Zcs) ;
    end

    if regb.Value
        pks = Xs(:,1:params.NmodesSax) ;
        reg = which_reg(pks) ; 
        reglbl.Text = ['Register ' num2str(reg)] ; 
        regb.Value = 0 ; 
    end

    if pplotb.Value
        pm = 2*sum(real(Xs(2:end,params.NmodesSax+1:params.Nmodes)),2) ; Dp = ps - pm ;
        xdot = Xs(2:end,params.Nmodes+1) ; xdotdot = Xs(2:end,params.Nmodes+2) ; 
        x = (Dp - gamma) - params.q_r*xdot/params.om_r - xdotdot/(params.om_r^2) ; 
        u = zeta.*heaviside(x+1).*sign(gamma - Dp).*sqrt(abs(gamma - Dp)) ; % zeta/2.*(x+1+sqrt((x+1).^2+params.epsi)).*(gamma-Dp)./(((gamma-Dp).^2+epsi).^(1/4)) ; 

        figure(21) ; clf ; 
        subplot(311)
        plot(tvec(1:end-1),ps) ; legend("$\hat p_s$","Interpreter", "Latex") ;
        subplot(312)
        plot(tvec(1:end-1),u) ; legend("$\hat u$","Interpreter", "Latex") ; xlabel("temps (s)") ; 
        subplot(313)
        plot(tvec(1:end-1),x) ; legend("$\hat x$","Interpreter", "Latex") ; xlabel("temps (s)") ;

        pplotb.Value = 0 ; 
    end
    drawnow limitrate ; % Check the value of UI elements, but not too often (to avoid messing with audio stream)
end