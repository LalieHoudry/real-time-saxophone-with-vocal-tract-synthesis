%% Build the user interface

% size 
dtiers = 500 ; % Width of the control area 
%heightui = 300 ; % Height of the control area 
margin = 40 ; dl = (dtiers-margin)/6 ; 

controlfig = uifigure ; controlfig.Name = "Real time synth saxophone with vocal tract - NL dynamic display" ; 
controlfig.Position = [10 200 3*dtiers dtiers] ;

%%
% buttons
playb = uibutton(controlfig,"state","Text","Play", "Icon","play.png", "IconAlignment","top", "Position",[dtiers+5*dl (margin+9*dl)/2 dl dl]);
% stopb = uibutton(controlfig,"state","Text","Stop synth", "Icon","stop.png", "IconAlignment","top", "Position",[dtiers+margin/2 margin/2 dl dl]);
% playlamp = uilamp(controlfig,"Position",[dtiers+margin/2 + 4*dl/3 margin/2+dl/3 dl/3 dl/3],"Color",[0 0 0]) ; 
regb = uibutton(controlfig,"state","Text","Which register ?","Interpreter",'Latex',"Position", [dtiers+margin/2+2*dl margin 2*dl margin]);
reglbl = uilabel(controlfig,"Position",[dtiers+margin+4*dl margin 2*dl margin],'Text', 'Register ',"Interpreter","Latex");
pplotb = uibutton(controlfig,"state","Text","Draw waveform","Interpreter",'Latex',"Position", [dtiers margin 2*dl margin]);

% impedance display
zmodax = uiaxes(controlfig) ; zmodax.Position = [2*dtiers+margin/2 margin/2+dl 6*dl 5*dl] ;
plot(zmodax,fvec(fmin:fmax),abs(Zmod_sax(fmin:fmax)/Zcs),fvec(fmin:fmax),abs(Zmod_CV(fmin:fmax)/Zcs)) ; legend(zmodax,{"$Z_s$", "$Z_m$"}, "Interpreter", "Latex");
title(zmodax,"Imp\'edances d'entr\'ee (approx modale)", "Interpreter", "Latex") ; grid(zmodax,"on") ;
xlabel(zmodax,"fr\'equence (Hz)","Interpreter","Latex") ; ylabel(zmodax,"$|Z_{mod}|$","Interpreter","Latex") ; % zmodax.XLabel.String = "fr\'equence (Hz)" ; zmodax.YLabel.String = "dB" ; 
refreshimpb = uibutton(controlfig,"state","Text","Refresh impedance plot","Interpreter",'Latex',"Position", [2*dtiers+margin/2+4*dl margin 2*dl margin]);
Zmslbl = uilabel(controlfig,"Position",[2*dtiers+dl margin 2*dl margin],'Text', ['$\frac{Z_{c,m}}{Z_{c,s}} = $ ' num2str(params.Zms)],"Interpreter","Latex") ;

% sliders
gammaslider = uislider(controlfig,"Limits",[0 2],"Position",[dtiers+margin/2 (margin+11*dl)/2 4*dl dl],"Value",gamma,'MajorTicks',[0 1/3 1/2 1 1.5 2],'MajorTickLabels',{'0' '1/3' '0.5' '1' '1.5' '2'},'MinorTicks',0:0.1:1.4,'ValueChangingFcn',@(gammaslider,event) assignvar(event,gammaslider,'gamma')) ;
uilabel(controlfig,"Position",[dtiers+margin/2+4*dl+15 (margin+11*dl)/2-7 margin/2 margin/2],'Text','$\gamma$',"Interpreter","Latex");
zetaslider = uislider(controlfig,"Limits",[0 2],"Position",[dtiers+margin/2 (margin+9*dl)/2 4*dl dl],"Value",zeta,'MajorTicks',0:0.5:2,'MinorTicks',0:0.1:2,'ValueChangingFcn',@(zetaslider,event) assignvar(event,zetaslider,'zeta')) ;
uilabel(controlfig,"Position",[dtiers+margin/2+4*dl+15 (margin+9*dl)/2-7 margin/2 margin/2],'Text','$\zeta$',"Interpreter","Latex");

f1slider = uislider(controlfig,"Limits",[1 1000],"Position",[dtiers+margin/2 (margin+7*dl)/2 2*dl dl],"Value",w1/(2*pi),'ValueChangingFcn',@(f1slider,event) changingf1(event,f1slider,xi1,A1)); 
uilabel(controlfig,"Position",[dtiers+margin/2+2*dl+15 (margin+7*dl)/2-7 margin/2 margin/2],'Text','$f_{CV}$',"Interpreter","Latex");
xi1slider = uislider(controlfig,"Limits",[0.001 1.5],"Position",[dtiers+margin/2+3*dl (margin+7*dl)/2 2*dl dl],"Value",xi1,'ValueChangingFcn',@(xi1slider,event) changingxi1(event,xi1slider,w1,A1)); 
uilabel(controlfig,"Position",[dtiers+margin/2+5*dl+15 (margin+7*dl)/2-7 margin margin/2],'Text','$\xi_{CV}$',"Interpreter","Latex");
A1slider = uislider(controlfig,"Limits",[1 2000],"Position",[dtiers+margin/2 (margin+5*dl)/2 2*dl dl],"Value",A1,'ValueChangingFcn',@(A1slider,event) changingA1(A1slider,event,w1,xi1)); 
uilabel(controlfig,"Position",[dtiers+margin/2+2*dl+15 (margin+5*dl)/2-7 dl margin/2],'Text','$A_{CV}$',"Interpreter","Latex");
Rmslider = uislider(controlfig,"Limits",[1 30],"Position",[dtiers+margin/2+3*dl (margin+5*dl)/2 2*dl dl],"Value",Rm*10^3,'ValueChangingFcn',@(Rmslider,event) changingRm(Rmslider,event,Zmslbl,params.R1)); 
uilabel(controlfig,"Position",[dtiers+margin/2+5*dl+15 (margin+5*dl)/2-7 dl margin/2],'Text','$R_{CV}$ (mm)',"Interpreter","Latex");

f0slider = uislider(controlfig,'Limits',[imag(sks(1,1))/(2*pi) imag(sks(1,end))/(2*pi)],"Position", [dtiers+margin/2 (margin+3*dl)/2 5*dl dl],"Value",imag(sk_sax(1))/(2*pi),'ValueChangedFcn',@(f0slider,event) changedf0(event,f0slider,params)) ; 
uilabel(controlfig,"Position",[dtiers+margin/2+5*dl+25 (margin+3*dl)/2-7 margin/2 margin/2],'Text','$f_0$',"Interpreter","Latex");

% bifurcation diagram
refreshBDb = uibutton(controlfig,"state","Text","Reset bifurcation diagram","Interpreter",'Latex',"Position",[margin/2 margin/2 2.5*dl margin]);
drawBDb = uibutton(controlfig,"state","Text","Draw bifurcation diagram","Interpreter",'Latex',"Position",[3*dl+margin/2 margin/2 2.5*dl margin]);

drawnow;
%% Fonctions utiles (T. Colinot : main_Modal_Clarinet_ATIAM_AM_BD.m)

% gamma et zeta
function assignvar(event,sld,varstr) 
assignin('base',varstr,event.Value);
end

% fréquence / longueur du résonateur + Ck_sax, sk_sax
function changedf0(event,sld,parameters)
fk = event.Value ;
[sk_s idx] = min(abs(imag(parameters.sks(1,:)/(2*pi))-fk)) ;
for im = 1:parameters.NmodesSax
    evalin('base',['Ck_sax(' num2str(im) ')=' num2str(parameters.Cks(im,idx)) ';']) ; 
    evalin('base',['sk_sax(' num2str(im) ')=' num2str(parameters.sks(im,idx)) ';']) ;
end
end

% Ck_CV et sk_CV
function changingf1(event,sld,xi1,A1)
w1 = event.Value*2*pi ; 
sk = -xi1*w1 + 1i*w1*sqrt(1-xi1^2) ;
assignin('base','w1',w1);
evalin('base',['sk_CV =' num2str(sk) ';']); 
end
function changingxi1(event,sld,w1,A1)
xi1 = event.Value ; 
sk = -xi1*w1 + 1i*w1*sqrt(1-xi1^2) ;
Ck = (A1/2)*(1+1i*w1*xi1/(w1*sqrt(1-xi1^2))) ;
assignin('base','xi1',xi1);
evalin('base',['sk_CV =' num2str(sk) ';']); 
evalin('base',['Ck_CV =' num2str(Ck) ';']); 
end
function changingA1(event,sld,w1,xi1)
A1 = event.Value ; 
Ck = (A1/2)*(1+1i*w1*xi1/(w1*sqrt(1-xi1^2))) ;
assignin('base','A1',A1);
evalin('base',['Ck_CV =' num2str(Ck) ';']); 
end
function changingRm(event,sld,uilbl,R1)
Rm = event.Value*10^(-3) ;
Zms = (R1^2)/(Rm^2) ;
assignin('base','Rm',Rm); 
evalin('base',['params.Zms =' num2str(Zms) ';']); 
end

%%

% % pressure display
% psim = uiaxes(controlfig) ; psim.Position = [margin/2 2*dl 5*dl 3*dl] ;
% psim.XGrid = "on" ; xlabel(psim,"temps (s)","Interpreter","Latex") ;
% plot(psim,tvec,zeros(length(tvec),1)) ;
% legend(psim,"$\hat p_s$","Interpreter", "Latex") ;