
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
n = 4;
beta = 40;
x_ini = [-2 -1 5 6]';
x_avg = (1/n) * (ones(1,n)*x_ini);
z_ini = (1-beta) * x_ini;
% z_ini = x_ini;
z_avg = (1/4) * (ones(1,4) * z_ini);
% d_ini = rand(n,1);
d_ini = zeros(n,1);
chi_ini = [ x_ini ; z_ini ; d_ini];
tspan = [1 1.0025];
options = odeset('RelTol',1e-7);
[t, chi] = ode45(@rslnt_iman, tspan, chi_ini, options);

%%%%%%%%%%%%%%%%%%%%%%% Another plot option %%%%%%%%%%%%%%%%%
% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 17;
opts.height     = 5;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure; clf

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;
axis tight
% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

% export to png
v = chi(:,1:4);
fig.PaperPositionMode   = 'auto';
plot(t, chi(:,1:4),'linewidth',1.5);
xlabel('Time (Seconds)');
ylabel('State infromation');
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position',[.45 .285 .25 .25])
box on 
indexOfInterest = (t < 1.0025) & (t > 1.0016); % range of t near consensus
plot(t(indexOfInterest),v(indexOfInterest,:),'linewidth',1);
ylabel('$x$','Interpreter','Latex');
grid on
axis tight
save2pdf('resilient_consensus_p_attack_hg',fig,600);
%%%%%%%%%%%%%%%%%%%%
function chi_dot = rslnt_iman(t, chi)
% beta1 = 40;
% beta2 = 1600;
% beta3 = 2 * sqrt(beta1 * beta2) + 20000000;

%%%%%%%%%%%%%%%%Gains used in ECC and EJC paper%%%%%%%%%%%
% beta1 = 50;
% beta2 = 100;
% beta3 = 100;
beta1 = 20;
beta2 = beta1^6;
beta3 = beta1^4;
% beta1 = 130;
% beta2 = 11000;
% beta3 = 4000;
%%% Unstable gains no attack
% beta1 = 2;
% beta2 = 2;
% beta3 = 4;
%%% stable gains no attack
% beta1 = 1;
% beta2 = 1;
% beta3 = 11;



% Ls = [1 -1 0 0;
%       -1 3 -1 -1;
%       0 -1 2 -1;
%       0 -1 -1 2];
Ls = [1 -1 0 0;
      0 1 -1 0;
      0 0 1 -1;
      -1 0 0 1];

  Fa = zeros(4,4);
  Fa(2,2) = -2;
  Fa(3,3) = -3;
  Ba = zeros(4,4);
  Ba(2,:) = -2 * [4 1 3 3];
%   Ba(3,:) = -1.5 * [3 2 3 4];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   x_dot = -  Ls * ((chi(1:4,:) - chi(9:12,:)) - (beta*chi(5:8,:))) ;
% alternative method   u = beta Ls z - beta^2 Ls x
%   x_dot = -beta0 * Ls *  (((chi(1:4,:))- chi(9:12,:) ) - (beta1*chi(5:8,:))) ; 
%%%%%%%%%%%%%%%%%%%%%%%%% Adding d to the whole systems
%   x_dot = - Ls *  (((chi(1:4,:)) ) - (beta1*chi(5:8,:))) + chi(9:12,:);
%   z_dot = - (beta2 * Ls *chi(1:4,:) ) - (beta3 * Ls * chi(5:8,:)) ;
%   d_dot = Fa * chi(9:12,:) + Ba * chi(1:4,:);
%   d_dot = Fa * chi(1:4,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%% Iqbal's virtual network%%%%%%%%%%%%%%%%%%
% beta1 = 2;
% beta2 = 2;
% beta3 = 2;
% beta2 = beta1^2;
% beta3 = beta1^2;
 x_dot = - Ls *  ((chi(1:4,:)) ) + (beta1 * Ls * (chi(5:8,:)) )  + Ls *  ((chi(9:12,:)) ); 
 z_dot = - beta2* ( Ls * chi(1:4,:)) - (beta3 * Ls * (chi(5:8,:))) ;
% z_dot = - beta2* ( Ls * (chi(1:4,:) )) - (beta3 * Ls * (chi(5:8,:) - chi(9:12,:))) ;
 d_dot = Fa * chi(9:12,:) + (Ba * chi(1:4,:));

%%%%%%%%%%%%%%%%%%%%%%%%%% Iqbal's virtual network not relative%%%%%%%%%%%%%%%%%%
% beta1 = 5;
% beta2 = beta1^6;
% beta3 = beta1^4;
%  x_dot = - Ls *  ((chi(1:4,:))- chi(9:12,:) ) + (beta1 * Ls * (chi(5:8,:)) )  + (Ls * chi(9:12,:)); 
%  z_dot = - ( Ls * ((beta2 * chi(1:4,:)) - chi(9:12,:) )) - (beta3 * Ls * (chi(5:8,:))) + (Ls * chi(9:12,:));
% % z_dot = - beta2* ( Ls * (chi(1:4,:) )) - (beta3 * Ls * (chi(5:8,:) - chi(9:12,:))) ;
%  d_dot = Fa * chi(9:12,:) + (Ba * chi(1:4,:));

%%%%%%%%%%%%%%%%%%%%%%%%%% Iman's virtual network%%%%%%%%%%%%%%%%%%
% beta0 = 1;
% beta1= 0.001;
% beta2 = beta1;
% beta3 =1;
%  x_dot = -beta0 * Ls *  ((chi(1:4,:))- chi(9:12,:) ) + (beta1 * Ls * (chi(5:8,:)) )  ; 
%  z_dot = - beta2* ( Ls * (chi(1:4,:) )) - (beta3 * Ls * (chi(5:8,:))) ;
%  d_dot = Fa * chi(9:12,:) + Ba * chi(1:4,:);
%%%%%%%%%%%%%%%%%%% Another way%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  x_dot = -k1 * Ls *  ((chi(1:4,:))- chi(9:12,:) ) + (Ls * chi(5:8,:)) ; 
%  z_dot = - ( Ls *chi(1:4,:) ) - (k2 * Ls * chi(5:8,:)) ;
%  d_dot = Fa * chi(9:12,:) + Ba * chi(1:4,:); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  chi_dot = [x_dot ; z_dot ; d_dot];
end