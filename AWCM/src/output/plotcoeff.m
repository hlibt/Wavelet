clc, clf, clear all

%-----------------------------------------------%
% This file plots the detail (wavelet)          %
% coefficients on a dyadic grid. After the      %
% thresholding process, many detail coeff's     %
% will be removed from the grid, leaving a      %
% sparse dyadic grid structure.                 %
% INPUTS:                                       %
%           - 'coeffi.dat' type files           %
% OUTPUTS:                                      %
%           - plot of the dyadic grid           %
% AUTHOR:                                       %
%           Brandon Gusto,                      %
%           Department of Scientific Computing, %
%           Florida State University            %
%-----------------------------------------------%

% maximum level j to search for
Jmax=13;

% check if input files exist, if so load them
chkexist=zeros(Jmax);
for j=0:Jmax
    if exist(sprintf('_coeff_files/coeff%d.dat',j),'file')
        data{j+1}=load(sprintf('_coeff_files/coeff%d.dat',j)); 
        chkexist(j+1)=true;
    end
end

% plot results
figure(1);
axis([-1 1 0 Jmax+1]); hold on
for j=0:Jmax
    if chkexist(j+1)==true
        plot(data{j+1}(:,1),data{j+1}(:,2),'ko'); hold on
        %set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k');
    end
end
hold off

% set fontsize
fs=16;

% set specific y tick marks
set(gca,'fontsize',fs)
set(gca,'ygrid','on')
set(gca,'Ytick',0:Jmax+1);

% name the x and y labels
ylabel('$j$','Interpreter','LaTex','fontsize',fs)
xlabel('$x_{k}^{j}$','Interpreter','Latex','fontsize',fs)
