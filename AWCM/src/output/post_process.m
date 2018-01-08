function post_process(numsteps,plot_soln,plot_coeffs,freq)

    %-----------------------------------------------%
    % This file reads the output data from awcm.cpp %
    % The data is the field/solution from the pde   %
    % solver as well as the magnitude of the        %
    % wavelet coefficients at each timestep.        %
    % The detail (wavelet) coefficients are plotted %
    % on a dyadic grid.                             %
    %                                               %
    % INPUTS:                                       %
    %           - number of expected timesteps      %
    %           - 'ui.dat' type files               %
    %           - 'coeffi.dat' type files           %
    % OUTPUTS:                                      %
    %           - plot of the solution in time      %
    %           - plot of the dyadic grid           %
    % AUTHOR:                                       %
    %           Brandon Gusto,                      %
    %           Department of Scientific Computing, %
    %           Florida State University            %
    %-----------------------------------------------%

    % set fontsize for all plots
    fs=16;

    if plot_coeffs == true

        % check if input files exist, if so load them
        chkexist = zeros(numsteps);
        for j=1:numsteps
            if exist( sprintf('_coeff_files/coeff%d.dat',j), 'file' )
                data{j} = load( sprintf('_coeff_files/coeff%d.dat', j) ); 
                Jmax = max( data{j}(:,2) );
                chkexist(j) = true;
            end
        end

        % plot the detail coefficients
        figure(1);
        for i=1:numsteps
            if chkexist(i) == true
                h = plot( data{i}(:,1), data{i}(:,2), 'ko' );
                set(h,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            end
            axis([-1 1 0 Jmax+1]);
            hold off
            pause(0.00001);
        end

    end

    % solution to burger's equation
    %vxn = 1000;
    %vx = linspace(-1,1,vxn);
    %vtn = n;
    %vt = linspace(0,1,vtn);
    %vu = burgers_solution(0.01,vxn,vx,vtn,vt);

    % input data and make movie
    figure(2)
    for i=0:numsteps-1
        if ( mod(i,freq)==0 )
            filename = sprintf('_soln_files/u%d.dat',i);
            U = load(filename);
            [x,I] = sort( U(:,1) );
            u = U(:,2);
            plot( x , u(I), 'b' ); grid on; hold on;
    %        plot( vx, vu(:,i+1) ,'r'); hold on;
            plot( x , u(I), '.b' );
            hold off;
            axis( [ -1 1 0 1.5 ] );
            pause(0.0000001);
        end
    end

end
