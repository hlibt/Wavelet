function convergence_analysis(minPts,numRuns)

    % method names
    methds = { 'upwind_difference', 
                'central_difference',
                'linear_upwind_difference',
                'muscl2',
                'piecewise_parabolic',
                'mono_piecewise_parabolic' };

    % set latex parameter for plotting
    set(0,'defaulttextinterpreter','latex');

    % vectors
    p1 = zeros(numRuns,1);
    p2 = zeros(numRuns,1);
    pinf = zeros(numRuns,1);
  
    % loop over all methods
    for i=1:length(methds)

        % check if output folder exists
        if exist( sprintf('%s', methds{i}), 'dir' ) == 7

            % number of points
            numPts = minPts;

            % loop through spatial convergence runs
            for j=1:numRuns

                % find last file
                k = 0;
                while (1)
                    
                    % name of norm file
                    fileName = sprintf('%s/_run_%dpts/soln%d.dat',methds{i},numPts,k);
                    
                    % check if exists
                    if exist( fileName, 'file' ) == 2
                        k = k + 1;
                    else
                        k = k - 1;
                        break;
                    end

                end

                % load data
                fileName = sprintf('%s/_run_%dpts/soln%d.dat',methds{i},numPts,k);
                Data = load( fileName );

                % record norms at this final time
                L1(j) = Data(1,3);
                L2(j) = Data(1,4);
                Linf(j) = Data(1,5);

                % increment
                nvec(j) = numPts; 
                numPts = numPts * 2;

                % compute the order of convergence
                if j>1
                    p1(j) = log( L1(j-1) / L1(j) ) / log(2);
                    p2(j) = log( L2(j-1) / L2(j) ) / log(2);
                    pinf(j) = log( Linf(j-1) / Linf(j) ) / log(2);
                end

            end

            % make table
            disp(methds{i});
            T = table( nvec', L1', L2', Linf', p1, p2, pinf )

        end

    end % methods loop

end
