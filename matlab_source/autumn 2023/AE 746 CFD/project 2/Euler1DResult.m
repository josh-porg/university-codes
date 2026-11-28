classdef Euler1DResult
    %Euler1DResult Class to store the results of 1D euler simulation
    %   I am extremely tempeted to put the simulation in here (as I would
    %   in java) but I will resist the temptation as matlab does badly with
    %   objects. That is to say to run the simulatioin in the contrtucutor.

    properties
        order
        limiter % 0 for none, 1 for minmod, 2 for barth
        cfl % work out this
        numCells
        gamma
        flow % defines problem 1 or 2
        isGlobal % use global time stepping
        useSimpleBC % use simple boundary conditions?
        useFullEigs % use all three egien values (doen't play well with characteristic BC)
        x_points % x locations of mesh points

        Q_exact % exact solution
        Q_converged % Q at the converged state
        resHistroy % history of the density residual

        error
    end

    methods
        function obj = Euler1DResult(order, limiter, numCells, cfl, gamma, flow, isGlobal, useSimpleBC, useFullEigs, x_points, Q_exact, Q_converged, resHistroy)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.order = order;
            obj.limiter = limiter; % 0 for none, 1 for minmod, 2 for barth
            obj.cfl = cfl; % work out this
            obj.numCells = numCells;
            obj.gamma = gamma;
            obj.flow = flow; % defines problem 1 or 2
            obj.isGlobal = isGlobal; % use global time stepping
            obj.useSimpleBC = useSimpleBC; % use simple boundary conditions?
            obj.useFullEigs = useFullEigs; % use all three egien values (doen't play well with characteristic BC)
            obj.x_points = x_points; % x locations of mesh points
            
            % store exact solution
            obj.Q_exact = Q_exact;

            % store the results
            obj.Q_converged = Q_converged;
            obj.resHistroy = resHistroy;

            % compute error
            obj.error = computeError(obj);
        end

        function name = getName(obj)
            name = sprintf("order %g, limiter %g, numCells %g, cfl %g, flow %g, isGlobal %g, useSimpleBC %g, useFullEigs %g", obj.order, obj.limiter, obj.numCells, obj.cfl, obj.flow, obj.isGlobal, obj.useSimpleBC, obj.useFullEigs);
        end

        function error = computeError(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Q_exact = obj.Q_exact;
            Q_converged = obj.Q_converged;

            error = sqrt( sum((Q_converged - Q_exact).^2, 2) ); % L2 norm error
        end


        function plotThermoState(obj)
            Q_exact = obj.Q_exact;
            Q_converged = obj.Q_converged;
            x_points = obj.x_points;
            gamma = obj.gamma;

            figure;
            obj.plotTStat(x_points, Q_exact, gamma,'r-.');

            hold on;
            ax = get(gcf,"Children");
            for a = ax(:)
                hold(a,"on");
            end
            obj.plotTStat(x_points, Q_converged, gamma, 'b')

            hold("off")
        end


        function plotState(obj)
            x_points = obj.x_points;
            Q = obj.Q_converged;

            % Plot the conserved variables
            subplot(3,1,1);
            plot(x_points,Q(1,:),'r-.');
            xlabel('x');
            ylabel('rho');
            %ylim([0.4,1.1]);

            subplot(3,1,2);
            plot(x_points,Q(2,:),'r-.');
            xlabel('x');
            ylabel('rho*u');
            %ylim([-0.,0.8]);

            subplot(3,1,3);
            plot(x_points,Q(3,:),'r-.');
            xlabel('x');
            ylabel('E_t');
            %ylim([0.8,2.0]);
        end

        function plotTStat(obj, x_points, Q, gamma, lineStyle)
            [u, p, ~, rho, ~] = getFlowProperties(Q, gamma);

            % Plot the conserved variables
            subplot(3,1,1);
            plot(x_points, rho, lineStyle);
            xlabel('x');
            ylabel('rho');
            %ylim([0.4,1.1]);

            subplot(3,1,2);
            plot(x_points, p, lineStyle);
            xlabel('x');
            ylabel('p');
            %ylim([-0.,0.8]);

            subplot(3,1,3);
            plot(x_points, u, lineStyle);
            xlabel('x');
            ylabel('u');
            %ylim([0.8,2.0]);
        end


        function plotThermo(obj)
            x_points = obj.x_points;
            Q = obj.Q_converged;
            gamma = obj.gamma;


            [u, p, E, rho, c] = getFlowProperties(Q, gamma);

            % Plot the conserved variables
            subplot(5,1,1);
            plot(x_points,u,'-.');
            xlabel('x');
            ylabel('u');
            %ylim([0.4,1.1]);

            subplot(5,1,2);
            plot(x_points,p,'-.');
            xlabel('x');
            ylabel('p');
            ylim([-0.,0.8]);

            subplot(5,1,3);
            plot(x_points,rho,'-.');
            xlabel('x');
            ylabel('rho');
            ylim([0.8,2.0]);

            subplot(5,1,4);
            plot(x_points,c,'-.');
            xlabel('x');
            ylabel('c');
            ylim([0.8,2.0]);

            tau = p./rho;
            subplot(5,1,5);
            plot(x_points,tau,'-.');
            xlabel('x');
            ylabel('T (if R = 1)');
            ylim([0.8,2.0]);
        end

        function plotWithExactSoln(obj)
            x_points = obj.x_points;
            Q = obj.Q_converged;
            q_exact = obj.Q_exact;


            % Plot the conserved variables
            %rho
            subplot(3,1,1);
            plot(x_points,Q(1,:),'r-.');
            xlabel('x');
            ylabel('rho');

            hold on;
            subplot(3,1,1);
            plot(x_points,q_exact(1,:),'b-.');
            xlabel('x');
            ylabel('rho');
            %ylim([0.4,1.1]);
            hold off;

            %rho u

            subplot(3,1,2);
            plot(x_points,Q(2,:),'r-.');
            xlabel('x');
            ylabel('rho*u');

            hold on;
            subplot(3,1,2);
            plot(x_points,q_exact(2,:),'b-.');
            xlabel('x');
            ylabel('rho*u');
            %ylim([-0.,0.8]);
            hold off;

            % E
            subplot(3,1,3);
            plot(x_points,Q(3,:),'r-.');
            xlabel('x');
            ylabel('E_t');

            hold on;
            subplot(3,1,3);
            plot(x_points,q_exact(3,:),'b-.');
            xlabel('x');
            ylabel('E_t');
            %ylim([0.8,2.0]);
            hold off;
        end

    end
end