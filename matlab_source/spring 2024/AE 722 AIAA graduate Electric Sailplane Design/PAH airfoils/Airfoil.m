classdef Airfoil
    %Airfoil An Object containing data from an airfoil file
    %   All units radians, 1\rad, or dimentionless

    properties
        name
        dataSource
        c_l_alpha
        alpha_0
        c_l_0
        alpha_star
        alpha_prime
        c_l_star
        alpha_c_l_max
        c_l_max
        alpha_post_stall
        c_l_post_stall
    end

    methods
        function obj = Airfoil(name, dataSource, c_l_alpha, alpha_0, c_l_0, alpha_star,alpha_prime, c_l_star, alpha_c_l_max, c_l_max, alpha_post_stall, c_l_post_stall)
            % Constructor for the Airfoil class
            % Initialize properties based on input arguments

            obj.name = name;
            obj.dataSource = dataSource;
            obj.c_l_alpha = c_l_alpha;
            obj.alpha_0 = alpha_0;
            obj.c_l_0 = c_l_0;
            obj.alpha_star = alpha_star;
            obj.alpha_prime = alpha_prime;
            obj.c_l_star = c_l_star;
            obj.alpha_c_l_max = alpha_c_l_max;
            obj.c_l_max = c_l_max;
            obj.alpha_post_stall = alpha_post_stall;
            obj.c_l_post_stall = c_l_post_stall;
        end

        function createAirfoilFile(obj, filename)
            % Create a text file with airfoil data
            
            % Extract properties from the Airfoil object
            name = obj.name;
            dataSrc = 'Comprehensive Reference Guide to Airfoil Sections for Light Aircraft';
            clalpha = obj.c_l_alpha;
            alpha_o = obj.alpha_0;
            c_l_o = obj.c_l_0;
            alpha_star = obj.alpha_star;
            alpha_prime = obj.alpha_prime;
            c_l_star = obj.c_l_star;
            alpha_c_l_max = obj.alpha_c_l_max;
            c_l_max = obj.c_l_max;
            alpha_post_S = obj.alpha_post_stall;
            c_l_post_S = obj.c_l_post_stall;
            
            % Create the formatted content
            content = sprintf('[Airfoil]\nName=%s\nData Source=%s\nclalpha [1/rad]=%.14f\nalpha_o [deg]=%.14f\nc_l_o [-]=%.14f\nalpha" [deg]=%.14f\nalpha* [deg]=%.14f\nc_l* [-]=%.14f\nalpha_c_l_max [deg]=%.14f\nc_l_max [-]=%.14f\nalpha_post_S [deg]=%.14f\nc_l_post_S [-]=%.14f\n', ...
                name, dataSrc, clalpha, rad2deg(alpha_o), c_l_o, rad2deg(alpha_prime), rad2deg(alpha_star), c_l_star, rad2deg(alpha_c_l_max), c_l_max, rad2deg(alpha_post_S), c_l_post_S);
            
            % Write content to the specified file
            fid = fopen(filename, 'w');
            fprintf(fid, content);
            fclose(fid);
        end

        % getPlotPoints -Returns the points to plot in degrees
        function [alpha, c_l] = getPlotPoints(obj)
            alpha = [obj.alpha_0, obj.alpha_star, obj.alpha_c_l_max, obj.alpha_post_stall];
            c_l = [0, obj.c_l_star, obj.c_l_max, obj.c_l_post_stall];
        end

        % get_C_l_alpha(alpha, airfoil) - Return C_l_alpha of the airfoil
        % for a given alpha
        function C_l_alpha = get_C_l_alpha(airfoil, alpha)
            % C_l_alpha_linear = (airfoil.c_l_star - 0) ./ (airfoil.alpha_star - airfoil.alpha_0);
            % C_l_alpha_nonlinear = (airfoil.c_l_max - airfoil.c_l_star) ./ (airfoil.alpha_c_l_max - airfoil.alpha_star);
            % C_l_alpha_postStall = (airfoil.c_l_post_stall - airfoil.c_l_max) ./ (airfoil.alpha_post_stall - airfoil.alpha_c_l_max);

            % produce C_l_alpha as a peicwise constant funtion
            % C_l_alpha = heaviside(-(alpha - airfoil.alpha_star)) .* C_l_alpha_linear + ...
            %     heaviside((alpha - airfoil.alpha_star)) .* heaviside(-(alpha - airfoil.alpha_c_l_max)) .* C_l_alpha_nonlinear + ...
            %     heaviside((alpha - airfoil.alpha_c_l_max)) .* heaviside(-(alpha - airfoil.alpha_post_stall)) .* C_l_alpha_postStall;

            c_l_spl = airfoil.CL_alpha_Smooother(0, 2);
            C_l_alpha = ppval(fnder(c_l_spl),alpha);
        end

         % get_C_l(airfoil, alpha) - Return C_l of the airfoil for a given
         % alpha
%         function C_l = get_C_l(airfoil, alpha)
%             % note there is a fractional offset of the 0 lift after the
%             % stall to allow the value post stall be be read acurately as
%             % nonzero
%             C_l = (airfoil.c_l_0 + integral(@(a) airfoil.get_C_l_alpha(a), 0, alpha)) * heaviside(-(alpha - (airfoil.alpha_post_stall+.000001)));
%         end

        function C_l = get_C_l(airfoil, alpha)
            % produce the smoothed cl's
            spl = airfoil.CL_alpha_Smooother(0, 2);
            C_l = fnval(spl,alpha);
            %C_l = (airfoil.c_l_0 + integral(@(a) airfoil.get_C_l_alpha(a), 0, alpha)) * heaviside(-(alpha -airfoil.alpha_post_stall));

            % Clamp to zero post stall
            C_l = C_l .* heaviside(-(alpha - airfoil.alpha_post_stall));
        end

        % this should be private
        % Note: I probably dont need to input four point matrix because I
        % can get it from airfoil
        function spl = CL_alpha_Smooother(airfoil, plot_yn, n)
            %   Example:
            %       points = [-.174533, 0; .1396, 2.1; .2443, 2.4; .2793, 2.2];
            %       spl = generate_spline(points, true);
            %       n between 1 and 5, 5 being most linear in linear region (adds spline points on line to force the function)
            [alphas, c_ls] = airfoil.getPlotPoints();
            four_points_matrix = [alphas', c_ls'];


            x_data = four_points_matrix(:,1)';
            y_data = four_points_matrix(:,2)';

            slope_at_point_one = (y_data(2) - y_data(1)) / (x_data(2) - x_data(1));

            %%%%%%%%%%%%%%%%%%%%%%%  Adjust linearity of linear region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the step size for inserting the new points
            step_size = 1 / (n + 1);

            % Initialize arrays to store the new points
            new_x = zeros(1, n);
            new_y = zeros(1, n);

            % Calculate the coordinates of the new points
            for k = 1:n
                new_x(k) = (1 - k * step_size) * x_data(1) + k * step_size * x_data(2);
                new_y(k) = (1 - k * step_size) * y_data(1) + k * step_size * y_data(2);
            end

            % Insert new points into x_data and y_data
            x_data = [x_data(1), new_x, x_data(2:end)];
            y_data = [y_data(1), new_y, y_data(2:end)];

            %%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            err = 1; m_end = 1; i = 0; signn = 1; d_maximum = 1;
            while abs(err) > 10^-3 % changes end slope until point 3 slope is zero
                i = i + 1;
                err_old = d_maximum;

                slope_at_point_four = (y_data(4) - m_end * y_data(3)) / (x_data(4) - x_data(3)); % KEEP SCALING FACTOR m_end
                spl = csape(x_data, y_data, 'complete', [slope_at_point_one, slope_at_point_four]);
                d_maximum = ppval(fnder(spl), x_data(end-1));


                if d_maximum > 0
                    m_end = m_end + signn*m_end * err * 10^-2;
                else
                    m_end = m_end - signn*m_end * err * 10^-2;
                end

                err = d_maximum;

                if i > 500
                    disp('error, no spline convergence')
                    break
                end

                if mod(i, 2) == 0
                    if abs(err_old) < abs(err)
                        signn = signn*-1;
                    end
                end



            end

            if nargin < 2
                plot_yn = true; % Default to plotting
            end

            if plot_yn
                x_interp = linspace(min(x_data), max(x_data), 1000);
                y_interp = fnval(spl, x_interp);

                % Plot the original data and the spline
                plot(x_data([1,end-2,end-1,end]), y_data([1,end-2,end-1,end]), 'ko', 'MarkerSize', 8);
                hold on;
                plot(x_data(2:end-3), y_data(2:end-3), 'r*', 'MarkerSize', 8);
                plot(x_interp, y_interp, '-', 'LineWidth', 1.5);        xlabel('alpha');
                ylabel('C_L');
            end
        end
    end
end