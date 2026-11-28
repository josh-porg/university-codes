classdef Untitled < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        LeftPanel                     matlab.ui.container.Panel
        BeginButton                   matlab.ui.control.Button
        Label                         matlab.ui.control.Label
        CenterPanel                   matlab.ui.container.Panel
        EnemyLaunchAngle45Label       matlab.ui.control.Label
        EnemyVelocity1735msLabel      matlab.ui.control.Label
        TextArea_2                    matlab.ui.control.TextArea
        UIAxes                        matlab.ui.control.UIAxes
        RightPanel                    matlab.ui.container.Panel
        Shot1AngleLabel               matlab.ui.control.Label
        Angle1Input                   matlab.ui.control.NumericEditField
        Shot2AngleLabel               matlab.ui.control.Label
        Angle2Input                   matlab.ui.control.NumericEditField
        Shot3AngleLabel               matlab.ui.control.Label
        Angle3Input                   matlab.ui.control.NumericEditField
        Timesincelaunch0secondsLabel  matlab.ui.control.Label
        Shot1Velocity0msLabel         matlab.ui.control.Label
        Shot2Velocity0msLabel         matlab.ui.control.Label
        Shot3Velocity0msLabel         matlab.ui.control.Label
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    
    methods (Access = private)
        
        function val = collide(app,xm,ym,x1i,y1i,x2i,y2i,x3i,y3i,i)
            
            dist1=sqrt((xm-x1i)^2+(ym-y1i)^2);
            dist2=sqrt((xm-x2i)^2+(ym-y2i)^2);
            dist3=sqrt((xm-x3i)^2+(ym-y3i)^2);
            
            if (i==1)
                if(dist1<=5 || dist2<=5 || dist3<=5)
                    val=1;
                else
                    val=0;
                end
            else
                
                 val=[dist1,dist2,dist3];
            
            end
            
        end
        
        function [x1i,y1i,x2i,y2i,x3i,y3i] = launch(app,lt,vii,a1,a2,a3,t)
            
            x1i=307E3-(t-lt)*vii*cosd(a1);
            y1i=(t-lt)*vii*sind(a1)-.5*9.8*(t-lt)^2;
            x2i=307E3-(t-lt)*vii*cosd(a2);
            y2i=(t-lt)*vii*sind(a2)-.5*9.8*(t-lt)^2;
            x3i=307E3-(t-lt)*vii*cosd(a3);
            y3i=(t-lt)*vii*sind(a3)-.5*9.8*(t-lt)^2;
            
        end
        
        function re = ending1(app,xm,ym)
            
            plot(app.UIAxes,xm,ym,'rx',xm,ym,'r+','MarkerSize' ,30);
            
            pause(1)
            
            set(app.TextArea_2,'Visible','on')
            
            pause(5)
            
            delete(app.UIFigure)
            
        end
        function re = ending2(app)
            
            plot(app.UIAxes,307E3,0,'rx',307E3,0,'r+','MarkerSize' ,30);
            
            pause(1)
            
            set(app.TextArea_2,'FontColor','r')
            set(app.TextArea_2,'FontSize',25)
            set(app.TextArea_2,'Value',[newline,newline,newline,newline,newline,newline,newline,newline,'Fail'])
            set(app.TextArea_2,'Visible','on')
            
            pause(5)
            
            delete(app.UIFigure)
            
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: BeginButton
        function BeginButtonPushed(app, event)
                        
            
%             [X1,map1] = imread('Untitled.png');
            
            vim=1735; %m/s
            la=45; %deg
            cityloc=307E3; %m
            lt=60; %seconds
            vii=2*vim;
            endd=0;
            
            t=0;
            
            x1i=-10;
            y1i=-10;
            x2i=-10;
            y2i=-10;
            x3i=-10;
            y3i=-10;
            
            a1=0; %28.643089147
            a2=0;
            a3=0;
            
            set(app.EnemyLaunchAngle45Label,'Visible','on')
            set(app.EnemyVelocity1735msLabel,'Visible','on')
            set(app.BeginButton,'Visible','off')
            
            while (endd==0)           
                
                tic;
                
                xm=t*vim*cosd(la);
                
                ym=t*vim*sind(la)-.5*9.8*t^2;
                
                if(mod(round(t,2),.1)==0)
                    
                     plot(app.UIAxes,xm,ym,'o',0,0,'d',cityloc,0,'p',x1i,y1i,'*',x2i,y2i,'*',x3i,y3i,'*','LineWidth',3);
                
%                      drawnow limitrate
                
                     axis(app.UIAxes,[-1,3.5E5,-1,1E5]);
                
                     set(app.Timesincelaunch0secondsLabel,'Text', ['Time since launch: ',num2str(round(t,2)),' seconds'])
                     
                     set(app.EnemyVelocity1735msLabel,'Text', ['Enemy Velocity: ',num2str(round(norm([vim*cosd(la),vim*sind(la)-9.8*t]),1)),' m/s'])
                     
                     if(y1i>=0)
                     
                         set(app.Shot1Velocity0msLabel,'Text', ['Shot 1 Velocity: ',num2str(round(norm([vii*cosd(a1),vii*sind(a1)-9.8*(t-lt)]),1)),' m/s'])
                     else
                         set(app.Shot1Velocity0msLabel,'Text', ['Shot 1 Velocity: ',num2str(0),' m/s'])
                     end
                     if(y2i>=0)
                         
                         set(app.Shot2Velocity0msLabel,'Text', ['Shot 2 Velocity: ',num2str(round(norm([vii*cosd(a2),vii*sind(a2)-9.8*(t-lt)]),1)),' m/s'])
                     else
                         set(app.Shot2Velocity0msLabel,'Text', ['Shot 2 Velocity: ',num2str(0),' m/s'])
                     end
                     if(y3i>=0)
                         
                         set(app.Shot3Velocity0msLabel,'Text', ['Shot 3 Velocity: ',num2str(round(norm([vii*cosd(a3),vii*sind(a3)-9.8*(t-lt)]),1)),' m/s'])
                     else
                         set(app.Shot3Velocity0msLabel,'Text', ['Shot 3 Velocity: ',num2str(0),' m/s'])
                     end
                
                end
                
                if(t>=60)
                
                    [x1i,y1i,x2i,y2i,x3i,y3i]=launch(app,lt,vii,a1,a2,a3,t);
                
                    if(collide(app,xm,ym,x1i,y1i,x2i,y2i,x3i,y3i,1)==1)
                        
                        endd=1; 
                        ending1(app,xm,ym);
                
                    end
                
                else
                    a1=app.Angle1Input.Value;
                    a2=app.Angle2Input.Value;
                    a3=app.Angle3Input.Value;
                end
                
                if(ym<=-.01)
                    
                    endd=1;
                    ending2(app);
                    
                end
                
                pause(.0001)
                
                d=collide(app,xm,ym,x1i,y1i,x2i,y2i,x3i,y3i,0);
                
                if(d(1)<=1000 || d(2)<=1000 || d(3)<=1000)%t>114.5 && t<114.8)
                
                   ta=.001;
                   t=t+ta;
                    
                else
                   
                   ta=toc;
                   t=t+ta;
                
                end
                
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {480, 480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {220, '1x', 220};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 860 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {220, '1x', 220};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create BeginButton
            app.BeginButton = uibutton(app.LeftPanel, 'push');
            app.BeginButton.ButtonPushedFcn = createCallbackFcn(app, @BeginButtonPushed, true);
            app.BeginButton.BackgroundColor = [1 0.6 0.851];
            app.BeginButton.FontSize = 30;
            app.BeginButton.Position = [32 60 156 113];
            app.BeginButton.Text = {'Begin'; ''};

            % Create Label
            app.Label = uilabel(app.LeftPanel);
            app.Label.HorizontalAlignment = 'center';
            app.Label.WordWrap = 'on';
            app.Label.Position = [16 196 188 269];
            app.Label.Text = 'The enemy has fired a missile! We need your help to protect our city. We are launching 3 interseptor missiles as soon as we can; we have found that we can launch 60 seconds after the enemy missile was launched. Can you find the right angle to shoot the missile down before we are destroyed? The interseptor missiles are twice as fast as our enemy''s at 3470 m/s. We also know the enemy missile started with a velocity of 1735 m/s and a launch angle of 45°. The enemy launch site is 307 km away.';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create EnemyLaunchAngle45Label
            app.EnemyLaunchAngle45Label = uilabel(app.CenterPanel);
            app.EnemyLaunchAngle45Label.HorizontalAlignment = 'center';
            app.EnemyLaunchAngle45Label.Visible = 'off';
            app.EnemyLaunchAngle45Label.Position = [19 47 144 22];
            app.EnemyLaunchAngle45Label.Text = 'Enemy Launch Angle: 45°';

            % Create EnemyVelocity1735msLabel
            app.EnemyVelocity1735msLabel = uilabel(app.CenterPanel);
            app.EnemyVelocity1735msLabel.HorizontalAlignment = 'center';
            app.EnemyVelocity1735msLabel.Visible = 'off';
            app.EnemyVelocity1735msLabel.Position = [10 26 161 22];
            app.EnemyVelocity1735msLabel.Text = 'Enemy Velocity: 1,735 m/s';

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.CenterPanel);
            app.TextArea_2.HorizontalAlignment = 'center';
            app.TextArea_2.FontSize = 21;
            app.TextArea_2.FontColor = [0 1 0];
            app.TextArea_2.BackgroundColor = [0 0 0];
            app.TextArea_2.Visible = 'off';
            app.TextArea_2.Position = [0 3 417 474];
            app.TextArea_2.Value = {''; ''; ''; ''; ''; ''; ''; ''; ''; 'Success'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.CenterPanel);
            title(app.UIAxes, 'Radar View')
            xlabel(app.UIAxes, 'Distance')
            ylabel(app.UIAxes, 'Height')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.PlotBoxAspectRatio = [1.20792079207921 1 1];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ColorOrder = [1 0 0;1 0.4118 0.1608;0 1 1;0 1 0;0 1 0;0 1 0;0.6353 0.0784 0.1843];
            app.UIAxes.Position = [3 62 413 357];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Create Shot1AngleLabel
            app.Shot1AngleLabel = uilabel(app.RightPanel);
            app.Shot1AngleLabel.HorizontalAlignment = 'right';
            app.Shot1AngleLabel.Position = [15 376 74 22];
            app.Shot1AngleLabel.Text = 'Shot 1 Angle';

            % Create Angle1Input
            app.Angle1Input = uieditfield(app.RightPanel, 'numeric');
            app.Angle1Input.Limits = [0 180];
            app.Angle1Input.ValueDisplayFormat = '%12.9g';
            app.Angle1Input.BackgroundColor = [1 0.6 0.851];
            app.Angle1Input.Position = [104 376 100 22];

            % Create Shot2AngleLabel
            app.Shot2AngleLabel = uilabel(app.RightPanel);
            app.Shot2AngleLabel.HorizontalAlignment = 'right';
            app.Shot2AngleLabel.Position = [15 287 74 22];
            app.Shot2AngleLabel.Text = 'Shot 2 Angle';

            % Create Angle2Input
            app.Angle2Input = uieditfield(app.RightPanel, 'numeric');
            app.Angle2Input.Limits = [0 180];
            app.Angle2Input.ValueDisplayFormat = '%12.9g';
            app.Angle2Input.BackgroundColor = [1 0.6 0.851];
            app.Angle2Input.Position = [104 287 100 22];

            % Create Shot3AngleLabel
            app.Shot3AngleLabel = uilabel(app.RightPanel);
            app.Shot3AngleLabel.HorizontalAlignment = 'right';
            app.Shot3AngleLabel.Position = [15 196 74 22];
            app.Shot3AngleLabel.Text = 'Shot 3 Angle';

            % Create Angle3Input
            app.Angle3Input = uieditfield(app.RightPanel, 'numeric');
            app.Angle3Input.Limits = [0 180];
            app.Angle3Input.ValueDisplayFormat = '%12.9g';
            app.Angle3Input.BackgroundColor = [1 0.6 0.851];
            app.Angle3Input.Position = [104 196 100 22];

            % Create Timesincelaunch0secondsLabel
            app.Timesincelaunch0secondsLabel = uilabel(app.RightPanel);
            app.Timesincelaunch0secondsLabel.HorizontalAlignment = 'center';
            app.Timesincelaunch0secondsLabel.FontSize = 14;
            app.Timesincelaunch0secondsLabel.Position = [-32 68 283 59];
            app.Timesincelaunch0secondsLabel.Text = 'Time since launch: 0 seconds';

            % Create Shot1Velocity0msLabel
            app.Shot1Velocity0msLabel = uilabel(app.RightPanel);
            app.Shot1Velocity0msLabel.HorizontalAlignment = 'center';
            app.Shot1Velocity0msLabel.Position = [15 342 189 22];
            app.Shot1Velocity0msLabel.Text = 'Shot 1 Velocity: 0 m/s';

            % Create Shot2Velocity0msLabel
            app.Shot2Velocity0msLabel = uilabel(app.RightPanel);
            app.Shot2Velocity0msLabel.HorizontalAlignment = 'center';
            app.Shot2Velocity0msLabel.Position = [15 251 189 22];
            app.Shot2Velocity0msLabel.Text = 'Shot 2 Velocity: 0 m/s';

            % Create Shot3Velocity0msLabel
            app.Shot3Velocity0msLabel = uilabel(app.RightPanel);
            app.Shot3Velocity0msLabel.HorizontalAlignment = 'center';
            app.Shot3Velocity0msLabel.Position = [15 160 189 22];
            app.Shot3Velocity0msLabel.Text = 'Shot 3 Velocity: 0 m/s';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end