for i = -1:1:1
    for j = -1:1:1
        for k = -1:1:1
            [t,x] = getBlackboxData2D(i,j,k);
            figure;
            plot3(x(:,1), x(:,2), x(:,3), 'w', 'LineWidth', 1.5);
            set(gca, 'color','k', 'xcolor','w', 'ycolor','w' , 'zcolor','w');
            set(gcf, 'color','k');
        end
    end
end

[t,x] = getBlackboxData2D(0, 1, 20);
plot3(x(:,1), x(:,2), x(:,3), 'w', 'LineWidth', 1.5);
set(gca, 'color','k', 'xcolor','w', 'ycolor','w' , 'zcolor','w');
set(gcf, 'color','k');