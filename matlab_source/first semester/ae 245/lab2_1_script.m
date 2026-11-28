% figure
% 
% for i=1:3
%    subplot(3,2,1+(i-1)*2)
%    plot(GPS(:,2),GPS(:,7+i))
% end
% hold on
% i = desired colum + (desired row - 1)*(number of rows)
% for i=1:3
%    subplot(3,2,2+(i-1)*2)
%    plot(ATT(:,2),ATT(:,2+2*i))
% end
figure
plot3(GPS(:,8), GPS(:,9), GPS(:,10))