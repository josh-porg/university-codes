load('combustorData150000_169999.mat');
imshow(mat2grey(p(1:100:end,1:100:end)')); grid off; title 'Data Snapshots'; 
%yticks
%set(gca,'YTick',)
dmd = DynamicModeDecomposer(p',"rank",-1);
figure;
imshow(mat2grey(real(dmd.reconstructedSnapshots(1:100:end,1:100:end)))); grid off; title 'reconstructed Snapshots';

visulization

%%%%%%%% use the profiler!!!!!!!!!!!!!!

% plot time dynamics of each mode
f= figure;

plot(real(dmd.timeDynamics)')
f.WindowState='maximized';




% plot eigen values vs frequency
figure
plot(dmd.continuousTimeEigenValues,'+r')

%%
% plot singul
figure;

%ar values
semilogy(dmd.singularValues./sum(dmd.singularValues),'r+'); title('log 10(Singular Values)');
hold on;
semilogy(dmd.truncatedSingularValues/sum(dmd.singularValues),'bo');

figure;
plot(dmd.singularValues/sum(dmd.singularValues),'r+'); title('Singular Values');
hold on
plot(dmd.truncatedSingularValues/sum(dmd.singularValues),'bo');

figure
semilogy(dmd.singularValues./sum(dmd.singularValues),'o'); title('log 10(Singular Values)');

%fft2();
