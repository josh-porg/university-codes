% tarrayfun = squeeze(zeros(10000,1));
% tcomparison = squeeze(zeros(10000,1));

for i = 1:1000
    tic;
    temp0 = trace(S);
    tarrayfun(i) = toc;
    tic;
    temp1 = sum(S,'all');
    tcomparison(i) = toc;
end

figure;
subplot(2,1,1);
histogram(tarrayfun),title("tarrayfun");
subplot(2,1,2);
histogram(tcomparison),title("tcomparison");

sv = tarrayfun(tarrayfun > 1);

figure;
subplot(2,1,1);
histogram(sv),title("svd");
subplot(2,1,2);
histogram(rsv),title("rsvd");