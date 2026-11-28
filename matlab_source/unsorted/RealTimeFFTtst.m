clear
close all
clc

sampfreq=100;
freq1=5;
freq2=10;
figure

for i=1:500

    x(i)=i/sampfreq;

    y(i)=sin(2*pi*freq2*x(i))+sin(2*pi*freq1*x(i))+.5*randn(1);

    f=fft(y);

    f=f.*conj(f)/length(x);

    freq=sampfreq/(length(x))*(0:length(x)-1);

    h=1:ceil(length(x)/2);
    subplot(1,3,1)
    plot(x,y)
    xlabel('Time')
    ylabel('Amplitude')
    subplot(1,3,2)
    scatter(freq(h),f(h))
    xlabel('Frequency')
    ylabel('Energy')
    subplot(1,3,3)
    scatter(log(freq(h)),f(h))
    xlabel('log(Frequency)')
    ylabel('Energy')

    drawnow
   
%     disp(num2str(i)+" "+num2str(freq(f(h)>mean(f(h))+5*std(f(h)))))

%     pause(1)

end