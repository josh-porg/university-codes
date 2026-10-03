
clear
clc
close all


name=input('Hello! What is your name? ','s');

disp(['Hi, ', name, ', it is nice to meet you! What is your favorite color?'])

color=input('','s');

if(color(1)=='g')
   reponse=', and it is actually my favorite!';
else
    reponse=', but my favorite is green!';
end

disp(['I really like ',color, reponse])

sibcount=input('How many sibilings do you have? ');

if(sibcount>=2)
   amt='a lot!';
else
    amt='not that many.';
end

disp(['Wow! ',num2str(sibcount),' is ',amt,' I have 6 siblings!'])
