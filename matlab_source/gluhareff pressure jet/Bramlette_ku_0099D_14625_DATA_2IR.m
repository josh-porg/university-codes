%Richard Bramlette
%University of Kansas
%Ph.D Research
%
%IMG2IR - Code to Read in an Image of an Object Heated to Incandescence and
%Convert the Information to Apparent Temperature Contours Scaled Using
%Measurements of the Actual Maximum and Minimum Temperatures

%User Input Values
Tmax = 2000;                                %Max Temperature in deg F
Tmin = 800;                                 %Min Temperature in deg F
dTContour = 10;                             %Contour Step Size, deg F

t_start = cputime;
%Read in User-Selected Image File
[filename pathname filterindex] = uigetfile({'*.png;*.jpg;*.tif;*.gif'});
[image map alpha] = imread(filename);       %Read Image at filename
[X Y N] = size(image);                      %Get Image Size
fprintf('%s was read into memory (%i Pixels Tall, %i Pixels Wide)\n',...
    filename,X,Y')

%Convert Red/Green/Blue Image to Black/White
BW = rgb2gray(image);                      %Note: BW is uint8 data
Imax = single(max(BW(:)));                 %Check Max Brightness
Imin = single(min(BW(:)));                 %Check Min Brightness

%Scale BW Data to Max and Min Temperature
Tapp = single(BW)*((Tmax - Tmin)/(Imax - Imin)) + Tmin;

%Plot Scaled Image Data Contours
figure(3)

subplot(2,1,1)
imshow(filename)                            %Show Original Image
title(sprintf('Original Image: %s',filename))

fprintf('Generating Contours of Apparent Temperature From Image File...\n')

subplot(2,1,2)
v = [Tmin:dTContour:Tmax];                  %Array of Temp. Contour Values
[Cf,hf] = contourf(Tapp,v);                 %Show Temperature Contours
if length(v) > 10                           %If more than 10 Contours,
    hf.LineColor = 'none';                  %Turn Off Contour Lines
end

%Limit Plot to Image Size with no skew
axis([0 Y 0 X],'equal')

%Flip Image Y-Axis and Remove Irrelevant Pixel Numbers on Axes
set(gca,'ydir','reverse','XTickLabel',[],'YTickLabel',[])

%Add Contour Plot Title, Colorbar, and Colorbar Title
title(sprintf('Apparent Temperature Contours, \\DeltaT = %2.0f \\circF',...
    dTContour))
Cb = colorbar('southoutside');
Cb.Label.String = 'Apparent Temperature (\circF)';

t_end = cputime;
t_run = t_end - t_start;
fprintf('Contour Image Processing Complete in %0.1f sec\n',t_run)