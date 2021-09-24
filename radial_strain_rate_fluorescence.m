clear all;close all;clc

DIRECTORY = dir('*gfp*.tif');
theta=linspace(0,2*pi,101);
se=strel('diamond',3);

%% Find mask for initial frame and calculate fluorescence signal
fr=1;
a=imread(DIRECTORY(fr).name);
a=imgaussfilt(a,2);
gr=imgradient(a);
edges=gr>2.5*mean2(gr);
edges=imfill(edges,'holes');
labels=bwlabel(edges);
for i =1:max(max(labels))
    temp=regionprops(labels==i,'Area');
    Area(i)=temp.Area;
end
edges=labels==find(Area==max(max(Area)));   
mask_I=imfill(edges,'holes');

a=double(a);  
cell_mask=double(mask_I);
cell_mask(cell_mask==0)=nan;
bg_mask=double(mask_I==0);
bg_mask(bg_mask==0)=nan;
BACKGROUND_SIGNAL=nanmean(nanmean(bg_mask.*a));
FLUORESCENCE_SIGNAL=nanmean(nanmean(cell_mask.*(a-BACKGROUND_SIGNAL)))
    
    
%% Loop through frames and calculate radial distance as a function of angle and time
for fr = 1:size(DIRECTORY,1)
    fr
%% Find Mask
a=imread(DIRECTORY(fr).name);
edges=a>.95*mean2(a);
edges=imfill(edges,'holes');
edges=imerode(edges,se);
edges=imerode(edges,se);
labels=bwlabel(edges);
clear Area
for i =1:max(max(labels))
	temp=regionprops(labels==i,'Area');
	Area(i)=temp.Area;
end
mask=labels==find(Area==max(max(Area))); 
mask=imdilate(mask,se);
mask=imdilate(mask,se);
edges=imgradient(mask);
edges=edges>0;
    
%% Find x,y coordinates of boundary and convert to polar coordinates
x=[];y=[];
for i = 1:size(a,1)
	for j = 1:size(a,2)
        if edges(i,j)==1
            x=[x j];
            y=[y i];
        end
	end
end
clear angle r
[xc yc R] = circfit(x,y);  
for i = 1:length(x)
	r(i)=norm([x(i)-xc,y(i)-yc]);
	if y(i)>yc
        angle(i) = acos(dot( ([x(i)-xc,y(i)-yc])/norm(([x(i)-xc,y(i)-yc])), [ 1 0]));
    else
        angle(i) = 2*pi-acos(dot( ([x(i)-xc,y(i)-yc])/norm(([x(i)-xc,y(i)-yc])), [ 1 0]));
	end
end
param=polyfit(angle,r,50);
x=linspace(0,2*pi,500);
y1=polyval(param,x);
r_theta_p=polyder(param);
r_2theta_p=polyder(r_theta_p);
r_theta=polyval(r_theta_p,x);
r_2theta=polyval(r_2theta_p,x);
for i = (length(theta)-1)/2 - (length(theta)-1)/4:(length(theta)-1)/2 + (length(theta)-1)/4
	RADIUS(i) = mean (  ((y1(x>theta(i) & x<theta(i+1)).^2 + r_theta(x>theta(i) & x<theta(i+1)).^2).^(3/2)) ./ abs(y1(x>theta(i) & x<theta(i+1)).^2 +2*r_theta(x>theta(i) & x<theta(i+1)).^2 - y1(x>theta(i) & x<theta(i+1)).*r_2theta(x>theta(i) & x<theta(i+1))));
end
        
RR(fr,126:375)=y1(126:375);
    
   
[angle2,I]=sort(angle);
r2=r(I);
angle2=angle2+pi;
angle2(angle2>2*pi)=angle2(angle2>2*pi)-2*pi;        

param=polyfit(angle2,r2,50);
x=linspace(0,2*pi,500);
y1=polyval(param,x);
r_theta_p=polyder(param);
r_2theta_p=polyder(r_theta_p);
r_theta=polyval(r_theta_p,x);
r_2theta=polyval(r_2theta_p,x);
x=x-pi;
x(x<0)=2*pi+x(x<0);  
y1=circshift(y1,250);
RR(fr,1:125)=y1(1:125);
RR(fr,376:500)=y1(376:500);

for i = 1:(length(theta)-1)/2 - (length(theta)-1)/4
	RADIUS(i) = mean (  ((y1(x>theta(i) & x<theta(i+1)).^2 + r_theta(x>theta(i) & x<theta(i+1)).^2).^(3/2)) ./ abs(y1(x>theta(i) & x<theta(i+1)).^2 +2*r_theta(x>theta(i) & x<theta(i+1)).^2 - y1(x>theta(i) & x<theta(i+1)).*r_2theta(x>theta(i) & x<theta(i+1))));
end
for i = (length(theta)-1)/2 + (length(theta)-1)/4:length(theta)-1
	RADIUS(i) = mean (  ((y1(x>theta(i) & x<theta(i+1)).^2 + r_theta(x>theta(i) & x<theta(i+1)).^2).^(3/2)) ./ abs(y1(x>theta(i) & x<theta(i+1)).^2 +2*r_theta(x>theta(i) & x<theta(i+1)).^2 - y1(x>theta(i) & x<theta(i+1)).*r_2theta(x>theta(i) & x<theta(i+1))));
end  
for i = 1:(length(theta)-1)/2 - (length(theta)-1)/4
	RADIUS(i) = mean (  ((y1(x>theta(i) & x<theta(i+1)).^2 + r_theta(x>theta(i) & x<theta(i+1)).^2).^(3/2)) ./ abs(y1(x>theta(i) & x<theta(i+1)).^2 +2*r_theta(x>theta(i) & x<theta(i+1)).^2 - y1(x>theta(i) & x<theta(i+1)).*r_2theta(x>theta(i) & x<theta(i+1))));
end
    
RADIUS_OF_CURVATURE(fr,:) = RADIUS;
  
x=linspace(0,2*pi,500);
XX=RR(fr,:).*cos(x)+xc;
YY=RR(fr,:).*sin(x)+yc;
imagesc(a);
hold on
plot(XX,YY,'r','Linewidth',2)
drawnow;
G(fr)=getframe;
end

%% Correct for intrinsic oocyte curvature
for i = 1:size(RADIUS_OF_CURVATURE,1)
	RADIUS_OF_CURVATURE(i,:)=RADIUS_OF_CURVATURE(i,:)-RADIUS_OF_CURVATURE(1,:);
end
figure
imshow(RADIUS_OF_CURVATURE(1:end,:)',[0 500])
colorbar

save('Variable.mat')

%% Smooth outer contour and calculate characteristic radial strain rate
for i = 1:size(RR,2);RR(:,i)=smooth(RR(:,i));end
characteristic_radial_strain_rate=abs(6*mean(min(diff(RR./mean(RR)))))
