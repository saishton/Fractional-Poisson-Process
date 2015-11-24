function [] = create_avi(data,startTime)

mapfilename = [startTime,'/create_avi-map-vid.avi'];
linksfilename = [startTime,'/create_avi-links-vid.avi'];
degdistfilename = [startTime,'/create_avi-degdist-vid.avi'];

contact_time = 20;
close all
scale = 20;

all_people = [data(:,2)' data(:,3)'];
num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max(all_people);
coords = zeros(num_people,2);
theta = 2*pi/num_people;

parfor n=1:num_people
    coords(n,:) = scale*[sin(n*theta) cos(n*theta)];
end

clustering = zeros(1,num_times);
map(num_times) = struct('cdata',[],'colormap',[]);
links(num_times) = struct('cdata',[],'colormap',[]);
degdist(num_times) = struct('cdata',[],'colormap',[]);
for m=1:num_times
    adj = zeros(num_people);
    current_time = (m-1)*contact_time;
    
    %== Create adjacency matrix and extract clustering coefficient ==%
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            adj(person1,person2) = 1;
            adj(person2,person1) = 1;
        end
    end
    adj2 = adj^2;
    adj3 = adj^3;
    adj2sum = sum(sum(adj2));
    contrip = adj2sum - trace(adj2);
    if contrip==0
        clustering(m) = 0;
    else
        clustering(m) = trace(adj3)/(adj2sum - trace(adj2));
    end
       
    %==Create Video Frames==%
    map_fig = figure();
    gplot(adj,coords,'-*');
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    map(m) = getframe(map_fig);
    close(map_fig);
    
%LINKS

    degdist_fig = figure();
    edges = -0.5:1:5.5;
    histogram(sum(adj),edges,'Normalization','cdf');
    str = sprintf('Time: %d', current_time);
    text(0.1,0.1,str);
    axis([-0.5 5.5 0.75 1]);
    ax = gca;
    ax.XTick = [0 1 2 3 4 5];
    xlabel('Degree');
    ylabel('Density');
    degdist(m) = getframe(degdist_fig);
    close(degdist_fig);
end

v = VideoWriter(mapfilename);
open(v)
writeVideo(v,map)
close(v)

% v = VideoWriter(linksfilename);
% open(v)
% writeVideo(v,links)
% close(v)

v = VideoWriter(degdistfilename);
open(v)
writeVideo(v,degdist)
close(v)

%==Plot clustering data==%
maxTime = (num_times-1)*contact_time;
T = linspace(0,maxTime,num_times);

hwin = 50;
Tmod = T;
Tmod((num_times+1-hwin):num_times) = [];
Tmod(1:hwin) = [];

MA = zeros(1,num_times-2*hwin);
parfor i=1:(num_times-2*hwin)
    upper = i+2*hwin;
    MA(i) = sum(clustering(i:upper))/(2*hwin+1);
end

figure()
hold on
plot(T,clustering)
plot(Tmod,MA,'LineWidth',4)
xlabel('Time (s)');
ylabel('Clustering Coefficient');
hold off

imagefilename = [startTime,'/create_avi-GCC-img.png'];
print(imagefilename,'-dpng')
close