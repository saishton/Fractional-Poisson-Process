function [] = create_avi(data,videofilename)

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
    coords(n,:) = 20*[sin(n*theta) cos(n*theta)];
end

M(num_times) = struct('cdata',[],'colormap',[]);
for m=1:num_times
    adj = zeros(num_people);
    current_time = (m-1)*contact_time;
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            adj(person1,person2) = 1;
            adj(person2,person1) = 1;
        end
    end
    this_fig = figure();
    gplot(adj,coords,'-*');
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    M(m) = getframe(this_fig);
    close(this_fig);
end

v = VideoWriter(videofilename);
open(v)
writeVideo(v,M)
close(v)

