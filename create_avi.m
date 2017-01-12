function [] = create_avi(data,dir_ref)

mapfilename = [dir_ref,'/create_avi-map-vid.avi'];
linksfilename = [dir_ref,'/create_avi-links-vid.avi'];
degdistfilename = [dir_ref,'/create_avi-degdist-vid.avi'];

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

step_nf = zeros(num_people,1);
line_freq = zeros(num_people,num_people);
clustering = zeros(1,num_times);
numlinks = zeros(1,num_times);
raw_data = zeros(num_times,num_people);

% links(num_times) = struct('cdata',[],'colormap',[]);
% map(num_times) = struct('cdata',[],'colormap',[]);

for m=1:num_times
    thisadj = zeros(num_people);
    current_time = (m-1)*contact_time;
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            thisadj(person1,person2) = 1;
            thisadj(person2,person1) = 1;
        end
    end
    adj2 = thisadj^2;
    adj3 = thisadj^3;
    adj2sum = sum(sum(adj2));
    contrip = adj2sum - trace(adj2);
    if contrip==0
        clustering(m) = 0;
    else
        clustering(m) = trace(adj3)/contrip;
    end
    step_nf = step_nf+sum(thisadj,2);
    adjsum = sum(sum(thisadj));
    numlinks(m) = adjsum/2;
    this_rd = sum(thisadj);
    raw_data(m,:) = this_rd;
    
%     %==Create Adj Frame==%
%     map_fig = figure();
%     gplot(thisadj,coords,'-*');
%     str = sprintf('Time: %d', current_time);
%     text(0,-1.2*scale,str);
%     axis([-1.5 1.5 -1.5 1.5]*scale);
%     axis off;
%     set(gcf,'color','w');
%     map(m) = getframe(map_fig);
%     close(map_fig);
%     %==End Frame Creation==%
    
%     this_rel_node_freq = step_nf/max(step_nf);
% 	line_freq = line_freq+thisadj;
%     thisRLF = line_freq/(max(max(line_freq)));
    
%     %==Create Activity Frame==%
% 	links_fig = figure();
%     link_size = this_rel_node_freq*50;
%     tempadj = logical(thisRLF);
%     [row,col] = find(tempadj);
%     tempcoords = zeros(num_people,2);
%     hold on
%     for i=1:length(row)
%         thisrow = row(i);
%         thiscol = col(i);
%         if thisrow >= thiscol
%             line_col = (1-thisRLF(thisrow,thiscol))*[1 1 1];
%             x1 = coords(thisrow,1);
%             y1 = coords(thisrow,2);
%             x2 = coords(thiscol,1);
%             y2 = coords(thiscol,2);
%             line([x1 x2],[y1 y2],'Color',line_col)
%         end
%         tempcoords(i,:) = coords(thisrow,:);
%     end
%     tempcoords = tempcoords(any(tempcoords,2),:);
%     tempcoords = unique(tempcoords,'rows','stable');
%     link_size(link_size==0) = [];
%     scatter(tempcoords(:,1),tempcoords(:,2),link_size,'filled')
%     hold off
%     str = sprintf('Time: %d', current_time);
%     text(0,-1.2*scale,str);
%     axis([-1.5 1.5 -1.5 1.5]*scale);
%     axis off;
%     set(gcf,'color','w');
%     links(m) = getframe(links_fig);
%     close(links_fig);
end

raw_data( :, ~any(raw_data,1) ) = [];
highX = max(max(raw_data))+2;
X = zeros(highX,num_times);
ccdf_data = zeros(highX,num_times);
amountPad = zeros(1,num_times);

fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);

lambda = zeros(1,num_times);
a = zeros(1,num_times);
b = zeros(1,num_times);
sigma = zeros(1,num_times);

% EX_vid(num_times) = struct('Parameters',[],'Statistics',[]);
% GM_vid(num_times) = struct('Parameters',[],'Statistics',[]);
% RL_vid(num_times) = struct('Parameters',[],'Statistics',[]);

parfor m=1:num_times
%     current_data = raw_data(m,:);
%     [thisF,thisX] = ecdf(current_data);
%     thisCCDF = 1-thisF;
    
%     %==Degree Distribution Fitting==%
%     [this_cf_ex,this_gof_ex] = fit(thisX,thisCCDF,ft_ex);
%     this_cv_ex = coeffvalues(this_cf_ex);
%     [this_cf_gm,this_gof_gm] = fit(thisX,thisCCDF,ft_gm);
%     this_cv_gm = coeffvalues(this_cf_gm);
%     [this_cf_rl,this_gof_rl] = fit(thisX,thisCCDF,ft_rl);
%     this_cv_rl = coeffvalues(this_cf_rl);
%     this_lambda = this_cv_ex(1);
%     this_a = this_cv_gm(1);
%     this_b = this_cv_gm(2);
%     this_sigma = this_cv_rl(1);
    
%     %==GoF Tests==%
%     this_test_data = sort(current_data)';
%    
%     z_ex = expcdf(this_test_data,this_lambda);
%     z_gm = gamcdf(this_test_data,this_a,this_b);
%     z_rl = raylcdf(this_test_data,this_sigma);
% 
%     stats_ex = testStatistics(this_test_data,z_ex);
%     stats_gm = testStatistics(this_test_data,z_gm);
%     stats_rl = testStatistics(this_test_data,z_rl);
% 
%     stats_ex.Root_MSE = this_gof_ex.rmse;
%     stats_gm.Root_MSE = this_gof_gm.rmse;
%     stats_rl.Root_MSE = this_gof_rl.rmse;
%     stats_ex.R_Squared = this_gof_ex.rsquare;
%     stats_gm.R_Squared = this_gof_gm.rsquare;
%     stats_rl.R_Squared = this_gof_rl.rsquare;
    
%     %==Save to Structure==%
%     struc_ex = struct('Rate',this_lambda);
%     struc_gm = struct('Shape',this_a,'Scale',this_b);
%     struc_rl = struct('Scale',this_sigma);
%     
%     EX_vid(m).Parameters = struc_ex;
%     GM_vid(m).Parameters = struc_gm;
%     RL_vid(m).Parameters = struc_rl;
%     EX_vid(m).Statistics = stats_ex;
%     GM_vid(m).Statistics = stats_gm;
%     RL_vid(m).Statistics = stats_rl;
    
%     amountPad(m) = highX - length(thisX);
%     thisCCDF_p = zeros(highX,1);
%     thisCCDF_p(1:length(thisCCDF),1) = thisCCDF;
%     thisX_p = zeros(highX,1);
%     thisX_p(1:length(thisX),1) = thisX;
%     X(:,m) = thisX_p;
%     ccdf_data(:,m) = thisCCDF_p;
%     lambda(m) = this_lambda;
%     a(m) = this_a;
%     b(m) = this_b;
%     sigma(m) = this_sigma;
end

% degdist(num_times) = struct('cdata',[],'colormap',[]);

for m=1:num_times
%     current_time = (m-1)*contact_time;
%     padding = amountPad(m);
%     thisX_p = X(:,m);
%     thisCCDF_p = ccdf_data(:,m);
%     thisX = thisX_p(1:end-padding);
%     thisCCDF = thisCCDF_p(1:end-padding);
    
%     %==BEST FIT CCDFs==%
%     ccdf_ex = expcdf(thisX,lambda(m),'upper');
%     ccdf_gm = gamcdf(thisX,a(m),b(m),'upper');
%     ccdf_rl = raylcdf(thisX,sigma(m),'upper');
       
%     %==Create Video Frames==%
%     degdist_fig = figure();
%     hold on
%     plot(thisX,thisCCDF,'o')
%     plot(thisX,ccdf_ex)
%     plot(thisX,ccdf_gm)
%     plot(thisX,ccdf_rl)
%     str = sprintf('Time: %d', current_time);
%     text(0.1,0.1,str);
%     axis([-1 6 0 1]);
%     ax = gca;
%     ax.XTick = [0 1 2 3 4 5];
%     xlabel('Degree');
%     ylabel('CCDF');
%     legend('Data','Exp','Gamma','Rayleigh');
%     hold off
%     degdist(m) = getframe(degdist_fig);
%     close(degdist_fig);
end

% v = VideoWriter(mapfilename);
% open(v)
% writeVideo(v,map)
% close(v)
% 
% v = VideoWriter(linksfilename);
% open(v)
% writeVideo(v,links)
% close(v)
% 
% v = VideoWriter(degdistfilename);
% open(v)
% writeVideo(v,degdist)
% close(v)

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

clusteringfig = figure();
hold on
plot(T,clustering)
plot(Tmod,MA,'LineWidth',4)
xlabel('Time (s)');
ylabel('Clustering Coefficient');
hold off
imagefilename = [dir_ref,'/create_avi-GCC-img.png'];
print(imagefilename,'-dpng')
close(clusteringfig);

%==LINKS ACTIVE FITTING==%
[F_links,X_links] = ecdf(numlinks);
ccdf_links = 1-F_links;

links_fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
links_ft_ex = fittype('expcdf(x,lambda,''upper'')','options',links_fo_ex);
[links_cf_ex,links_gof_ex] = fit(X_links,ccdf_links,links_ft_ex);
links_cv_ex = coeffvalues(links_cf_ex);

links_fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
links_ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',links_fo_gm);
[links_cf_gm,links_gof_gm] = fit(X_links,ccdf_links,links_ft_gm);
links_cv_gm = coeffvalues(links_cf_gm);

links_fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
links_ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',links_fo_rl);
[links_cf_rl,links_gof_rl] = fit(X_links,ccdf_links,links_ft_rl);
links_cv_rl = coeffvalues(links_cf_rl);

%==BEST FIT CCDFs==%
links_lambda = links_cv_ex(1);
links_ccdf_ex = expcdf(X_links,links_lambda,'upper');

links_a = links_cv_gm(1);
links_b = links_cv_gm(2);
links_ccdf_gm = gamcdf(X_links,links_a,links_b,'upper');

links_sigma = links_cv_rl(1);
links_ccdf_rl = raylcdf(X_links,links_sigma,'upper');

%==GoF Testing==%
links_test_data = sort(numlinks)';

links_z_ex = expcdf(links_test_data,links_lambda);
links_z_gm = gamcdf(links_test_data,links_a,links_b);
links_z_rl = raylcdf(links_test_data,links_sigma);

links_stats_ex = testStatistics(links_test_data,links_z_ex);
links_stats_gm = testStatistics(links_test_data,links_z_gm);
links_stats_rl = testStatistics(links_test_data,links_z_rl);

links_stats_ex.Root_MSE = links_gof_ex.rmse;
links_stats_gm.Root_MSE = links_gof_gm.rmse;
links_stats_rl.Root_MSE = links_gof_rl.rmse;
links_stats_ex.R_Squared = links_gof_ex.rsquare;
links_stats_gm.R_Squared = links_gof_gm.rsquare;
links_stats_rl.R_Squared = links_gof_rl.rsquare;

%==Plotting==%
linksactivefig = figure();
hold on
plot(X_links,ccdf_links,'o');
plot(X_links,links_ccdf_ex);
plot(X_links,links_ccdf_gm);
plot(X_links,links_ccdf_rl);
xlabel('Number of Active Links');
ylabel('CCDF');
legend('Data','Exp','Gamma','Rayleigh');
imagefilename = [dir_ref,'/create_avi-LAD-img.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Save Relevant Data==%
links_struc_ex = struct('Rate',links_lambda);
links_struc_gm = struct('Shape',links_a,'Scale',links_b);
links_struc_rl = struct('Scale',links_sigma);

EX = struct('Parameters',links_struc_ex,'Statistics',links_stats_ex);
GM = struct('Parameters',links_struc_gm,'Statistics',links_stats_gm);
RL = struct('Parameters',links_struc_rl,'Statistics',links_stats_rl);

datafilename = [dir_ref,'/create_avi-data.mat'];
%save(datafilename,'EX_vid','GM_vid','RL_vid','EX','GM','RL')
save(datafilename,'EX','GM','RL')