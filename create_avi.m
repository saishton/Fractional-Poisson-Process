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

node_freq = zeros(num_people,1);
line_freq = zeros(num_people);
clustering = zeros(1,num_times);
numlinks = zeros(1,num_times);
map(num_times) = struct('cdata',[],'colormap',[]);
links(num_times) = struct('cdata',[],'colormap',[]);
degdist(num_times) = struct('cdata',[],'colormap',[]);
ML_vid(num_times) = struct('Parameters',[],'Statistics',[]);
GP_vid(num_times) = struct('Parameters',[],'Statistics',[]);
WB_vid(num_times) = struct('Parameters',[],'Statistics',[]);
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
    node_freq = node_freq+sum(adj,2);
    rel_node_freq = node_freq/max(node_freq);
    line_freq = line_freq+adj;
    rel_line_freq = line_freq/(max(max(line_freq)));
    adjsum = sum(sum(adj));
    numlinks(m) = adjsum/2;
    
    %==Degree Distribution Fitting==%
    raw_data = sum(adj);
    [F,X] = ecdf(raw_data);
    ccdf_data = 1-F;
    
    fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
    ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
    [cf_ml,gof_ml] = fit(X,ccdf_data,ft_ml);
    cv_ml = coeffvalues(cf_ml);

    fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
    ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
    [cf_gp,gof_gp] = fit(X,ccdf_data,ft_gp);
    cv_gp = coeffvalues(cf_gp);

    fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
    ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
    [cf_wb,gof_wb] = fit(X,ccdf_data,ft_wb);
    cv_wb = coeffvalues(cf_wb);
    
    beta = cv_ml(1);
    gamma = cv_ml(2);
    ccdf_ml = mlf(beta,1,-gamma*X.^beta,6);

    k = cv_gp(1);
    sigma = cv_gp(2);
    theta = cv_gp(3);
    ccdf_gp = gpcdf(X,k,sigma,theta,'upper');

    a = cv_wb(1);
    b = cv_wb(2);
    ccdf_wb = wblcdf(X,a,b,'upper');
    
    %==GoF Testing==%
    test_data = sort(raw_data)';

    z_ml = ones(length(test_data),1)-mlf(beta,1,-gamma*test_data.^beta,6);
    z_gp = gpcdf(test_data,k,sigma,theta);
    z_wb = wblcdf(test_data,a,b);

    stats_ml = testStatistics(raw_data,z_ml);
    stats_gp = testStatistics(raw_data,z_gp);
    stats_wb = testStatistics(raw_data,z_wb);

    stats_ml.Chi_Squared = sum(rdivide((ccdf_data-ccdf_ml).^2,ccdf_ml));
    stats_gp.Chi_Squared = sum(rdivide((ccdf_data-ccdf_ml).^2,ccdf_gp));
    stats_wb.Chi_Squared = sum(rdivide((ccdf_data-ccdf_ml).^2,ccdf_wb));

    stats_ml.Root_MSE = gof_ml.rmse;
    stats_gp.Root_MSE = gof_gp.rmse;
    stats_wb.Root_MSE = gof_wb.rmse;
    stats_ml.R_Squared = gof_ml.rsquare;
    stats_gp.R_Squared = gof_gp.rsquare;
    stats_wb.R_Squared = gof_wb.rsquare;
    
    %==Save to Structure==%
    struc_ml = struct('Stability',beta,'Scale',gamma);
    struc_gp = struct('Shape',k,'Scale',sigma,'Location',theta);
    struc_wb = struct('Scale',a,'Shape',b);
    
    ML_vid(m).Parameters = struc_ml;
    GP_vid(m).Parameters = struc_gp;
    WB_vid(m).Parameters = struc_wb;
    ML_vid(m).Statistics = stats_ml;
    GP_vid(m).Statistics = stats_gp;
    WB_vid(m).Statistics = stats_wb;
       
    %==Create Video Frames==%
    map_fig = figure();
    gplot(adj,coords,'-*');
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    map(m) = getframe(map_fig);
    close(map_fig);
    
    links_fig = figure();
    link_size = rel_node_freq*50;
    tempadj = logical(rel_line_freq);
    [row,col] = find(tempadj);
    tempcoords = zeros(num_people,2);
    hold on
    parfor i=1:length(row)
        thisrow = row(i);
        thiscol = col(i);
        if thisrow >= thiscol
            line_col = (1-rel_line_freq(thisrow,thiscol))*[1 1 1];
            x1 = coords(thisrow,1);
            y1 = coords(thisrow,2);
            x2 = coords(thiscol,1);
            y2 = coords(thiscol,2);
            line([x1 x2],[y1 y2],'Color',line_col)
        end
        tempcoords(i,:) = coords(thisrow,:);
    end
    tempcoords = tempcoords(any(tempcoords,2),:);
    tempcoords = unique(tempcoords,'rows','stable');
    link_size(link_size==0) = [];
    scatter(tempcoords(:,1),tempcoords(:,2),link_size,'filled')
    hold off
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    links(m) = getframe(links_fig);
    close(links_fig);
    
    degdist_fig = figure();
    plot(X,ccdf_data,'o');
    plot(X,ccdf_ml);
    plot(X,ccdf_gp);
    plot(X,ccdf_wb);
    str = sprintf('Time: %d', current_time);
    text(0.1,0.1,str);
    axis([-1 6 0 1]);
    ax = gca;
    ax.XTick = [0 1 2 3 4 5];
    xlabel('Degree');
    ylabel('CCDF');
    legend('Data','ML','Gen. Pareto','Weibull');
    degdist(m) = getframe(degdist_fig);
    close(degdist_fig);
end

v = VideoWriter(mapfilename);
open(v)
writeVideo(v,map)
close(v)

v = VideoWriter(linksfilename);
open(v)
writeVideo(v,links)
close(v)

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

fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,gof_ex] = fit(X_links,ccdf_links,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
[cf_gm,gof_gm] = fit(X_links,ccdf_links,ft_gm);
cv_gm = coeffvalues(cf_gm);

fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
[cf_rl,gof_rl] = fit(X_links,ccdf_links,ft_rl);
cv_rl = coeffvalues(cf_rl);

%==BEST FIT CCDFs==%
lambda = cv_ex(1);
ccdf_ex = expcdf(X_links,lambda,'upper');

a = cv_gm(1);
b = cv_gm(2);
ccdf_gm = gamcdf(X,a,b,'upper');

sigma = cv_rl(1);
ccdf_rl = raylcdf(x,sigma,'upper');

%==GoF Testing==%
test_data = sort(numlinks)';

z_ex = expcdf(test_data,lambda);
z_gm = gamcdf(test_data,a,b);
z_rl = raylcdf(test_data,sigma);

stats_ex = testStatistics(numlinks,z_ex);
stats_gm = testStatistics(numlinks,z_gm);
stats_rl = testStatistics(numlinks,z_rl);

stats_ex.Root_MSE = gof_ex.rmse;
stats_gm.Root_MSE = gof_gm.rmse;
stats_rl.Root_MSE = gof_rl.rmse;
stats_ex.R_Squared = gof_ex.rsquare;
stats_gm.R_Squared = gof_gm.rsquare;
stats_rl.R_Squared = gof_rl.rsquare;

%==Plotting==%
linksactivefig = figure();
hold on
plot(X_links,ccdf_links,'o');
plot(X_links,ccdf_ex);
plot(X_links,ccdf_gm);
plot(X_links,ccdf_rl);
xlabel('Number of Active Links');
ylabel('CCDF');
legend('Data','Exp','Gamma','Rayleigh');
imagefilename = [dir_ref,'/create_avi-LAD-img.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Save Relevant Data==%
struc_ex = struct('Rate',lambda);
struc_gm = struct('Shape',a,'Scale',b);
struc_rl = struct('Scale',sigma);

EX = struct('Parameters',struc_ex,'Statistics',stats_ml);
GM = struct('Parameters',struc_gm,'Statistics',stats_gm);
RL = struct('Parameters',struc_rl,'Statistics',stats_rl);

datafilename = [dir_ref,'/create_avi-data.mat'];
save(datafilename,'ML_vid','GP_vid','WB_vid''EX','GM','RL')