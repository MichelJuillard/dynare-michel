function []=graph_decomp(z,varlist,initial_period,freq)
    global M_ 
    exo_nbr = M_.exo_nbr;

    gend = size(z,3);
    x = initial_period-1/freq:(1/freq):initial_period+(gend-1)/freq;
    
    [i_var,nvar] = varlist_indices(varlist);
    for j=1:nvar;
        z1 = squeeze(z(i_var(j),:,:));
        xmin = x(1);
        xmax = x(end);
        ix = z1 > 0;
        ymax = max(sum(z1.*ix));
        ix = z1 < 0;
        ymin = min(sum(z1.*ix));
        if ymax-ymin < 1e-6
            continue
        end
        figure('Name',M_.endo_names(i_var(j),:));
        ax=axes('Position',[0.1 0.1 0.6 0.8]);
        axis(ax,[xmin xmax ymin ymax]);
        plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
        hold on;
        for i=1:gend
            i_1 = i-1;
            yp = 0;
            ym = 0;
            for k = 1:exo_nbr+1 
                zz = z1(k,i);
                if zz > 0
                    fill([x(i) x(i) x(i+1) x(i+1)],[yp yp+zz yp+zz yp],k);
                    yp = yp+zz;
                else
                    fill([x(i) x(i) x(i+1) x(i+1)],[ym ym+zz ym+zz ym],k);
                    ym = ym+zz;
                end
                hold on;
            end
        end
        plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
        hold off;

        axes('Position',[0.75 0.1 0.2 0.8]);
        axis([0 1 0 1]);
        axis off;
        hold on;
        y1 = 0;
        height = 1/(exo_nbr+1);
        labels = strvcat(M_.exo_names,'Initial values');
        
        for j=1:exo_nbr+1
            fill([0 0 0.2 0.2],[y1 y1+0.7*height y1+0.7*height y1],j);
            hold on
            text(0.3,y1+0.3*height,labels(j,:),'Interpreter','none');
            hold on
            y1 = y1 + height;
        end
        hold off
    end