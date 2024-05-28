function [images,ns,nmat] = isoimages(n1,n2,grinprofile,shapeparameters,discretization,pixelsize)

%  this script calculates the discretized RIs,
%  as well as the corresponding images

%% parameters
% n1 is the highest refractive index
% n2 is the lowest refractive index
% pixelsize is the pixel size of the printer in um

sps = shapeparameters;
nn = discretization;
ps = pixelsize/1000;

switch grinprofile
    case 'spherical'
        r = sps(1)/2;
        ro = sps(2)/2;
        
        gr = sqrt((n1-n2)*2/n1)/r;
        gr2 = -gr*r/(ro-r);

        np = ceil(ro*2/ps);

        y = (1:np)-0.5;
        x = (1:np)-0.5;
        
        [X,Y] = meshgrid(x,y);

        c = [np np]/2;

        Y = Y-c(1);
        X = X-c(2);
        
        rmat = sqrt(Y.^2+X.^2)*ps;

        nmat = zeros(np,np);

        grin = rmat <= r;
        shell = rmat > r & rmat <= ro;

        nmat(grin) = n1*(1-rmat(grin).^2*gr*gr/2);
        nmat(shell) = n1*(1-(ro-rmat(shell)).^2*gr2*gr2/2);

        rsgrin = (0:nn)*r/nn; rsgrin(1) = -1;
        rsshell = ro-(0:nn)*(ro-r)/nn;

        images = zeros(np,np,nn);
        ns = zeros(nn,1);

        for i = 1:nn
            nsi = nmat(rmat>rsgrin(i)&rmat<=rsgrin(i+1));
            ns(i) = mean(nsi(:));

            pixelsi = (rmat>rsgrin(i)&rmat<=rsgrin(i+1))...
                    |(rmat<=rsshell(i)&rmat>rsshell(i+1));
            
            m = zeros(size(nmat));
            m(pixelsi) = 1;
            images(:,:,i) = m;
        end
    
    case 'aspheric negative'
    case 'cylindrical'
        l = sps(1);
        w = sps(2);
        lo = sps(3);
        wo = sps(4);

        gr = sqrt((n1-n2)*2/n1)/(w/2);
        gr2 = -gr*w/(wo-w);

        np1 = ceil(wo/ps);
        np2 = ceil(l/ps);

        y = (1:np1)-0.5;
        x = (1:np2)-0.5;

        [~,Y] = meshgrid(x,y);

        Y = Y-np1/2;

        ymat = abs(Y)*ps;

        grin = ymat<=w/2;
        shell = ymat>w/2 & ymat<=wo/2;

        nmat = zeros(np1,np2);

        nmat(grin) = n1*(1-ymat(grin).^2*gr*gr/2);
        nmat(shell) = n1*(1-(wo/2-ymat(shell)).^2*gr2*gr2/2);

        ysgrin = (0:nn)*w/2/nn;
        ysshell = wo/2-(0:nn)*(wo-w)/2/nn;

        np3 = ceil(lo/ps);
        d1 = floor((np3-np2)/2);
        d2 = np3-np2-d1;
        
        images = zeros(np1,np3,nn);
        ns = zeros(nn,1);

        for i = 1:nn
            nsi = nmat(ymat>ysgrin(i)&ymat<=ysgrin(i+1));
            ns(i) = mean(nsi(:));
            
            pixelsi = (ymat>ysgrin(i)&ymat<=ysgrin(i+1))...
                    |(ymat<=ysshell(i)&ymat>ysshell(i+1));

            m = zeros(size(nmat));
            m(pixelsi) = 1;
            
            if isequal(i,1)
                m1 = [ones(np1,d1) m ones(np1,d2)];
            else
                m1 = [zeros(np1,d1) m zeros(np1,d2)];
            end
            images(:,:,i) = m1;
        end

    case 'sinusoidal'
    case 'tent'
end

% hf = figure;
% for i = 1:10:nn
%     imshow(uint8(images(:,:,i)*255));
%     pause;
% end
% close(hf);
