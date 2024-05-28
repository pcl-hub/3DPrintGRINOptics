function fullimage = placetheimages(fullimagesize,image,qty,spacings,offcenter,pixelsize)

% spacings is a 1x2 variable, unit: mm
% where the 1st term is the separation in Dimension 1 (y) direction
% where the 2nd term is the separation in Dimension 2 (x) direction

if isempty(offcenter); offcenter = [0 0]; end % unit: mm
if isempty(pixelsize); pixelsize = 64.8; end % unit: um
if isempty(fullimagesize); fullimagesize = [1080 1920]; end

fullimage = zeros(fullimagesize);
cim = floor(fullimagesize/2);
spacings = ceil(spacings/pixelsize*1000);
offcenter = ceil(offcenter/pixelsize*1000);

dy = ceil(spacings(1)/2);
dx = ceil(spacings(2)/2);
y0 = cim(1)-ceil(size(image,1)/2);
x0 = cim(2)-ceil(size(image,2)/2);
yrange = (y0+1:y0+size(image,1))+offcenter(1);
xrange = (x0+1:x0+size(image,2))+offcenter(2);


switch qty
    case 1
        fullimage(yrange,xrange) = image;

    case 2
        yrange1 = yrange-dy;
        xrange1 = xrange-dx;
        yrange2 = yrange+dy;
        xrange2 = xrange+dx;

        fullimage(yrange1,xrange1) = image;
        fullimage(yrange2,xrange2) = image;

    case 3
        yrange1 = yrange-dy;
        xrange1 = xrange;
        yrange2 = yrange+dy;
        xrange2 = xrange-dx;
        yrange3 = yrange+dy;
        xrange3 = xrange+dx;

        fullimage(yrange1,xrange1) = image;
        fullimage(yrange2,xrange2) = image;
        fullimage(yrange3,xrange3) = image;
        
    case 4
        yrange1 = yrange-dy;
        xrange1 = xrange-dx;
        yrange2 = yrange+dy;
        xrange2 = xrange-dx;
        yrange3 = yrange-dy;
        xrange3 = xrange+dx;
        yrange4 = yrange+dy;
        xrange4 = xrange+dx;
        
        fullimage(yrange1,xrange1) = image;
        fullimage(yrange2,xrange2) = image;
        fullimage(yrange3,xrange3) = image;
        fullimage(yrange4,xrange4) = image;
end




