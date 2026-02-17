function Vq = interp2_wx(varargin)
%MATLAB Code Generation Library Function

%   Limitations and Notes:
%   1. Xq and Yq must be the same size. Use MESHGRID to evaluate on on
%      a grid.
%   2. For best performance supply X and Y as vectors.
%   3. The 'cubic' method requires a uniformly-spaced grid, otherwise an
%      error is issued. Use the 'spline' method explicitly for an
%      unevenly-spaced grid.
%   4. The 'spline' method is optimized for Xq and Yq being outputs of
%      MESHGRID and for a small number of interpolation points relative to
%      the dimensions of V. Using the 'spline' method to interpolate over a
%      large set of scattered points may be inefficient.

%   Copyright 1984-2014 The MathWorks, Inc.
%#codegen

narginchk(1,5);
coder.internal.prefer_const(varargin{:});
if nargin == 5

	% interp2(X,Y,V,Xq,Yq)
	X = varargin{1};
	Y = varargin{2};
	V = varargin{3};
	Xq = varargin{4};
	Yq = varargin{5};
	METHOD = LINEAR;
	EXTRAP = default_extrap(METHOD);
	extrapval = cast(coder.internal.nan,'like',V);
	interp2_validate(V,Xq,Yq,METHOD,extrapval,X,Y);
	Vq = interp2_dispatch(V,Xq,Yq,METHOD,EXTRAP,extrapval, ...
		unmeshgrid(X,2),unmeshgrid(Y,1));
end

%--------------------------------------------------------------------------

function interp2_validate(V,Xq,Yq,METHOD,extrapval,X,Y)
AUTOGRID = nargin <= 5;
eml_invariant(ismatrix(V), ...
    'Coder:toolbox:interp2_VMustBe2D');
eml_invariant(size(V,2) >= 2 && size(V,1) >= 2, ...
    'Coder:toolbox:NotEnoughPoints');
eml_invariant(isfloat(Xq) && isfloat(Yq), ...
    'Coder:toolbox:interp2_invalidXYClass');
eml_invariant(isreal(Xq) && isreal(Yq), ...
    'Coder:toolbox:interp2_nonrealXY');
eml_invariant(isa(Xq,class(Yq)), ...
    'Coder:toolbox:interp2_classXIandYIMustMatch');
eml_invariant(isequal(size(Xq),size(Yq)), ...
    'Coder:toolbox:interp2_sizeXIandYIMustMatch');
if ~AUTOGRID
    eml_invariant(isfloat(X) && isfloat(Y), ...
        'Coder:toolbox:interp2_invalidXYClass');
    eml_invariant(isreal(X) && isreal(Y), ...
        'Coder:toolbox:interp2_nonrealXY');
    eml_invariant(isa(X,class(Y)), ...
        'Coder:toolbox:interp2_mixedTypeGrid');
    if eml_is_const(isvector(X)) && isvector(X) && ...
            eml_is_const(isvector(Y)) && isvector(Y)
        nx = numel(X);
        ny = numel(Y);
    else
        eml_invariant(isplaid('meshgrid',X,Y), ...
            'Coder:toolbox:GridMustBeVectorsOrPlaid');
        nx = size(X,2);
        ny = size(Y,1);
    end
    eml_invariant(ismatrix(V) && ...
        size(V,1) == ny && ...
        size(V,2) == nx, ...
        'MATLAB:xyzchk:lengthXAndYDoNotMatchSizeZ');
end
if METHOD ~= NEAREST
    eml_invariant(isfloat(V), ...
        'Coder:toolbox:interp2_ZMustBeFloat');
    eml_invariant(coder.internal.isBuiltInNumeric(extrapval), ...
        'MATLAB:interp2:extrapvalNotNumeric');
end
% The extrapolation value cannot be complex unless z is complex.
eml_invariant(isreal(extrapval) || ~isreal(V), ...
    'Coder:toolbox:interp2_extrapvalNotReal');

%--------------------------------------------------------------------------

function Vq = interp2_ntimes(V,N,METHOD)
coder.internal.prefer_const(N,METHOD);
eml_invariant(ismatrix(V), ...
    'Coder:toolbox:interp2_VMustBe2D');
eml_invariant(isscalar(N) && isnumeric(N), ...
    'Coder:toolbox:NMustBeANumericScalar');
eml_invariant(ismatrix(V), ...
    'Coder:toolbox:interp2_VMustBe2D');
eml_invariant(size(V,2) >= 2 && size(V,1) >= 2, ...
    'Coder:toolbox:NotEnoughPoints');
eml_invariant(METHOD == NEAREST || isfloat(V), ...
    'Coder:toolbox:interp2_ZMustBeFloat');
ny = size(V,1);
nx = size(V,2);
if isa(V,'single')
    d = 1/pow2(single(floor(N)));
else
    d = 1/pow2(double(floor(N)));
end
[Xq,Yq] = meshgrid(1:d:nx,1:d:ny);
Vq = interp2_dispatch(V,Xq,Yq,METHOD,false,cast(0,'like',V));

%--------------------------------------------------------------------------

function p = is_unit_colon(x)
% Determine if x is a *constant* colon expression of the form 1:numel(x).
coder.internal.prefer_const(x);
if eml_is_const(x)
    p = coder.const(feval('isequal',x,feval('colon',1,numel(x))));
else
    p = false;
end

%--------------------------------------------------------------------------

function Vq = interp2_dispatch(V,Xq,Yq,METHOD,EXTRAP,extrapval,X,Y)
coder.internal.prefer_const(METHOD,EXTRAP);
if nargin <= 6
    AUTOGRID = true;
else
    coder.internal.prefer_const(X,Y);
    if coder.const(is_unit_colon(X) && is_unit_colon(Y))
        AUTOGRID = true;
    else
        AUTOGRID = false;
    end
end
if METHOD == SPLINE
    % Since SPLINE does not really support an automatic grid, allow the
    % X and Y inputs to pass through, even if they are of the form 1:n.
    if AUTOGRID && nargin <= 6
        Vq = interp2_spline(V,Xq,Yq,EXTRAP,extrapval);
    else
        Vq = interp2_spline(V,Xq,Yq,EXTRAP,extrapval,X,Y);
    end
elseif METHOD == CUBIC
    if AUTOGRID
        Vq = interp2_cubic(V,Xq,Yq,EXTRAP,extrapval);
    else
        Vq = interp2_cubic(V,Xq,Yq,EXTRAP,extrapval,X,Y);
    end
else
    if AUTOGRID
        Vq = interp2_local(METHOD,V,Xq,Yq,EXTRAP,extrapval);
    else
        Vq = interp2_local(METHOD,V,Xq,Yq,EXTRAP,extrapval,X,Y);
    end
end

%--------------------------------------------------------------------------

function Vq = interp2_local(method,V,Xq,Yq,EXTRAP,extrapval,X,Y)
% The local interpolation methods:  linear and nearest.
% This function assumes numel(x)>=2, numel(y)>=2, and numel(z)>=2 and that
% the sizes are consistent with size(v).
coder.internal.prefer_const(method,EXTRAP);
AUTOGRID = nargin <= 6;
if method == LINEAR
    interpf = @scalar_bilinear_interp;
else
    interpf = @scalar_nearest_interp;
end
Vq = coder.nullcopy(zeros(size(Xq),'like',V));
n = coder.internal.indexInt(numel(Vq));
if coder.internal.useParforConst('interp2_orders_0_1',n)
    if AUTOGRID
        ixmax = coder.internal.indexInt(size(V,2) - 1);
        iymax = coder.internal.indexInt(size(V,1) - 1);
        parfor k = 1:n
            Vq(k) = local_autogrid_loop_body(interpf,Xq(k),Yq(k),V, ...
                ixmax,iymax,EXTRAP,extrapval);
        end
    else
        parfor k = 1:n
            Vq(k) = local_loop_body(interpf,Xq(k),Yq(k),X,Y,V, ...
                EXTRAP,extrapval);
        end
    end
else
    if AUTOGRID
        ixmax = coder.internal.indexInt(size(V,2) - 1);
        iymax = coder.internal.indexInt(size(V,1) - 1);
        for k = 1:n
            Vq(k) = local_autogrid_loop_body(interpf,Xq(k),Yq(k),V, ...
                ixmax,iymax,EXTRAP,extrapval);
        end
    else
        for k = 1:n
            Vq(k) = local_loop_body(interpf,Xq(k),Yq(k),X,Y,V, ...
                EXTRAP,extrapval);
        end
    end
end

%--------------------------------------------------------------------------

function Vqk = local_loop_body(interpf,Xqk,Yqk,X,Y,V,EXTRAP,extrapval)
coder.inline('always');
coder.internal.prefer_const(interpf,EXTRAP);
if EXTRAP || ...
        Xqk >= X(1) && Xqk <= X(end) && ...
        Yqk >= Y(1) && Yqk <= Y(end)
    ix = coder.internal.bsearch(X,Xqk);
    iy = coder.internal.bsearch(Y,Yqk);
    x1 = X(ix);
    x2 = X(ix + 1);
    y1 = Y(iy);
    y2 = Y(iy + 1);
    zx1y1 = V(iy,ix);
    zx2y1 = V(iy,ix + 1);
    zx1y2 = V(iy + 1,ix);
    zx2y2 = V(iy + 1,ix + 1);
    Vqk = interpf(x1,x2,y1,y2,zx1y1,zx2y1,zx1y2,zx2y2,Xqk,Yqk);
else
    Vqk = extrapval;
end

%--------------------------------------------------------------------------

function Vqk = local_autogrid_loop_body(interpf, ...
    Xqk,Yqk,V,ixmax,iymax,EXTRAP,extrapval)
coder.inline('always');
coder.internal.prefer_const(interpf,ixmax,iymax,EXTRAP);
if EXTRAP || ...
        Xqk >= 1 && Xqk <= size(V,2) && ...
        Yqk >= 1 && Yqk <= size(V,1)
    ix = autogrid_lookup(Xqk,ixmax);
    iy = autogrid_lookup(Yqk,iymax);
    x1 = cast(ix,'like',Xqk);
    x2 = x1 + 1;
    y1 = cast(iy,'like',Yqk);
    y2 = y1 + 1;
    zx1y1 = V(iy,ix);
    zx2y1 = V(iy,ix + 1);
    zx1y2 = V(iy + 1,ix);
    zx2y2 = V(iy + 1,ix + 1);
    Vqk = interpf(x1,x2,y1,y2,zx1y1,zx2y1,zx1y2,zx2y2,Xqk,Yqk);
else
    Vqk = extrapval;
end

%--------------------------------------------------------------------------

function idx = autogrid_lookup(x,maxidx)
% Substitute for bsearch when the grid is automatic.
coder.inline('always')
if x <= 1
    idx = coder.internal.indexInt(1);
elseif x <= maxidx
    idx = coder.internal.indexInt(floor(x));
else
    idx = coder.internal.indexInt(maxidx);
end

%--------------------------------------------------------------------------

function zi = scalar_bilinear_interp( ...
    x1,x2,y1,y2, ...
    zx1y1,zx2y1,zx1y2,zx2y2, ...
    xi,yi)
coder.inline('always');
T = zeros('like',zx1y1);
% zi = (1 - ry)*(onemrx*zx1y1 + rx*zx2y1) + ry*(onemrx*zx1y2 + rx*zx2y2);
if xi == x1
    qx1 = cast(zx1y1,'like',T);
    qx2 = cast(zx1y2,'like',T);
elseif xi == x2
    qx1 = cast(zx2y1,'like',T);
    qx2 = cast(zx2y2,'like',T);
else
    rx = (xi - x1)/(x2 - x1);
    onemrx = 1 - rx;
    if zx1y1 == zx2y1
        qx1 = cast(zx1y1,'like',T);
    else
        qx1 = cast(onemrx*zx1y1 + rx*zx2y1,'like',T);
    end
    if zx1y2 == zx2y2
        qx2 = cast(zx1y2,'like',T);
    else
        qx2 = cast(onemrx*zx1y2 + rx*zx2y2,'like',T);
    end
end
if yi == y1 || qx1 == qx2
    zi = qx1;
elseif yi == y2
    zi = qx2;
else
    ry = (yi - y1)/(y2 - y1);
    zi = cast((1 - ry)*qx1 + ry*qx2,'like',T);
end

%--------------------------------------------------------------------------

function zi = scalar_nearest_interp( ...
    x1,x2,y1,y2, ...
    zx1y1,zx2y1,zx1y2,zx2y2, ...
    xi,yi)
coder.inline('always');
dx1 = xi - x1;
dx2 = x2 - xi;
dy1 = yi - y1;
dy2 = y2 - yi;
if dx1 < dx2
    if dy1 < dy2
        zi = zx1y1;
    else
        zi = zx1y2;
    end
else
    if dy1 < dy2
        zi = zx2y1;
    else
        zi = zx2y2;
    end
end

%--------------------------------------------------------------------------

function Vq = interp2_cubic(V,Xq,Yq,EXTRAP,extrapval,X,Y)
coder.inline('always');
T = zeros('like',V);
RT = zeros('like',real(V));
AUTOGRID = nargin <= 5;
if AUTOGRID
    nx = coder.internal.indexInt(size(V,2));
    ny = coder.internal.indexInt(size(V,1));
    ixmax = nx - 1;
    iymax = ny - 1;
else
    dx = cast(EqualSpacing(X),'like',RT);
    dy = cast(EqualSpacing(Y),'like',RT);
    eml_invariant(dx > 0 && dy > 0, ...
        'Coder:toolbox:CubicGridMustBeUniform');
    nx = coder.internal.indexInt(numel(X));
    ny = coder.internal.indexInt(numel(Y));
    xmin = X(1);
    ymin = Y(1);
    xmax = X(nx);
    ymax = Y(ny);
end
eml_invariant(ny >= 3 && nx >= 3, ...
    'MATLAB:interp2:cubic:Vsize')
Vq = coder.nullcopy(zeros(size(Xq),'like',T));
% Expand V so interpolation is valid at the boundaries.
VV = zeros([ny+2,nx+2],'like',T);
for ix = 1:nx
    for iy = 1:ny
        VV(iy+1,ix+1) = V(iy,ix);
    end
end
for ix = 1:nx+2
    VV(1,ix) = 3*VV(2,ix) - 3*VV(3,ix) + VV(4,ix); % Y edges
    VV(ny+2,ix)  = 3*VV(ny+1,ix) - 3*VV(ny,ix) + VV(ny-1,ix);
end
for iy = 1:ny+2
    VV(iy,1) = 3*VV(iy,2) - 3*VV(iy,3) + VV(iy,4); % X edges
    VV(iy,nx+2)  = 3*VV(iy,nx+1) - 3*VV(iy,nx) + VV(iy,nx-1);
end
nout = coder.internal.indexInt(numel(Xq));
if coder.internal.useParforConst('interp2_orders_2_3',nout)
    if AUTOGRID
        parfor k = 1:nout
            Vq(k) = cubic_autogrid_loop_body(Xq(k),Yq(k),VV, ...
                ixmax,iymax,EXTRAP,extrapval);
        end
    else
        parfor k = 1:nout
            Vq(k) = cubic_loop_body(Xq(k),Yq(k),X,Y,VV, ...
                xmin,xmax,ymin,ymax,dx,dy,EXTRAP,extrapval);
        end
    end
else
    if AUTOGRID
        for k = 1:nout
            Vq(k) = cubic_autogrid_loop_body(Xq(k),Yq(k),VV, ...
                ixmax,iymax,EXTRAP,extrapval);
        end
    else
        for k = 1:nout
            Vq(k) = cubic_loop_body(Xq(k),Yq(k),X,Y,VV, ...
                xmin,xmax,ymin,ymax,dx,dy,EXTRAP,extrapval);
        end
    end
end

function X = localevaluate(x,iter)
coder.inline('always');
coder.internal.prefer_const(iter);
if iter == 0
    X = ((2-x).*x-1).*x;
elseif iter == 1
    X = (3*x-5).*x.*x+2;
elseif iter == 2
    X = ((4-3*x).*x+1).*x;
else
    X = (x-1).*x.*x;
end

%--------------------------------------------------------------------------

function Vqk = cubic_loop_body(Xqk,Yqk,X,Y,VV, ...
    xmin,xmax,ymin,ymax,dx,dy,EXTRAP,extrapval)
coder.inline('always');
coder.internal.prefer_const(EXTRAP);
realtype = zeros('like',real(VV));
if EXTRAP || ( ...
        Xqk >= xmin && Xqk <= xmax && ...
        Yqk >= ymin && Yqk <= ymax)
    ix = coder.internal.bsearch(X,Xqk);
    iy = coder.internal.bsearch(Y,Yqk);
    s = (cast(Xqk,'like',realtype) - cast(X(ix),'like',realtype))/dx;
    t = (cast(Yqk,'like',realtype) - cast(Y(iy),'like',realtype))/dy;
    zik = zeros('like',VV);
    for is = coder.unroll(0:3)
        ss = localevaluate(s,is);
        for it = coder.unroll(0:3)
            tt = localevaluate(t,it);
            zik = zik + VV(iy+it,ix+is)*ss*tt;
        end
    end
    Vqk = zik/4;
else
    Vqk = extrapval;
end

%--------------------------------------------------------------------------

function Vqk = cubic_autogrid_loop_body(Xqk,Yqk,VV, ...
    ixmax,iymax,EXTRAP,extrapval)
coder.inline('always');
coder.internal.prefer_const(EXTRAP);
realtype = zeros('like',real(VV));
if EXTRAP || ( ...
        Xqk >= 1 && Xqk <= (ixmax + 1) && ...
        Yqk >= 1 && Yqk <= (iymax + 1))
    ix = autogrid_lookup(Xqk,ixmax);
    iy = autogrid_lookup(Yqk,iymax);
    s = cast(Xqk,'like',realtype) - cast(ix,'like',realtype);
    t = cast(Yqk,'like',realtype) - cast(iy,'like',realtype);
    zik = zeros('like',VV);
    for is = coder.unroll(0:3)
        ss = localevaluate(s,is);
        for it = coder.unroll(0:3)
            tt = localevaluate(t,it);
            zik = zik + VV(iy+it,ix+is)*ss*tt;
        end
    end
    Vqk = zik/4;
else
    Vqk = extrapval;
end

%--------------------------------------------------------------------------

function Vq = interp2_spline(V,Xq,Yq,EXTRAP,extrapval,X,Y)
%2-D spline interpolation
if nargin < 7
    X = ones('like',real(V)):size(V,2);
    Y = ones('like',real(V)):size(V,1);
end
Vq = TensorInterp23(@spline,Y,X,V,Yq,Xq);
if ~EXTRAP
    Vq = MaskOutOfRangeValues(extrapval,Vq,X,Y,Xq,Yq);
end

%--------------------------------------------------------------------------

function m = LINEAR
coder.inline('always');
m = uint8(0);

function m = NEAREST
coder.inline('always');
m = uint8(1);

function m = CUBIC
coder.inline('always');
m = uint8(3);

function m = SPLINE
coder.inline('always');
m = uint8(4);

%--------------------------------------------------------------------------

function m = StringToMethodID(method)
coder.internal.prefer_const(method);
eml_invariant(eml_is_const(method), ...
    'Coder:toolbox:MethodMustBeConstant');
n = coder.const(max(1,length(method)));
islinear = coder.const(strncmpi(method,'linear',n));
isnearest = coder.const(strncmpi(method,'nearest',n));
iscubic = coder.const(strncmpi(method,'cubic',n));
isspline = coder.const(strncmpi(method,'spline',n));
eml_invariant(islinear || isnearest || iscubic || isspline, ...
    'Coder:toolbox:interp2_unsupportedMethod');
if isnearest
    m = NEAREST;
elseif iscubic
    m = CUBIC;
elseif isspline
    m = SPLINE;
else
    m = LINEAR;
end

%--------------------------------------------------------------------------

function p = default_extrap(METHOD)
if METHOD == SPLINE
    p = true;
else
    p = false;
end

%--------------------------------------------------------------------------
