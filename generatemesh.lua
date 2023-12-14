print("Mesh dimension:")
local dim = io.read()
dim = tonumber(dim)
if (dim > 2) then
  print("only applied for dim <= 2.")
  dim = 2
end
local path = "./"..dim.."D.mesh"
path = tostring(path)

local XL, XR, YL, YR, ZL, ZR
YL = 0; YR = 0; ZL = 0; ZR = 0
local numele, npoint
local nx, ny, nz
nx = 0; ny = 0; nz = 0;

::label::print("uniform mesh [y/n]?")
local ifuni = io.read()
ifuni = tostring(ifuni)
if (ifuni == "y") then
  XL = 0; XR = 1;
  YL = 0; YR = 1;
elseif (ifuni == "n") then
  print("X range:")
  XL = io.read()
  XR = io.read()
  XL = tonumber(XL)
  XR = tonumber(XR)
  if (dim>=2) then
    print("Y range:")
    YL = io.read()
    YR = io.read()
    YL = tonumber(YL)
    YR = tonumber(YR)
    if (dim == 3) then
      print("Z range:")
      ZL = io.read()
      XR = io.read()
      ZL = tonumber(ZL)
      ZR = tonumber(ZR)
    end
  end
else goto label
end

print("number of mesh in X")
nx = io.read()
nx = tonumber(nx)
if (dim >= 2) then
print("number of mesh in Y")
ny = io.read()
ny = tonumber(ny)
  if (dim == 3) then
    print("number of mesh in Z")
    nz = io.read()
    nz = tonumber(nz)
  end
end

::periodic::print("If periodic?[y/n]")
local ifperi = io.read()
ifperi = tostring(ifperi)
if (ifperi == "y") then
  dom = {}
  dom.xl = XL
  dom.xr = XR
  dom.yl = YL
  dom.yr = YR
  periodic_square(dim, dom, nx, ny, nz)
  dom = nil
  os.exit()
elseif (ifperi == "n") then
else goto periodic
end


local npx = nx + 1;
local npy = ny + 1;

-- coordinates
local coor = {}
npoint = 0
if (dim == 1) then
  for ix = 0, nx do
    local x = XL + (XR - XL) * ix / nx
    coor[npoint+1] = {x}
    npoint = npoint+1
  end
elseif (dim == 2) then
  for ix=0,nx do
  for iy=0,ny do
    local x = XL + (XR - XL) * ix / nx
    local y = YL + (YR - YL) * iy / ny
    coor[npoint+1] = {x,y}
    npoint = npoint+1
  end
  end
end

if (dim == 1) then
  numele = nx
elseif (dim == 2) then
  numele = nx * ny
end

-- write mesh
local file = io.open(path,"w+")
if (file == nil) then
  print("file dosen't exist")
  os.exit()
end
-- io.output(file)
file:write([[
MFEM mesh v1.0

#
# MFEM Geometry Types (see mesh/geom.hpp):
#
# POINT       = 0
# SEGMENT     = 1
# TRIANGLE    = 2
# SQUARE      = 3
# TETRAHEDRON = 4
# CUBE        = 5
#

dimension
]],
dim,
"\n\n"
)

file:write("elements\n",numele,"\n")

-- write element vertices
local eattr = 1
local etype = 3
if (dim == 1) then
  etype = 1
end

if (dim == 1) then
  for ie = 0, nx-1 do
    local v1 = ie
    local v2 = v1 + 1
      file:write (eattr," ",etype," ",v1," ",v2,"\n")
  end
elseif (dim == 2) then
  for ie = 0,nx-1 do
    for je = 0, ny-1 do
      local v1 = npx*ie + je
      local v2 = v1 + npx
      local v3 = v2 + 1
      local v4 = v1 + 1
      file:write (eattr," ",etype," ",v1," ",v2," ",v3," ",v4,"\n")
    end
  end
end

file:write("\n")
file:write("boundary\n")

local nb = 2 * (nx + ny)
if (dim == 1) then
  nb = 2
end

file:write(nb, "\n")
if (dim == 1) then
  file:write(1," ",0," ",0,"\n")
  file:write(1," ",0," ",nx,"\n")
elseif (dim == 2) then
  for ie = 0, nx-1 do
    file:write(1," ",1," ",npx*ie," ",npx*(ie+1),"\n")
  end
  for ie = 0, nx-1 do
    file:write(1," ",1," ",npx*ie+npy-1," ",npx*(ie+1)+npy-1,"\n")
  end
  for je = 0, ny-1 do
    file:write(2," ",1," ",je," ",je+1,"\n")
  end
  for je = 0, ny-1 do
    file:write(2," ",1," ",(npx-1)*npx+je," ",(npx-1)*npx+je+1,"\n")
  end
end
file:write("\n")
file:write("vertices\n",npoint,"\n")
file:write(dim,"\n")
for ip = 0, npoint-1 do
  if (dim == 1) then
    file:write(coor[ip+1][1],"\n")
  elseif (dim == 2) then
    file:write(coor[ip+1][1]," ",coor[ip+1][2],"\n")
  end
end

file:close()

function periodic_square(dim, dom, nx, ny, nz)
  -- write periodic_square.mesh
  if (dim ~= 2) then
    print("need dim == 2");
    os.exit()
  end

  local vertices, nodes
  os.exit()
end
