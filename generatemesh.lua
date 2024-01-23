print("Mesh dimension:")
local dim = io.read()
dim = tonumber(dim)
if (dim > 2) then
  print("only applied for dim <= 2.")
  dim = 2
end
local path = "./"..dim.."D.mesh"
path = tostring(path)

local XL, XR, YL, YR
YL = 0; YR = 0; ZL = 0; ZR = 0
local numele, nvertex, numnode, numboundary
local nx, ny
nx = 0; ny = 0;

-- Get domain
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
  end
else goto label
end

-- Get mesh size
print("number of mesh in X")
nx = io.read()
nx = tonumber(nx)
if (dim >= 2) then
print("number of mesh in Y")
ny = io.read()
ny = tonumber(ny)
end

if (dim == 1) then
  numele = nx
elseif (dim == 2) then
  numele = nx * ny
end

local numnodex = nx + 1;
numnode = numnodex
local numnodey = 0
if (dim >= 2) then
numnodey = ny + 1;
numnode = numnodex * numnodey
end

-- Periodic
::periodic::print("If periodic?[y/n]")
local x_periodic, y_periodic
local peri = io.read()
peri = tostring(peri)
if (peri == "y") then
  peri = true
  if(dim == 1) then
  else
    ::xp:: print("x periodic?[y/n]")
    x_periodic = io.read()
    x_periodic = tostring(x_periodic)
    if (x_periodic == 'y') then
      x_periodic = true
    elseif (x_periodic == 'n') then
      x_periodic = false
    else goto xp
    end
    ::yp:: print("y periodic?[y/n]")
    y_periodic = io.read()
    y_periodic = tostring(y_periodic)
    if (y_periodic == 'y') then
      y_periodic = true
    elseif (y_periodic == 'n') then
      y_periodic = false
    else goto yp
    end
  end
elseif (peri == "n") then
  peri = false
else goto periodic
end

local nv_x, nv_y
if (dim == 1) then
  if (peri) then
    nvertex = nx
  else
    nvertex = nx + 1
  end
end
if (dim == 2) then
if (x_periodic) then
  nv_y = ny
else
  nv_y = ny + 1
end
if (y_periodic) then
  nv_x = nx
else
  nv_x = nx + 1
end
nvertex = nv_x * nv_y
end

-- coordinates
local coor = {}
numnode = 0
if (dim == 1) then
  for ix=0,numnodex-1 do
    local x = XL + (XR - XL) * ix / (numnodex - 1)
    coor[numnode+1] = {x}
    numnode = numnode + 1
  end
elseif (dim == 2) then
  for ix=0,numnodex-1 do
  for iy=0,numnodey-1 do
    local x = XL + (XR - XL) * ix / (numnodex - 1)
    local y = YL + (YR - YL) * iy / (numnodey - 1)
    coor[numnode+1] = {x,y}
    numnode = numnode + 1
  end
  end
end

-- vertex matrix
local vmt = {}
for i = 1, nx+1 do
  if (dim == 1) then
    vmt[i] = i-1
    if (peri and i == nx+1)then
      vmt[i] = 0
    end
  elseif (dim == 2) then
    vmt[i] = {}
    for j = 1, ny+1 do
      if (y_periodic and i == nx+1) then
        vmt[i][j] = vmt[1][j]
      elseif (x_periodic and j == ny+1) then
        vmt[i][j] = vmt[i][1]
      else
        vmt[i][j] = nv_y * (i-1) + j - 1
      end
    end
  end
end

-- node matrix, use this to define the mapping between vertices and nodes
local nodemt = {}
for i = 1, nx+1 do
  if (dim == 1) then
    nodemt[i] = i-1
  elseif (dim == 2) then
    nodemt[i] = {}
    for j = 1, ny+1 do
      nodemt[i][j] = (ny + 1) * (i-1) + j-1
    end
  end
end


-- element vertex
local el_v = {}
local el_n = {}
for i = 1, nx do
  if (dim == 1) then
    el_v[i] = {vmt[i],vmt[i+1]}
    el_n[i] = {nodemt[i], nodemt[i+1]}
  elseif (dim == 2) then
    el_v[i] = {}
    el_n[i] = {}
    for j = 1, ny do
      el_v[i][j] = {vmt[i][j],vmt[i+1][j], vmt[i+1][j+1], vmt[i][j+1]}
      el_n[i][j] = {nodemt[i][j],nodemt[i+1][j], nodemt[i+1][j+1], nodemt[i][j+1]}
    end
  end
end

-- Boundary
numboundary = 0
if (dim == 1) then
  if (not peri) then
    numboundary = 2
  end
elseif (dim == 2) then
  if (not x_periodic) then
    numboundary = numboundary + 2 * nx
  end
  if (not y_periodic) then
    numboundary = numboundary + 2 * ny
  end
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
  for ie = 1, nx do
      file:write (eattr," ",etype," ",el_v[ie][1]," ",el_v[ie][2],"\n")
  end
elseif (dim == 2) then
  for ie = 1,nx do
    for je = 1, ny do
      file:write (eattr," ",etype," ",el_v[ie][je][1]," ",el_v[ie][je][2]," ",el_v[ie][je][3]," ",el_v[ie][je][4],"\n")
    end
  end
end

file:write("\n")
file:write("boundary\n")

file:write(numboundary, "\n")
if (dim == 1) then
  if (not peri) then
    file:write(1," ",0," ",0,"\n")
    file:write(1," ",0," ",nx,"\n")
  end
elseif (dim == 2) then
  if (not x_periodic) then
    for ie = 1, nx do
      file:write(1," ",1," ",el_v[ie][1][1]," ",el_v[ie][1][2],"\n")
    end
    for ie = 1, nx do
      file:write(1," ",1," ",el_v[ie][ny][3]," ",el_v[ie][ny][4],"\n")
    end
  end
  if (not y_periodic) then
    for je = 1, ny do
      file:write(2," ",1," ",el_v[1][je][4]," ",el_v[1][je][1],"\n")
    end
    for je = 1, ny do
      file:write(2," ",1," ",el_v[nx][je][2]," ",el_v[nx][je][3],"\n")
    end
  end
end

file:write("\n")
file:write("vertices\n",numnode,"\n")
if (not peri) then
  file:write(dim,"\n")
  for ip = 0, numnode-1 do
    if (dim == 1) then
      file:write(coor[ip+1][1],"\n")
    elseif (dim == 2) then
      file:write(coor[ip+1][1]," ",coor[ip+1][2],"\n")
    end
  end
else
  file:write("nodes\n")
  if (dim == 1) then
    file:write("FiniteElementSpace\n","FiniteElementCollection: L2_T1_1D_P1\n")
  elseif(dim == 2) then
    file:write("FiniteElementSpace\n","FiniteElementCollection: L2_T1_2D_P1\n")
  end
  file:write("VDim: ", dim, "\n")
  file:write("Ordering: 1\n")

  if (dim == 1) then
    for ie = 1, nx do
      file:write("\n",coor[el_n[ie][1]+1][1],' ',coor[el_n[ie][2]+1][1])
    end
  elseif (dim == 2) then
    for ie = 1, nx do
      for je = 1, ny do
        file:write("\n")
        file:write("\n",coor[el_n[ie][je][1]+1][1], ' ', coor[el_n[ie][je][1]+1][2])
        file:write("\n",coor[el_n[ie][je][2]+1][1], ' ', coor[el_n[ie][je][2]+1][2])
        file:write("\n",coor[el_n[ie][je][4]+1][1], ' ', coor[el_n[ie][je][4]+1][2])
        file:write("\n",coor[el_n[ie][je][3]+1][1], ' ', coor[el_n[ie][je][3]+1][2])
      end
    end
  end

end

file:close()