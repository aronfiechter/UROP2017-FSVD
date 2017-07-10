----------------------------------------------------------------------
-- CGAL test ipelet description
----------------------------------------------------------------------

label = "Test"

about = [[
This ipelet is based on the CGAL_ipelet package. See www.cgal.org.
]]

-- this variable will store the C++ ipelet when it has been loaded
ipelet = false

function run(model, num)
  if not ipelet then ipelet = assert(ipe.Ipelet(dllname)) end
  model:runIpelet(methods[num].label, ipelet, num)
end

methods = {
  { label="two points euclidean bisector" },
  { label="two points L_inf bisector" },
  { label="point/segment L_inf-parabola" },
  { label="two sites L_inf bisector" },
  { label="Linf 2D Voronoi Diagram" },
  { label="L2 farthest Voronoi Diagram" },
  { label="Hausdorff Voronoi Diagram" },
  { label="L2 FVD polygonal input" },
  { label="farthest color Voronoi diagram" },
  { label="L2 NVD polygonal input" },
  { label="Help" },
}

----------------------------------------------------------------------
