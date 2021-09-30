#include "deform_tool.h"
#include <iostream>


void info(const deform_tool& dt)
  {
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "DEFORM TOOL" << std::endl;
  std::cout << "Triangles: " << dt.triangles.size() << std::endl;
  std::cout << "Vertices: " << dt.vertices.size() << std::endl;
  std::cout << "Coordinate system: " << std::endl;
  for (int i = 0; i < 4; ++i)
    {
    for (int j = 0; j < 4; ++j)
      {
      std::cout << dt.cs[i + 4 * j] << " ";
      }
    std::cout << std::endl;
    } 
  std::cout << "Visible: " << (dt.visible ? "Yes" : "No") << std::endl;
  std::cout << "Decay factor: " << dt.decay_factor << std::endl;
  std::cout << "Signed distance: " << (dt.signed_distance ? "Yes" : "No") << std::endl;
  std::cout << "Discretization: " << dt.discretization << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  }